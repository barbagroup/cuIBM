#include "NSWithBody.h"
#include "kernels/calculateForce.h"

/// Initialise the bodies
template <typename memoryType>
void NSWithBody<memoryType>::initialiseBodies()
{
	NavierStokesSolver<memoryType>::logger.startTimer("initialiseBodies");
	
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	B.initialise(db, *NavierStokesSolver<memoryType>::domInfo);

	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/forces";
	int timeStep(db["simulation"]["startStep"].get<int>());
	if (timeStep != 0)
	{
		forceFile.open(out.str().c_str(), std::ofstream::app);
	}
	else
	{
		forceFile.open(out.str().c_str());
	}
	
	NavierStokesSolver<memoryType>::logger.stopTimer("initialiseBodies");
}
	
/// Update the body information at each time step during motion
template <typename memoryType>
void NSWithBody<memoryType>::updateBodies()
{
	NavierStokesSolver<memoryType>::logger.startTimer("updateBodies");
	
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	real dt   = db["simulation"]["dt"].get<real>();
	real Time = dt*(NavierStokesSolver<memoryType>::timeStep+1);
	B.update(db, *NavierStokesSolver<memoryType>::domInfo, Time);
	
	NavierStokesSolver<memoryType>::logger.stopTimer("updateBodies");
};

template <>
void NSWithBody<host_memory>::calculateForce()
{
}

template <>
void NSWithBody<device_memory>::calculateForce()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	real dt = (*paramDB)["simulation"]["dt"].get<real>(),
	     nu = (*paramDB)["flow"]["nu"].get<real>();
	
	real *q_r      = thrust::raw_pointer_cast(&q[0]),
		 *qOld_r   = thrust::raw_pointer_cast(&qOld[0]),
		 *lambda_r = thrust::raw_pointer_cast(&lambda[0]),
		 *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
		 *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	// Calculating drag
	cusp::array1d<real, device_memory>
		FxX(B.numCellsY[0]),
		FxY(B.numCellsX[0]+1),
		FxU((B.numCellsX[0]+1)*B.numCellsY[0]);
	
	real *FxX_r = thrust::raw_pointer_cast(&FxX[0]),
	     *FxY_r = thrust::raw_pointer_cast(&FxY[0]),
	     *FxU_r = thrust::raw_pointer_cast(&FxU[0]);

	const int blockSize = 256;
	dim3 dimGrid( int((B.numCellsX[0]+B.numCellsY[0]+1-0.5)/blockSize)+1, 1 );
	dim3 dimBlock(blockSize, 1);
	
	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxX_r, q_r, lambda_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxY_r, q_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridX( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
	kernels::dragUnsteady <<<dimGridX, dimBlock>>> (FxU_r, q_r, qOld_r, dxD, dyD, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	forceX = thrust::reduce(FxX.begin(), FxX.end()) + thrust::reduce(FxY.begin(), FxY.end()) + thrust::reduce(FxU.begin(), FxU.end());
	
	// Calculating lift
	cusp::array1d<real, device_memory>
		FyX(B.numCellsY[0]+1),
		FyY(B.numCellsX[0]),
		FyU((B.numCellsX[0]+1)*B.numCellsY[0]);
	
	real *FyX_r = thrust::raw_pointer_cast(&FyX[0]),
	     *FyY_r = thrust::raw_pointer_cast(&FyY[0]),
	     *FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, q_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, q_r, lambda_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );
	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, q_r, qOld_r, dxD, dyD, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	forceY = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());
}

template <typename memoryType>
void NSWithBody<memoryType>::writeCommon()
{	
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;
	std::string  caseFolder = db["inputs"]["caseFolder"].get<std::string>();
	int          nsave = db["simulation"]["nsave"].get<int>();
	int          timeStep = NavierStokesSolver<memoryType>::timeStep;

	NavierStokesSolver<memoryType>::writeCommon();
	
	// write the coordinates of the body points
	if (timeStep % nsave == 0)
		B.writeToFile(caseFolder, timeStep);

	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}

template <typename memoryType>
void NSWithBody<memoryType>::shutDown()
{
	io::printTimingInfo(NavierStokesSolver<memoryType>::logger);
	forceFile.close();
	NavierStokesSolver<memoryType>::iterationsFile.close();
}

template class NSWithBody<host_memory>;
template class NSWithBody<device_memory>;
