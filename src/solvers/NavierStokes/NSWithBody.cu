/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

#include <solvers/NavierStokes/NSWithBody.h>
#include <solvers/NavierStokes/kernels/calculateForce.h>

/// Initialise the bodies
template <typename memoryType>
void NSWithBody<memoryType>::initialiseBodies()
{
	NavierStokesSolver<memoryType>::logger.startTimer("initialiseBodies");
	
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	B.initialise(db, *NavierStokesSolver<memoryType>::domInfo);

	std::string folderName = db["inputs"]["folderName"].get<std::string>();
	std::stringstream out;
	out << folderName << "/forces";
	forceFile.open(out.str().c_str());
	
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
void NSWithBody<memoryType>::writeData()
{
	NavierStokesSolver<memoryType>::logger.startTimer("output");
	
	// write the velocity and pressure
	NavierStokesSolver<memoryType>::writeCommon();
	
	// calculate and write the force on the body
	calculateForce();
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	real dt = db["simulation"]["dt"].get<real>();
	forceFile << NavierStokesSolver<memoryType>::timeStep*dt << '\t' << forceX << '\t' << forceY << '\n';
	
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
