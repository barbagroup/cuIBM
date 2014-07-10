/***************************************************************************//**
* \file DirectForcingSolver.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the methods of the class \c DirectForcingSolver
*/

#include <solvers/NavierStokes/DirectForcingSolver.h>
#include <sys/stat.h>

/**
* \brief To be documented
*/
template <typename memoryType>
void DirectForcingSolver<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	
	NSWithBody<memoryType>::initialiseBodies();
	
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
/*	
	tags.resize(numUV);
	tagsD.resize(numUV);
	coeffs.resize(numUV);
	coeffsD.resize(numUV);
*/
	tagsX.resize(numUV);
	tagsXD.resize(numUV);
	tagsY.resize(numUV);
	tagsYD.resize(numUV);
	coeffsX.resize(numUV);
	coeffsXD.resize(numUV);
	coeffsY.resize(numUV);
	coeffsYD.resize(numUV);
	uvX.resize(numUV);
	uvY.resize(numUV);
	uvXD.resize(numUV);
	uvYD.resize(numUV);
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	
	tagPoints();
	std::cout << "Done tagging points!" << std::endl;
	
	NavierStokesSolver<memoryType>::assembleMatrices();
}

/**
* \brief To be documented
*/
template <typename memoryType>
void DirectForcingSolver<memoryType>::updateSolverState()
{
	if (NSWithBody<memoryType>::B.bodiesMove)
	{
		// update the locations of the body points
		NSWithBody<memoryType>::updateBodies();
		
		// retag points
		tagPoints();
		
		// assemble the matrices generated using new tags
		NavierStokesSolver<memoryType>::assembleMatrices();
		
/*		logger.startTimer("assembleMatrices");
		
		generateL();
		generateA(intgSchm.alphaImplicit[subStep]);
		PC1 = new preconditioner< cusp::coo_matrix<int, real, memoryType> >(A, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>());
		generateBN();	
	
		logger.stopTimer("assembleMatrices");

		generateQT();
		generateC(); // QT*BN*Q
	
		logger.startTimer("preconditioner2");
		PC2 = new preconditioner< cusp::coo_matrix<int, real, memoryType> >(C, (*paramDB)["PoissonSolve"]["preconditioner"].get<preconditionerType>());
		logger.stopTimer("preconditioner2");*/
	}
}

/**
* \brief To be documented
*/
template <typename memoryType>
void DirectForcingSolver<memoryType>::assembleRHS1()
{
	NavierStokesSolver<memoryType>::assembleRHS1();
	
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
	updateRHS1();
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
}

/**
* \brief To be documented
*/
template <typename memoryType>
void DirectForcingSolver<memoryType>::writeData()
{	
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;
	real         dt  = db["simulation"]["dt"].get<real>();
	int          timeStep = NavierStokesSolver<memoryType>::timeStep;

	NSWithBody<memoryType>::writeCommon();
	
	// Print forces calculated using the CV approach
	NSWithBody<memoryType>::calculateForce();
	NSWithBody<memoryType>::forceFile << timeStep*dt << '\t' << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << std::endl;
	
	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}

/**
* \brief To be documented
*/
template <typename memoryType>
DirectForcingSolver<memoryType>::DirectForcingSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

#include "DirectForcing/tagPoints.inl"
#include "DirectForcing/generateL.inl"
#include "DirectForcing/generateA.inl"
#include "DirectForcing/updateRHS1.inl"
#include "DirectForcing/generateQT.inl"

// specialization of the class
template class DirectForcingSolver<host_memory>;
template class DirectForcingSolver<device_memory>;
