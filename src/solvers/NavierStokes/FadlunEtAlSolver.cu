#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <sys/stat.h>

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialise()
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

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateSolverState()
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

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::assembleRHS1()
{
	NavierStokesSolver<memoryType>::assembleRHS1();
	
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
	updateRHS1();
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
}

#include "FadlunEtAl/tagPoints.inl"
#include "FadlunEtAl/generateL.inl"
//
#include "FadlunEtAl/generateA.inl"
//
#include "FadlunEtAl/updateRHS1.inl"
#include "FadlunEtAl/generateQT.inl"
#include "FadlunEtAl/calculateForce.inl"

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
