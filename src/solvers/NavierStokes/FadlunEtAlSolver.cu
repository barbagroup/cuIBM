#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <sys/stat.h>

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	initialiseBodies();
	
	NavierStokesSolver<memoryType>::initialiseCommon();
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
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	
	tagPoints();
	std::cout << "Done tagging points!" << std::endl;
	
	NavierStokesSolver<memoryType>::assembleMatrices();
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateBodies()
{
}

/*template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateSolverState()
{	
	NavierStokesSolver<memoryType>::updateBoundaryConditions();
	if (B.bodiesMove)
	{
		updateBodies();
		//updateA();
	}
}*/

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateRN()
{
	NavierStokesSolver<memoryType>::generateRNFull();
	updateRN();
	
	/**
	* Added for testing. REMOVE SOON!!!
	*/
	//cusp::blas::fill(NavierStokesSolver<memoryType>::rn, 0.0);
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateBC1()
{
	NavierStokesSolver<memoryType>::generateBC1Full(NavierStokesSolver<memoryType>::intgSchm.alphaImplicit[NavierStokesSolver<memoryType>::subStep]);
	updateBC1();
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateBC1()
{
}

#include "FadlunEtAl/initialiseBodies.inl"
#include "FadlunEtAl/tagPoints.inl"
#include "FadlunEtAl/generateL.inl"
#include "FadlunEtAl/updateRN.inl"
#include "FadlunEtAl/generateQT.inl"

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
