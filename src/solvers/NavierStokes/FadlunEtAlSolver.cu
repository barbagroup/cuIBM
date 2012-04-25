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
	
	tags.resize(numUV);
	tagsD.resize(numUV);
	coeffs.resize(numUV);
	coeffsD.resize(numUV);
	tagPoints();
	std::cout << "Tagged points!" << std::endl;
	
	NavierStokesSolver<memoryType>::assembleMatrices();
	std::cout << "Assembled matrices!" << std::endl;
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateBodies()
{
}

/*template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateSolverState(int i)
{	
	NavierStokesSolver<memoryType>::updateBoundaryConditions();
	if (B.bodiesMove)
	{
		updateBodies();
		//updateA();
	}
}*/

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateRN(int i)
{
	NavierStokesSolver<memoryType>::generateRNFull(i);
	updateRN();
	
	/**
	* Added for testing. REMOVE SOON!!!
	*/
	//cusp::blas::fill(NavierStokesSolver<memoryType>::rn, 0.0);
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateBC1(int i)
{
	NavierStokesSolver<memoryType>::generateBC1Full(NavierStokesSolver<memoryType>::intgSchm.alphaImplicit[i]);
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

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
