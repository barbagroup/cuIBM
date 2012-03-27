#include <solvers/NavierStokes/TairaColoniusSolver.h>

template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateBodies()
{
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialise()
{
	initialiseBodies();
	printf("NS initalising\n");
	NavierStokesSolver<memoryType>::timeStep = NavierStokesSolver<memoryType>::simPar->startStep;
		
	//io::createDirectory
	mkdir(NavierStokesSolver<memoryType>::opts->folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	io::writeGrid(NavierStokesSolver<memoryType>::opts->folderName, *NavierStokesSolver<memoryType>::domInfo);
	
	printf("TC: initalising\n");
	initialiseArrays();
	
	std::cout << "Initialised bodies!" << std::endl;
	NavierStokesSolver<memoryType>::assembleMatrices();
	std::cout << "Assembled matrices!" << std::endl;
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialiseArrays()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
	    ny = NavierStokesSolver<memoryType>::domInfo->ny;
		
	int numU  = (nx-1)*ny;
	int numUV = numU + nx*(ny-1);
	int numP  = numU + ny;
//	int numB  = 0;
//	for(int k=0; k<flowDesc->numBodies; k++)
//		numB += flowDesc->B[k].numPoints;
	
	NavierStokesSolver<memoryType>::q.resize(numUV);
	NavierStokesSolver<memoryType>::qStar.resize(numUV);
	NavierStokesSolver<memoryType>::rn.resize(numUV);
	NavierStokesSolver<memoryType>::H.resize(numUV);
	NavierStokesSolver<memoryType>::bc1.resize(numUV);
	NavierStokesSolver<memoryType>::rhs1.resize(numUV);
	NavierStokesSolver<memoryType>::temp1.resize(numUV);
	
	cusp::blas::fill(NavierStokesSolver<memoryType>::rn, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::H, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::bc1, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::rhs1, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::temp1, 0.0);
	
	//lambda.resize(numP+2*numB);
	//rhs2.resize(numP+2*numB);
	NavierStokesSolver<memoryType>::lambda.resize(numP+2*B.totalPoints);
	NavierStokesSolver<memoryType>::bc2.resize(numP+2*B.totalPoints);
	NavierStokesSolver<memoryType>::bc2Host.resize(numP+2*B.totalPoints);
	NavierStokesSolver<memoryType>::rhs2.resize(numP+2*B.totalPoints);
	NavierStokesSolver<memoryType>::temp2.resize(numP+2*B.totalPoints);
	
	cusp::blas::fill(NavierStokesSolver<memoryType>::lambda, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::bc2, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::bc2Host, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::rhs2, 0.0);
	cusp::blas::fill(NavierStokesSolver<memoryType>::temp2, 0.0);
	
	NavierStokesSolver<memoryType>::initialiseFluxes();
	NavierStokesSolver<memoryType>::initialiseBoundaryArrays();
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateSolverState()
{
	NavierStokesSolver<memoryType>::updateBoundaryConditions();
	if (B.bodiesMove) {
		updateBodies();
		updateQT();
		NavierStokesSolver<memoryType>::generateC();
	}
}

#include "TairaColonius/generateQT.inl"
#include "TairaColonius/generateBC2.inl"
#include "TairaColonius/initialiseBodies.inl"

template class TairaColoniusSolver<host_memory>;
template class TairaColoniusSolver<device_memory>;