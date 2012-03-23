#include <solvers/NavierStokes/TairaColoniusSolver.h>

template <typename Matrix, typename Vector>
void TairaColoniusSolver<Matrix, Vector>::updateBodies()
{
}
template <typename Matrix, typename Vector>
void TairaColoniusSolver<Matrix, Vector>::initialise()
{
	initialiseBodies();
	printf("NS initalising\n");
	NavierStokesSolver<Matrix, Vector>::timeStep = NavierStokesSolver<Matrix, Vector>::simPar->startStep;
		
	//io::createDirectory
	mkdir(NavierStokesSolver<Matrix, Vector>::opts->folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	io::writeGrid(NavierStokesSolver<Matrix, Vector>::opts->folderName, *NavierStokesSolver<Matrix, Vector>::domInfo);
	
	printf("TC: initalising\n");
	initialiseArrays();
	
	std::cout << "Initialised bodies!" << std::endl;
	NavierStokesSolver<Matrix, Vector>::assembleMatrices();
	std::cout << "Assembled matrices!" << std::endl;
}

template <typename Matrix, typename Vector>
void TairaColoniusSolver<Matrix, Vector>::initialiseArrays()
{
	int nx = NavierStokesSolver<Matrix, Vector>::domInfo->nx,
	    ny = NavierStokesSolver<Matrix, Vector>::domInfo->ny;
		
	int numU  = (nx-1)*ny;
	int numUV = numU + nx*(ny-1);
	int numP  = numU + ny;
//	int numB  = 0;
//	for(int k=0; k<flowDesc->numBodies; k++)
//		numB += flowDesc->B[k].numPoints;
	
	NavierStokesSolver<Matrix, Vector>::q.resize(numUV);
	NavierStokesSolver<Matrix, Vector>::qStar.resize(numUV);
	NavierStokesSolver<Matrix, Vector>::rn.resize(numUV);
	NavierStokesSolver<Matrix, Vector>::H.resize(numUV);
	NavierStokesSolver<Matrix, Vector>::bc1.resize(numUV);
	NavierStokesSolver<Matrix, Vector>::rhs1.resize(numUV);
	NavierStokesSolver<Matrix, Vector>::temp1.resize(numUV);
	
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::rn, 0.0);
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::H, 0.0);
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::bc1, 0.0);
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::rhs1, 0.0);
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::temp1, 0.0);
	
	//lambda.resize(numP+2*numB);
	//rhs2.resize(numP+2*numB);
	NavierStokesSolver<Matrix, Vector>::lambda.resize(numP+2*B.totalPoints);
	NavierStokesSolver<Matrix, Vector>::bc2.resize(numP+2*B.totalPoints);
	NavierStokesSolver<Matrix, Vector>::rhs2.resize(numP+2*B.totalPoints);
	NavierStokesSolver<Matrix, Vector>::temp2.resize(numP+2*B.totalPoints);
	
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::lambda, 0.0);
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::bc2, 0.0);
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::rhs2, 0.0);
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::temp2, 0.0);
	
	/// Initialise velocity fluxes
	int i;
	for(i=0; i < numU; i++)
	{
		NavierStokesSolver<Matrix, Vector>::q[i] = NavierStokesSolver<Matrix, Vector>::flowDesc->initialU * NavierStokesSolver<Matrix, Vector>::domInfo->dy[i/(nx-1)];
	}
	for(; i < numUV; i++)
	{
		NavierStokesSolver<Matrix, Vector>::q[i] = NavierStokesSolver<Matrix, Vector>::flowDesc->initialV * NavierStokesSolver<Matrix, Vector>::domInfo->dx[(i-numU)%nx];
	}
	NavierStokesSolver<Matrix, Vector>::qStar = NavierStokesSolver<Matrix, Vector>::q;
	
	NavierStokesSolver<Matrix, Vector>::initialiseBoundaryArrays();
}

template <typename Matrix, typename Vector>
void TairaColoniusSolver<Matrix, Vector>::updateSolverState()
{
	NavierStokesSolver<Matrix, Vector>::updateBoundaryConditions();
	if (B.bodiesMove) {
		updateBodies();
		updateQT();
		NavierStokesSolver<Matrix, Vector>::generateC();
	}
}

#include "TairaColonius/generateQT.inl"
#include "TairaColonius/generateBC2.inl"
#include "TairaColonius/initialiseBodies.inl"

template class TairaColoniusSolver<cooH, vecH>;
template class TairaColoniusSolver<cooD, vecD>;