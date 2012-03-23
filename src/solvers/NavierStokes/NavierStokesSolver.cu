#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::initialise()
{
	printf("NS initalising\n");
	timeStep = simPar->startStep;
	
	//io::createDirectory
	mkdir(opts->folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	io::writeGrid(opts->folderName, *domInfo);
	
	initialiseArrays();
	assembleMatrices();
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::assembleMatrices()
{
	std::cout << "Entered assembleMatrices" << std::endl;
	generateA();
	std::cout << "Assembled A!" << std::endl;
	generateBN();
	std::cout << "Assembled BN!" << std::endl;
	generateQT();
	std::cout << "Assembled QT!" << std::endl;
	generateC(); // QT*BN*Q
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::initialiseArrays()
{
	int nx = domInfo->nx,
	    ny = domInfo->ny;
		
	int numU  = (nx-1)*ny;
	int numUV = numU + nx*(ny-1);
	int numP  = numU + ny;
	
	q.resize(numUV);
	qStar.resize(numUV);
	rn.resize(numUV);
	H.resize(numUV);
	bc1.resize(numUV);
	rhs1.resize(numUV);
	temp1.resize(numUV);
	
	cusp::blas::fill(rn, 0.0);
	cusp::blas::fill(H, 0.0);
	cusp::blas::fill(bc1, 0.0);
	cusp::blas::fill(rhs1, 0.0);
	cusp::blas::fill(temp1, 0.0);
	
	//lambda.resize(numP+2*numB);
	//rhs2.resize(numP+2*numB);
	lambda.resize(numP);
	bc2.resize(numP);
	rhs2.resize(numP);
	temp2.resize(numP);
	
	cusp::blas::fill(lambda, 0.0);
	cusp::blas::fill(bc2, 0.0);
	cusp::blas::fill(rhs2, 0.0);
	cusp::blas::fill(temp2, 0.0);
	
	/// Initialise velocity fluxes
	int i;
	for(i=0; i < numU; i++)
	{
		q[i] = flowDesc->initialU * domInfo->dy[i/(nx-1)];
	}
	for(; i < numUV; i++)
	{
		q[i] = flowDesc->initialV * domInfo->dx[(i-numU)%nx];
	}
	qStar = q;
	
	initialiseBoundaryArrays();
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::initialiseBoundaryArrays()
{
	int nx = domInfo->nx,
		ny = domInfo->ny;
		
	bcN[XMINUS].resize(2*ny-1);   
	bcN[XPLUS].resize(2*ny-1);
	bcN[YMINUS].resize(2*nx-1);
	bcN[YPLUS].resize(2*nx-1);

	/// Top and Bottom
	for(int i=0; i<nx-1; i++)
	{
		bcN[YMINUS][i] = flowDesc->bcInfo[YMINUS][0].second;
		bcN[YPLUS][i]  = flowDesc->bcInfo[YPLUS][0].second;
		bcN[YMINUS][i+nx-1]	= flowDesc->bcInfo[YMINUS][1].second;
		bcN[YPLUS][i+nx-1]	= flowDesc->bcInfo[YPLUS][1].second;
	}
	bcN[YMINUS][2*nx-2]	= flowDesc->bcInfo[YMINUS][1].second;
	bcN[YPLUS][2*nx-2]	= flowDesc->bcInfo[YPLUS][1].second;
	
	/// Left and Right
	for(int i=0; i<ny-1; i++)
	{
		bcN[XMINUS][i] = flowDesc->bcInfo[XMINUS][0].second;
		bcN[XPLUS][i]  = flowDesc->bcInfo[XPLUS][0].second;
		bcN[XMINUS][i+ny] = flowDesc->bcInfo[XMINUS][1].second;
		bcN[XPLUS][i+ny]  = flowDesc->bcInfo[XPLUS][1].second;
	}
	bcN[XMINUS][ny-1] = flowDesc->bcInfo[XMINUS][0].second;
	bcN[XPLUS][ny-1]  = flowDesc->bcInfo[XPLUS][0].second;
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateA()
{
	generateM();
	generateL();
	cusp::wrapped::subtract(M, L, A);
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateBN()
{
	BN = Minv;
}

/*
template <typename Matrix, typename Vector>
template <>
void NavierStokesSolver<Matrix, Vector>::generateBN<3>()
{
	Matrix	temp1, temp2;
	cusp::multiply(Minv, L, temp1);
	cusp::multiply(temp1, Minv, BN);
	cusp::add(Minv, BN, BN);
	cusp::multiply(temp1, BN, temp2);
	cusp::add(Minv, temp2, BN);
}*/

template <>
void NavierStokesSolver<cooD, vecD>::generateC()
{
	// Should this temp matrix be created each time step?
	cooD temp;
	cusp::wrapped::multiply(QT, BN, temp);
	cusp::wrapped::multiply(temp, Q, C);
	C.values[0] += C.values[0];
}

template <>
void NavierStokesSolver<cooH, vecH>::generateC()
{
	cooH temp;
	cusp::wrapped::multiply(QT, BN, temp);
	cusp::wrapped::multiply(temp, Q, C);
	C.sort_by_row_and_column();
	C.values[0] += C.values[0];
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::stepTime()
{
	//for(int i=0; i<simPar->intSch.substeps; i++)
	//std::cout << timeStep << ", ";
	for(int i=0; i<1; i++)
	{
		updateSolverState();

		generateRN();
		generateBC1();
		assembleRHS1();

		solveIntermediateVelocity();

		generateBC2();
		assembleRHS2();

		solvePoisson();

		projectionStep();
	}
	timeStep++;
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::assembleRHS1()
{
	cusp::blas::axpby(rn, bc1, rhs1, 1.0, 1.0);
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::solveIntermediateVelocity()
{
	cusp::default_monitor<real> sys1Mon(rhs1, 10000);
	cusp::krylov::cg(A, qStar, rhs1, sys1Mon);//, PC1);
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::assembleRHS2()
{
	cusp::wrapped::multiply(QT, qStar, temp2);
	cusp::blas::axpby(temp2, bc2, rhs2, 1.0, -1.0 );
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::solvePoisson()
{
	cusp::default_monitor<real> sys2Mon(rhs2, 20000);
	cusp::krylov::cg(C, lambda, rhs2, sys2Mon);//, PC2);
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::projectionStep()
{
	cusp::wrapped::multiply(Q, lambda, temp1);
	cusp::wrapped::multiply(BN, temp1, q);
	cusp::blas::axpby(qStar, q, q, 1.0, -1.0 );
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::writeData()
{
	if (timeStep % simPar->nsave == 0)
	{
		io::writeData(opts->folderName, timeStep, q, lambda, *domInfo);
	}
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::updateBoundaryConditions()
{
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::updateSolverState()
{
	updateBoundaryConditions();
}

template <typename Matrix, typename Vector>
bool NavierStokesSolver<Matrix, Vector>::finished()
{
	return (timeStep < simPar->nt) ? false : true;
}

/**
* \brief Factory method to select the required IBM solver
* \param a Description
* \return Pointer to an instance of the required dervied class.
*/
template <typename Matrix, typename Vector>
NavierStokesSolver<Matrix, Vector>* NavierStokesSolver<Matrix, Vector>::createSolver(options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
{
	NavierStokesSolver<Matrix, Vector> *solver = 0;
	switch(sim_par.ibmSch)
	{
		case SAIKI_BIRINGEN:
			break;
		case TAIRA_COLONIUS:
			solver = new TairaColoniusSolver<Matrix, Vector>;
			break;
		case NAVIER_STOKES:
			solver = new NavierStokesSolver<Matrix, Vector>;
			break;
		case FADLUN_ET_AL:
			solver = new FadlunEtAlSolver<Matrix, Vector>;
			break;
	}
	solver->opts = &opts;
	solver->flowDesc = &flow_desc;
	solver->simPar = &sim_par;
	solver->domInfo = &dom_info;
	std::cout << "Selected solver: " << solver->name() << std::endl;
	return solver;
}

#include "NavierStokes/generateM.inl"
#include "NavierStokes/generateL.inl"
#include "NavierStokes/generateQT.inl"
#include "NavierStokes/generateRN.inl"
#include "NavierStokes/generateBC1.inl"
#include "NavierStokes/generateBC2.inl"

template class NavierStokesSolver<cooH, vecH>;
template class NavierStokesSolver<cooD, vecD>;