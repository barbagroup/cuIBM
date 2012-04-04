#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <sys/stat.h>

template <typename memoryType>
void NavierStokesSolver<memoryType>::initialise()
{
	printf("NS initalising\n");
	timeStep = (*paramDB)["simulation"]["startStep"].get<int>();
	
	//io::createDirectory
  std::string folderName = (*paramDB)["inputs"]["folderName"].get<std::string>();
	mkdir(folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	io::writeGrid(folderName, *domInfo);
	
	initialiseArrays();
	assembleMatrices();
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleMatrices()
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

template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseArrays()
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
	//bc2Host.resize(numP);
	rhs2.resize(numP);
	temp2.resize(numP);
	
	cusp::blas::fill(lambda, 0.0);
	cusp::blas::fill(bc2, 0.0);
	//cusp::blas::fill(bc2Host, 0.0);
	cusp::blas::fill(rhs2, 0.0);
	cusp::blas::fill(temp2, 0.0);
	
	initialiseFluxes();
	initialiseBoundaryArrays();
}

template <>
void NavierStokesSolver <host_memory>::initialiseFluxes()
{
	int nx = domInfo->nx,
	    ny = domInfo->ny;
	int numU  = (nx-1)*ny;
	int numUV = numU + nx*(ny-1);
	int i;
  real uInitial, vInitial;
  uInitial = (*paramDB)["flow"]["uInitial"].get<real>();
  vInitial = (*paramDB)["flow"]["vInitial"].get<real>();
	for(i=0; i < numU; i++)
	{
		q[i] = uInitial * domInfo->dy[i/(nx-1)];
	}
	for(; i < numUV; i++)
	{
		q[i] = vInitial * domInfo->dx[(i-numU)%nx];
	}
	qStar = q;
}

template<>
void NavierStokesSolver <device_memory>::initialiseFluxes()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	vecH qHost(numUV);
	int  i;
  real uInitial, vInitial;
  uInitial = (*paramDB)["flow"]["uInitial"].get<real>();
  vInitial = (*paramDB)["flow"]["vInitial"].get<real>();
	for(i=0; i < numU; i++)
	{
		qHost[i] = uInitial * domInfo->dy[i/(nx-1)];
	}
	for(; i < numUV; i++)
	{
		qHost[i] = vInitial * domInfo->dx[(i-numU)%nx];
	}
	q = qHost;
	qStar = q;
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseBoundaryArrays()
{
	int nx = domInfo->nx,
		ny = domInfo->ny;

  boundaryCondition **bcInfo = (*paramDB)["flow"]["boundaryConditions"].get<boundaryCondition **>();
	
	bc[XMINUS].resize(2*ny-1);
	bc[XPLUS].resize(2*ny-1);
	bc[YMINUS].resize(2*nx-1);
	bc[YPLUS].resize(2*nx-1);
	bcHost[XMINUS].resize(2*ny-1);   
	bcHost[XPLUS].resize(2*ny-1);
	bcHost[YMINUS].resize(2*nx-1);
	bcHost[YPLUS].resize(2*nx-1);

	/// Top and Bottom
	for(int i=0; i<nx-1; i++)
	{
		bcHost[YMINUS][i] = bcInfo[YMINUS][0].value;
		bcHost[YPLUS][i]  = bcInfo[YPLUS][0].value;
		bcHost[YMINUS][i+nx-1]	= bcInfo[YMINUS][1].value;
		bcHost[YPLUS][i+nx-1]	= bcInfo[YPLUS][1].value;
	}
	bcHost[YMINUS][2*nx-2]	= bcInfo[YMINUS][1].value;
	bcHost[YPLUS][2*nx-2]	= bcInfo[YPLUS][1].value;
	
	/// Left and Right
	for(int i=0; i<ny-1; i++)
	{
		bcHost[XMINUS][i] = bcInfo[XMINUS][0].value;
		bcHost[XPLUS][i]  = bcInfo[XPLUS][0].value;
		bcHost[XMINUS][i+ny] = bcInfo[XMINUS][1].value;
		bcHost[XPLUS][i+ny]  = bcInfo[XPLUS][1].value;
	}
	bcHost[XMINUS][ny-1] = bcInfo[XMINUS][0].value;
	bcHost[XPLUS][ny-1]  = bcInfo[XPLUS][0].value;
	
	bc[XMINUS] = bcHost[XMINUS];
	bc[XPLUS]  = bcHost[XPLUS];
	bc[YMINUS] = bcHost[YMINUS];
	bc[YPLUS]  = bcHost[YPLUS];
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::generateA()
{
	generateM();
	generateL();
	cusp::wrapped::subtract(M, L, A);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::generateBN()
{
	BN = Minv;
}

/*
template <typename memoryType>
template <>
void NavierStokesSolver<memoryType>::generateBN<3>()
{
	Matrix	temp1, temp2;
	cusp::multiply(Minv, L, temp1);
	cusp::multiply(temp1, Minv, BN);
	cusp::add(Minv, BN, BN);
	cusp::multiply(temp1, BN, temp2);
	cusp::add(Minv, temp2, BN);
}*/

template <>
void NavierStokesSolver<device_memory>::generateC()
{
	// Should this temp matrix be created each time step?
	cooD temp;
	cusp::wrapped::multiply(QT, BN, temp);
	cusp::wrapped::multiply(temp, Q, C);
	C.values[0] += C.values[0];
}

template <>
void NavierStokesSolver<host_memory>::generateC()
{
	cooH temp;
	cusp::wrapped::multiply(QT, BN, temp);
	cusp::wrapped::multiply(temp, Q, C);
	C.sort_by_row_and_column();
	C.values[0] += C.values[0];
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::stepTime()
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

template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleRHS1()
{
	cusp::blas::axpby(rn, bc1, rhs1, 1.0, 1.0);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::solveIntermediateVelocity()
{
	cusp::default_monitor<real> sys1Mon(rhs1, 10000);
	cusp::krylov::cg(A, qStar, rhs1, sys1Mon);//, PC1);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleRHS2()
{
	cusp::wrapped::multiply(QT, qStar, temp2);
	cusp::blas::axpby(temp2, bc2, rhs2, 1.0, -1.0 );
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::solvePoisson()
{
	cusp::default_monitor<real> sys2Mon(rhs2, 20000);
	cusp::krylov::cg(C, lambda, rhs2, sys2Mon);//, PC2);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::projectionStep()
{
	cusp::wrapped::multiply(Q, lambda, temp1);
	cusp::wrapped::multiply(BN, temp1, q);
	cusp::blas::axpby(qStar, q, q, 1.0, -1.0 );
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::writeData()
{
  int nsave = (*paramDB)["simulation"]["nsave"].get<int>();
  std::string folderName = (*paramDB)["inputs"]["folderName"].get<std::string>();
	if (timeStep % nsave == 0)
	{
		io::writeData(folderName, timeStep, q, lambda, *domInfo);
	}
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::updateBoundaryConditions()
{
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::updateSolverState()
{
	updateBoundaryConditions();
}

template <typename memoryType>
bool NavierStokesSolver<memoryType>::finished()
{
  int nt = (*paramDB)["simulation"]["nt"].get<int>();
	return (timeStep < nt) ? false : true;
}

/**
* \brief Factory method to select the required IBM solver
* \param a Description
* \return Pointer to an instance of the required dervied class.
*/
template <typename memoryType>
NavierStokesSolver<memoryType>* NavierStokesSolver<memoryType>::createSolver(parameterDB &paramDB, domain &dom_info)
{
	ibmScheme ibm = paramDB["simulation"]["ibmScheme"].get<ibmScheme>();
	NavierStokesSolver<memoryType> *solver = 0;
	switch(ibm)
	{
		case SAIKI_BIRINGEN:
			break;
		case TAIRA_COLONIUS:
			solver = new TairaColoniusSolver<memoryType>;
			break;
		case NAVIER_STOKES:
			solver = new NavierStokesSolver<memoryType>;
			break;
		case FADLUN_ET_AL:
			solver = new FadlunEtAlSolver<memoryType>;
			break;
	}
  solver->paramDB = &paramDB;
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

template class NavierStokesSolver<host_memory>;
template class NavierStokesSolver<device_memory>;
