#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <sys/stat.h>

template <typename memoryType>
void NavierStokesSolver<memoryType>::initialise()
{
	printf("NS initalising\n");
	
	int nx = domInfo->nx,
	    ny = domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	initialiseCommon();
	initialiseArrays(numUV, numP);
	assembleMatrices();
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseCommon()
{
	QCoeff = 1.0;
	
	timeScheme convScheme = (*paramDB)["simulation"]["convTimeScheme"].get<timeScheme>(),
	           diffScheme = (*paramDB)["simulation"]["diffTimeScheme"].get<timeScheme>();
	intgSchm.initialise(convScheme, diffScheme);
	
	// initial values of timeStep
	timeStep = (*paramDB)["simulation"]["startStep"].get<int>();
	
	// create directory 
	std::string folderName = (*paramDB)["inputs"]["folderName"].get<std::string>();
	mkdir(folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	// write the grids information to a file
	io::writeGrid(folderName, *domInfo);
	
	std::cout << "Initialised common stuff!" << std::endl;
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseArrays(int numQ, int numLambda)
{	
	q.resize(numQ);
	qStar.resize(numQ);
	rn.resize(numQ);
	H.resize(numQ);
	bc1.resize(numQ);
	rhs1.resize(numQ);
	temp1.resize(numQ);
	
	cusp::blas::fill(rn, 0.0);
	cusp::blas::fill(H, 0.0);
	cusp::blas::fill(bc1, 0.0);
	cusp::blas::fill(rhs1, 0.0);
	cusp::blas::fill(temp1, 0.0);
	
	lambda.resize(numLambda);
	bc2.resize(numLambda);
	rhs2.resize(numLambda);
	temp2.resize(numLambda);
	
	cusp::blas::fill(lambda, 0.0);
	cusp::blas::fill(bc2, 0.0);
	cusp::blas::fill(rhs2, 0.0);
	cusp::blas::fill(temp2, 0.0);
	
	initialiseFluxes();
	initialiseBoundaryArrays();
	
	generateRNFull(0);
	cusp::blas::scal(H, 1.0/intgSchm.gamma[0]);
	std::cout << "Initialised arrays!" << std::endl;
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

	/// Top and Bottom
	for(int i=0; i<nx-1; i++)
	{
		bc[YMINUS][i] = bcInfo[YMINUS][0].value;
		bc[YPLUS][i]  = bcInfo[YPLUS][0].value;
		bc[YMINUS][i+nx-1]	= bcInfo[YMINUS][1].value;
		bc[YPLUS][i+nx-1]	= bcInfo[YPLUS][1].value;
	}
	bc[YMINUS][2*nx-2]	= bcInfo[YMINUS][1].value;
	bc[YPLUS][2*nx-2]	= bcInfo[YPLUS][1].value;
	
	/// Left and Right
	for(int i=0; i<ny-1; i++)
	{
		bc[XMINUS][i] = bcInfo[XMINUS][0].value;
		bc[XPLUS][i]  = bcInfo[XPLUS][0].value;
		bc[XMINUS][i+ny] = bcInfo[XMINUS][1].value;
		bc[XPLUS][i+ny]  = bcInfo[XPLUS][1].value;
	}
	bc[XMINUS][ny-1] = bcInfo[XMINUS][0].value;
	bc[XPLUS][ny-1]  = bcInfo[XPLUS][0].value;
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleMatrices()
{
	std::cout << "Entered assembleMatrices" << std::endl;
	generateM();
	generateL();
	generateA(intgSchm.alphaImplicit[0]);
	std::cout << "Assembled A!" << std::endl;
	////
	//cusp::print(L);
	//cusp::print(A);

//	PC1 = cusp::precond::diagonal<real, memoryType>(A);
	
	generateBN();
	std::cout << "Assembled BN!" << std::endl;
	
	generateQT();
	std::cout << "Assembled QT!" << std::endl;
	//cusp::print(QT);
	
	generateC(); // QT*BN*Q
	std::cout << "Generated C!" << std::endl;
	//cusp::print(C);

//	PC2 = cusp::precond::smoothed_aggregation<int, real, memoryType>(C);

	std::cout << "Assembled matrices!" << std::endl;
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
	cooD temp; // Should this temp matrix be created each time step?
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
	//std::cout << timeStep << std::endl;
	for(int i=0; i<intgSchm.subSteps; i++)
	{
		updateSolverState(i);

		generateRN(i);
		generateBC1(i);
		assembleRHS1();
		
		//cusp::print(rn);
//		cusp::print(rhs1);

		solveIntermediateVelocity();

		generateBC2();
		assembleRHS2();

		solvePoisson();

		projectionStep();
//		cusp::print(q);
	}
	
	timeStep++;
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::updateSolverState(int i)
{
//	generateA(intgSchm.alphaImplicit[i]);
//	updateQ(intgSchm.gamma[i]);
//	updateBoundaryConditions();
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::generateRN(int i)
{
	generateRNFull(i);
	/**
	* Does this include the pressure term on the RHS?
	*/
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::generateBC1(int i)
{
	generateBC1Full(intgSchm.alphaImplicit[i]);
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
	cusp::krylov::bicgstab(A, qStar, rhs1, sys1Mon);//, PC1);
	//cusp::krylov::cg(A, qStar, rhs1, sys1Mon);//, PC1);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleRHS2()
{
	cusp::wrapped::multiply(QT, qStar, temp2);
	cusp::blas::axpby(temp2, bc2, rhs2, 1.0, -1.0);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::solvePoisson()
{
	cusp::default_monitor<real> sys2Mon(rhs2, 30000);
	//cusp::krylov::gmres(C, lambda, rhs2, 50, sys2Mon);//, PC2);
	//cusp::krylov::bicgstab(C, lambda, rhs2, sys2Mon);//, PC2);
	cusp::krylov::cg(C, lambda, rhs2, sys2Mon);//, PC2);
	if (!sys2Mon.converged())
	{
		std::cout << "ERROR: Solve for Lambda failed at time step " << timeStep << std::endl;
		std::exit(-1);
	}
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::projectionStep()
{
	cusp::wrapped::multiply(Q, lambda, temp1);
	cusp::wrapped::multiply(BN, temp1, q);
	cusp::blas::axpby(qStar, q, q, 1.0, -1.0);
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
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
//	calculateForce();
//	io::writeForce(folderName, timeStep*dt, forceX, forceY);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::updateQ(real gamma)
{
	cusp::blas::scal(Q.values, gamma/QCoeff);
	QCoeff = gamma;
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::updateBoundaryConditions()
{
}

template <typename memoryType>
bool NavierStokesSolver<memoryType>::finished()
{
	int nt = (*paramDB)["simulation"]["nt"].get<int>();
	return (timeStep < nt) ? false : true;
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::calculateForce()
{
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
#include "NavierStokes/generateA.inl"
#include "NavierStokes/generateQT.inl"
#include "NavierStokes/generateRN.inl"
#include "NavierStokes/generateBC1.inl"
#include "NavierStokes/generateBC2.inl"

template class NavierStokesSolver<host_memory>;
template class NavierStokesSolver<device_memory>;
