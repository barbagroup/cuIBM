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

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <sys/stat.h>

//##############################################################################
//                              INITIALISE
//##############################################################################

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
	logger.startTimer("initialiseCommon");
	
	QCoeff = 1.0;
	subStep = 0;
	
	timeScheme convScheme = (*paramDB)["simulation"]["convTimeScheme"].get<timeScheme>(),
	           diffScheme = (*paramDB)["simulation"]["diffTimeScheme"].get<timeScheme>();
	intgSchm.initialise(convScheme, diffScheme);
	
	/// initial values of timeStep
	timeStep = (*paramDB)["simulation"]["startStep"].get<int>();
	
	/// create directory 
	std::string folderName = (*paramDB)["inputs"]["folderName"].get<std::string>();
	mkdir(folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	/// write the grids information to a file
	io::writeGrid(folderName, *domInfo);
	
	/// open the required files
	std::stringstream out;
	out << folderName << "/forces";
	forceFile.open(out.str().c_str());
	out.str("");
	out << folderName << "/iterations";
	iterationsFile.open(out.str().c_str());
	
	/// write the plot information to a file
	
	std::cout << "Initialised common stuff!" << std::endl;
	
	logger.stopTimer("initialiseCommon");
}

/**
* \brief Initialises all required arrays
* \param numQ Total number velocity variables (u and v)
* \param numLambda Number of pressure variables + twice the number of body force variables
*/
template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseArrays(int numQ, int numLambda)
{	
	logger.startTimer("initialiseArrays");
	
	q.resize(numQ);
	qStar.resize(numQ);
		qOld.resize(numQ);
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
	
	generateRNFull();
	cusp::blas::scal(H, 1.0/intgSchm.gamma[subStep]);
	std::cout << "Initialised arrays!" << std::endl;
	
	logger.stopTimer("initialiseArrays");
}

/**
* \brief Sets the initial value of all the fluxes in the flow field
*/
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
	real uInitial, vInitial;
	int  i;
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

/**
* \brief Sets the initial values of the the boundary velocities
*/
template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseBoundaryArrays()
{
	int nx = domInfo->nx,
		ny = domInfo->ny;

	boundaryCondition 
		**bcInfo
	     = (*paramDB)["flow"]["boundaryConditions"].get<boundaryCondition **>();
	
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

//##############################################################################
//                          ASSEMBLE MATRICES
//##############################################################################

/**
* \brief Assembles all the required matrices
*/
template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleMatrices()
{
	logger.startTimer("assembleMatrices");
	
	generateM();
	generateL();
	generateA(intgSchm.alphaImplicit[subStep]);
	PC1 = new preconditioner< cusp::coo_matrix<int, real, memoryType>, cusp::array1d<real, memoryType> >(A, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>());
	generateBN();
	generateQT();
	generateC(); // QT*BN*Q
	PC2 = new preconditioner< cusp::coo_matrix<int, real, memoryType>, cusp::array1d<real, memoryType> >(C, (*paramDB)["PoissonSolve"]["preconditioner"].get<preconditionerType>());

	std::cout << "Assembled matrices!" << std::endl;	
	logger.stopTimer("assembleMatrices");
}

/**
* \brief Generates the approximate inverse of the matrix A
*/
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

//##############################################################################
//                            TIME STEPPING
//##############################################################################

template <typename memoryType>
void NavierStokesSolver<memoryType>::stepTime()
{
	qOld = q;
	for(subStep=0; subStep < intgSchm.subSteps; subStep++)
	{
		updateSolverState();

		// Set up and solve the first system for the intermediate velocity
		generateRN();
		generateBC1();
		assembleRHS1();
		solveIntermediateVelocity();

		// Set up and solve the Poisson system
		generateBC2();
		assembleRHS2();
		solvePoisson();

		// Projection step
		projectionStep();
	}
	
	timeStep++;
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::updateSolverState()
{
//	generateA(intgSchm.alphaImplicit[subStep]);
//	updateQ(intgSchm.gamma[i]);
//	updateBoundaryConditions();
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::updateQ(real gamma)
{
//	cusp::blas::scal(Q.values, gamma/QCoeff);
//	QCoeff = gamma;
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

//##############################################################################
//                          GENERATE VECTORS
//##############################################################################

template <typename memoryType>
void NavierStokesSolver<memoryType>::generateRN()
{
	generateRNFull();
	/**
	* Does this include the pressure term on the RHS?
	*/
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::generateBC1()
{
	generateBC1Full(intgSchm.alphaImplicit[subStep]);
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleRHS1()
{
	logger.startTimer("assembleRHS1");
	
	cusp::blas::axpby(rn, bc1, rhs1, 1.0, 1.0);
	
	logger.stopTimer("assembleRHS1");
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleRHS2()
{
	cusp::wrapped::multiply(QT, qStar, temp2);
	cusp::blas::axpby(temp2, bc2, rhs2, 1.0, -1.0);
}

//##############################################################################
//                           LINEAR SOLVES
//##############################################################################

template <typename memoryType>
void NavierStokesSolver<memoryType>::solveIntermediateVelocity()
{
	logger.startTimer("solveIntermediateVel");
	
	int  maxIters = (*paramDB)["velocitySolve"]["maxIterations"].get<int>();
	real relTol = (*paramDB)["velocitySolve"]["tolerance"].get<real>();

	cusp::default_monitor<real> sys1Mon(rhs1, maxIters, relTol);
	cusp::krylov::bicgstab(A, qStar, rhs1, sys1Mon, *PC1);
	//cusp::krylov::cg(A, qStar, rhs1, sys1Mon);//, PC1);
	iterationCount1 = sys1Mon.iteration_count();
	if (!sys1Mon.converged())
	{
		std::cout << "ERROR: Solve for q* failed at time step " << timeStep << std::endl;
		std::cout << "Iterations   : " << iterationCount1 << std::endl;          
		std::cout << "Residual norm: " << sys1Mon.residual_norm() << std::endl;
		std::cout << "Tolerance    : " << sys1Mon.tolerance() << std::endl;
		std::exit(-1);
	}

	logger.stopTimer("solveIntermediateVel");
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::solvePoisson()
{
	logger.startTimer("solvePoisson");
	
	int  maxIters = (*paramDB)["PoissonSolve"]["maxIterations"].get<int>();
	real relTol = (*paramDB)["PoissonSolve"]["tolerance"].get<real>();
	
	cusp::default_monitor<real> sys2Mon(rhs2, maxIters, relTol);
	//cusp::krylov::gmres(C, lambda, rhs2, 50, sys2Mon);//, PC2);
	//cusp::krylov::bicgstab(C, lambda, rhs2, sys2Mon);//, PC2);
	cusp::krylov::cg(C, lambda, rhs2, sys2Mon, *PC2);
	iterationCount2 = sys2Mon.iteration_count();
	if (!sys2Mon.converged())
	{
		std::cout << "ERROR: Solve for Lambda failed at time step " << timeStep << std::endl;
		std::cout << "Iterations   : " << iterationCount2 << std::endl;          
		std::cout << "Residual norm: " << sys2Mon.residual_norm() << std::endl;
		std::cout << "Tolerance    : " << sys2Mon.tolerance() << std::endl;
		std::exit(-1);
	}
	
	logger.stopTimer("solvePoisson");
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::projectionStep()
{
	logger.startTimer("projectionStep");

	cusp::wrapped::multiply(Q, lambda, temp1);
	cusp::wrapped::multiply(BN, temp1, q);
	cusp::blas::axpby(qStar, q, q, 1.0, -1.0);

	logger.stopTimer("projectionStep");
}

//##############################################################################
//                               OUTPUT
//##############################################################################

template <typename memoryType>
void NavierStokesSolver<memoryType>::writeData()
{
	logger.startTimer("output");
	
	int nsave = (*paramDB)["simulation"]["nsave"].get<int>();
	std::string folderName = (*paramDB)["inputs"]["folderName"].get<std::string>();
	if (timeStep % nsave == 0)
	{
		io::writeData(folderName, timeStep, q, lambda, *domInfo);
	}
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	calculateForce();
//	io::writeForce(folderName, timeStep*dt, forceX, forceY);
	forceFile << timeStep*dt << '\t' << forceX << '\t' << forceY << std::endl;
	iterationsFile << timeStep << '\t' << iterationCount1 << '\t' << iterationCount2 << std::endl;
//	io::writeIterations(iterationsFile, folderName, timeStep, iterationCount1, iterationCount2);
	
	logger.stopTimer("output");
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::calculateForce()
{
}

template <typename memoryType>
void NavierStokesSolver<memoryType>::shutDown()
{
	io::printTimingInfo(logger);
//	forceFile.close();
	iterationsFile.close();
}

template <typename memoryType>
NavierStokesSolver<memoryType>* NavierStokesSolver<memoryType>::createSolver(parameterDB &paramDB, domain &domInfo)
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
	solver->domInfo = &domInfo;
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
