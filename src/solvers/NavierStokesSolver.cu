/**
 * \file NavierStokesSolver.cu
 * \brief Implementation of the methods of the class \c NavierStokesSolver.
 */


#include <sys/stat.h>

#include "NavierStokesSolver.h"
#include "io/io.h"


//##############################################################################
//                              INITIALISE
//##############################################################################

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 *
 * \param pDB database that contains all the simulation parameters
 * \param dInfo information related to the computational grid
 */
template <typename memoryType>
NavierStokesSolver<memoryType>::NavierStokesSolver(parameterDB *pDB, domain *dInfo)
{
	paramDB = pDB;
	domInfo = dInfo;
} // NavierStokesSolver


/**
 * \brief Initializes parameters, arrays and matrices required for the simulation.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::initialise()
{
	printf("Initializing Navier-Stokes solver ...\n");
	
	int nx = domInfo->nx,
	    ny = domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1),
	    numP = nx*ny;
	
	initialiseCommon();
	initialiseArrays(numUV, numP);
	assembleMatrices();
} // initialise


/**
 * \brief Initializes parameters common to all Navier-Stokes solvers.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseCommon()
{
	printf("Initializing common parts ...\n");
	logger.startTimer("initialiseCommon");
	
	QCoeff = 1.0;
	subStep = 0;
	
	timeScheme convScheme = (*paramDB)["simulation"]["convTimeScheme"].get<timeScheme>(),
	           diffScheme = (*paramDB)["simulation"]["diffTimeScheme"].get<timeScheme>();
	intgSchm.initialise(convScheme, diffScheme);
	
	// set initial timeStep
	timeStep = (*paramDB)["simulation"]["startStep"].get<int>();
	// get folder path 
	std::string folder = (*paramDB)["inputs"]["caseFolder"].get<std::string>();

	// writes the grids information to a file
	io::writeGrid(folder, *domInfo);

	// opens the file to which the number of iterations at every step is written
	std::stringstream out;
	out << folder << "/iterations";
	if (timeStep == 0)
	{
		iterationsFile.open(out.str().c_str());
	}
	else
	{
		iterationsFile.open(out.str().c_str(), std::ofstream::app);
	}
	
	logger.stopTimer("initialiseCommon");
} // initialiseCommon


/**
 * \brief Initializes all arrays required to solve the Navier-Stokes equations.
 *
 * \param numQ total number velocity (or flux) unknowns (x- and y- directions)
 * \param numLambda number of pressure unknowns (plus number of body force unknowns)
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseArrays(int numQ, int numLambda)
{	
	printf("Initializing arrays ...\n");
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
	
	generateRN();
	cusp::blas::scal(H, 1.0/intgSchm.gamma[subStep]);
	
	logger.stopTimer("initialiseArrays");
} // initialiseArrays


/**
 * \brief Initializes velocity flux vectors (on the host).
 *
 * It creates a raw pointer before calling a method to initialize the flux vector.
 *
 */
template <>
void NavierStokesSolver<host_memory>::initialiseFluxes()
{
	real *q_r = thrust::raw_pointer_cast(&(q[0]));
	initialiseFluxes(q_r);
	qStar = q;
} // initialiseFluxes


/**
 * \brief Initializes velocity flux vectors (on the device).
 *
 * It creates a raw pointer before calling a method to initialize the flux vector.
 *
 */
template<>
void NavierStokesSolver <device_memory>::initialiseFluxes()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	vecH qHost((nx-1)*ny+nx*(ny-1));
	
	// creating raw pointers
	real *qHost_r = thrust::raw_pointer_cast(&(qHost[0]));
	initialiseFluxes(qHost_r);
	q = qHost;
	qStar=q;
} // initialiseFluxes


/**
 * \brief Initializes velocity flux vectors.
 *
 * \param q the velocity flux vector
 */
template <typename memoryType>
void NavierStokesSolver <memoryType>::initialiseFluxes(real *q)
{
	if (timeStep != 0)
	{
		// case directory
		std::string caseFolder = (*paramDB)["inputs"]["caseFolder"].get<std::string>();
		// read velocity fluxes from file
		io::readData(caseFolder, timeStep, q, "q");
		return;
	}

	int  nx = domInfo->nx,
	     ny = domInfo->ny,
	     numU = (nx-1)*ny;
	
	real xmin = domInfo->x[0],
	     xmax = domInfo->x[nx],
	     ymin = domInfo->y[0],
	     ymax = domInfo->y[ny];
	
	real uInitial = (*paramDB)["flow"]["uInitial"].get<real>(),
	     uPerturb = (*paramDB)["flow"]["uPerturb"].get<real>(),
	     vInitial = (*paramDB)["flow"]["vInitial"].get<real>(),
	     vPerturb = (*paramDB)["flow"]["vPerturb"].get<real>();
	
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			q[j*(nx-1) + i] = ( uInitial + uPerturb * cos( 0.5*M_PI*(2*domInfo->xu[i]-xmax-xmin)/(xmax-xmin) ) * sin( M_PI * (2*domInfo->yu[j]-ymax-ymin)/(ymax-ymin) ) ) * domInfo->dy[j];
		}
	}
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			q[j*nx + i + numU] = ( vInitial + vPerturb * cos( 0.5*M_PI*(2*domInfo->yv[j]-ymax-ymin)/(ymax-ymin) ) * sin( M_PI * (2*domInfo->xv[i]-xmax-xmin)/(xmax-xmin) ) ) * domInfo->dx[i];
		}
	}
} // initialiseFluxes


/**
 * \brief Initializes boundary velocity arrays with values stored in the database.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::initialiseBoundaryArrays()
{
	int nx = domInfo->nx,
	    ny = domInfo->ny;

	boundaryCondition **bcInfo = (*paramDB)["flow"]["boundaryConditions"].get<boundaryCondition **>();

	// resize boundary arrays by the number of velocity points on boundaries (u and v points)
	bc[XMINUS].resize(2*ny-1);
	bc[XPLUS].resize(2*ny-1);
	bc[YMINUS].resize(2*nx-1);
	bc[YPLUS].resize(2*nx-1);

	/// Top and Bottom
	for(int i=0; i<nx-1; i++)
	{
		bc[YMINUS][i] = bcInfo[YMINUS][0].value;
		bc[YPLUS][i] = bcInfo[YPLUS][0].value;
		bc[YMINUS][i+nx-1] = bcInfo[YMINUS][1].value;
		bc[YPLUS][i+nx-1]	= bcInfo[YPLUS][1].value;
	}
	bc[YMINUS][2*nx-2] = bcInfo[YMINUS][1].value;
	bc[YPLUS][2*nx-2]	= bcInfo[YPLUS][1].value;
	
	/// Left and Right
	for(int i=0; i<ny-1; i++)
	{
		bc[XMINUS][i] = bcInfo[XMINUS][0].value;
		bc[XPLUS][i] = bcInfo[XPLUS][0].value;
		bc[XMINUS][i+ny] = bcInfo[XMINUS][1].value;
		bc[XPLUS][i+ny] = bcInfo[XPLUS][1].value;
	}
	bc[XMINUS][ny-1] = bcInfo[XMINUS][0].value;
	bc[XPLUS][ny-1] = bcInfo[XPLUS][0].value;
} // initialiseBoundaryArrays


//##############################################################################
//                            TIME STEPPING
//##############################################################################

/**
 * \brief Calculates the variables at the next time step.
 */
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
} // stepTime


/**
 * \brief Doing nothing. Used in immersed boundary methods when the body moves.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::updateSolverState()
{
//	generateA(intgSchm.alphaImplicit[subStep]);
//	updateQ(intgSchm.gamma[i]);
//	updateBoundaryConditions();
} // updateSolverState


/**
 * \brief Doing nothing.
 *
 * \param gamma coefficient of the convection term at the current time step
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::updateQ(real gamma)
{
//	cusp::blas::scal(Q.values, gamma/QCoeff);
//	QCoeff = gamma;
} // updateQ


/**
 * \brief Doing nothing.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::updateBoundaryConditions()
{
} // updateBoundaryConditions


/**
 * \brief Evaluates the condition required to stop the simulation.
 *
 * \return a Boolean to continue or stop the simulation
 */
template <typename memoryType>
bool NavierStokesSolver<memoryType>::finished()
{
	int startStep = (*paramDB)["simulation"]["startStep"].get<int>();
	int nt = (*paramDB)["simulation"]["nt"].get<int>();
	return (timeStep < startStep+nt) ? false : true;
} // finished


//##############################################################################
//                          ASSEMBLE MATRICES
//##############################################################################

/**
 * \brief Assembles matrices of the intermediate flux solver and the Poisson solver.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleMatrices()
{
	printf("Initializing matrices ...\n");
	logger.startTimer("assembleMatrices");
	
	generateM();
	generateL();
	generateA(intgSchm.alphaImplicit[subStep]);
	PC1 = new preconditioner< cusp::coo_matrix<int, real, memoryType> >(A, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>());
	generateBN();	
	
	logger.stopTimer("assembleMatrices");

	generateQT();
	generateC(); // QT*BN*Q
	
	logger.startTimer("preconditioner2");
	PC2 = new preconditioner< cusp::coo_matrix<int, real, memoryType> >(C, (*paramDB)["PoissonSolve"]["preconditioner"].get<preconditionerType>());
	logger.stopTimer("preconditioner2");
} // assembleMatrices


/**
 * \brief Generates approximate inverse of the matrix resulting from implicit velocity terms.
 *
 * It computes the N-th order Taylor expansion of the inverse matrix.
 * Currently, the order is N=1.
 *
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::generateBN()
{
	BN = Minv; // 1st-order
} // generateBN


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


/**
 * \brief Generates the matrix of the Poisson solver (on the device).
 */
template <>
void NavierStokesSolver<device_memory>::generateC()
{
	logger.startTimer("generateC");
	
	cooD temp; // Should this temp matrix be created each time step?
	cusp::multiply(QT, BN, temp);
	cusp::multiply(temp, Q, C);
	C.values[0] += C.values[0];
	
	logger.stopTimer("generateC");
} // generateC


/**
 * \brief Generates the matrix of the Poisson solver (on the host).
 */
template <>
void NavierStokesSolver<host_memory>::generateC()
{
	logger.startTimer("generateC");
	
	cooH temp;
	cusp::multiply(QT, BN, temp);
	cusp::multiply(temp, Q, C);
	C.sort_by_row_and_column();
	C.values[0] += C.values[0];
	
	logger.stopTimer("generateC");
} // generateC


//##############################################################################
//                          GENERATE VECTORS
//##############################################################################

/**
 * \brief Assembles the right hand-side of the system for the intermediate flux.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleRHS1()
{
	logger.startTimer("assembleRHS1");
	
	cusp::blas::axpby(rn, bc1, rhs1, 1.0, 1.0);
	
	logger.stopTimer("assembleRHS1");
} // assembleRHS1


/**
 * \brief Assembles the right hand-side of the Poisson system.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::assembleRHS2()
{
	logger.startTimer("assembleRHS2");
	
	cusp::multiply(QT, qStar, temp2);
	cusp::blas::axpby(temp2, bc2, rhs2, 1.0, -1.0);
	
	logger.stopTimer("assembleRHS2");
} // assembleRHS2


//##############################################################################
//                           LINEAR SOLVES
//##############################################################################

/**
 * \brief Solves for the intermediate flux velocity.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::solveIntermediateVelocity()
{
	logger.startTimer("solveIntermediateVel");
	
	int maxIte = (*paramDB)["velocitySolve"]["maxIterations"].get<int>();
	real rTol = (*paramDB)["velocitySolve"]["rTol"].get<real>();
	real aTol = (*paramDB)["velocitySolve"]["aTol"].get<real>();

	cusp::monitor<real> sys1Mon(rhs1, maxIte, rTol, aTol, false);
	cusp::krylov::bicgstab(A, qStar, rhs1, sys1Mon, *PC1);
	//cusp::krylov::cg(A, qStar, rhs1, sys1Mon, *PC1);
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
} // solveIntermediateVelocity


/**
 * \brief Solves the Poisson system.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::solvePoisson()
{
	logger.startTimer("solvePoisson");
	
	int  maxIte = (*paramDB)["PoissonSolve"]["maxIterations"].get<int>();
	real rTol = (*paramDB)["PoissonSolve"]["rTol"].get<real>();
	real aTol = (*paramDB)["PoissonSolve"]["aTol"].get<real>();
	
	cusp::monitor<real> sys2Mon(rhs2, maxIte, rTol, aTol, false);
	if((*paramDB)["simulation"]["ibmScheme"].get<ibmScheme>() != DF_IMPROVED)
	{
		cusp::krylov::cg(C, lambda, rhs2, sys2Mon, *PC2);
	}
	else
	{
		cusp::krylov::bicgstab(C, lambda, rhs2, sys2Mon, *PC2);
		//cusp::krylov::gmres(C, lambda, rhs2, 50, sys2Mon, *PC2);
	}
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
} // solvePoisson


/**
 * \brief Projects the flux onto the divergence-free field.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::projectionStep()
{
	logger.startTimer("projectionStep");

	cusp::multiply(Q, lambda, temp1);
	cusp::multiply(BN, temp1, q);
	cusp::blas::axpby(qStar, q, q, 1.0, -1.0);

	logger.stopTimer("projectionStep");
} // projectionStep


//##############################################################################
//                               OUTPUT
//##############################################################################

/**
 * \brief Writes numerical solution at current time-step,
 *        as well as the number of iterations performed in each solver.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::writeCommon()
{
	int nsave = (*paramDB)["simulation"]["nsave"].get<int>();
	std::string folder = (*paramDB)["inputs"]["caseFolder"].get<std::string>();
	
	// write the velocity fluxes and the pressure values
	if (timeStep % nsave == 0)
	{
		io::writeData(folder, timeStep, q, lambda, *domInfo);
	}
	
	// write the number of iterations for each solve
	iterationsFile << timeStep << '\t' << iterationCount1 << '\t' << iterationCount2 << std::endl;
} // writeCommon


/**
 * \brief Writes data into files.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::writeData()
{
	logger.startTimer("output");

	writeCommon();	
	
	logger.stopTimer("output");
} // writeData


/**
 * \brief Prints timing information and closes the different files.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::shutDown()
{
	io::printTimingInfo(logger);
	iterationsFile.close();
} // shutDown


// include inline files
#include "NavierStokes/generateM.inl"
#include "NavierStokes/generateL.inl"
#include "NavierStokes/generateA.inl"
#include "NavierStokes/generateQT.inl"
#include "NavierStokes/generateRN.inl"
#include "NavierStokes/generateBC1.inl"
#include "NavierStokes/generateBC2.inl"


// specialization of the class NavierStokesSolver
template class NavierStokesSolver<host_memory>;
template class NavierStokesSolver<device_memory>;
