#include <solvers/NavierStokes/SLL0Solver.h>
#include <sys/stat.h>
#include <cusp/io/matrix_market.h>

//##############################################################################
//                          GENERATE VECTORS
//##############################################################################

template <typename memoryType>
void SLL0Solver<memoryType>::assembleRHS1()
{
	NavierStokesSolver<memoryType>::logger.startTimer("assembleRHS1");
	
	cusp::blas::axpby(NavierStokesSolver<memoryType>::rn, NavierStokesSolver<memoryType>::bc1, NavierStokesSolver<memoryType>::rhs1, 1.0, 1.0);
	cusp::multiply(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::lambda, NavierStokesSolver<memoryType>::temp1);
	cusp::blas::axpby(NavierStokesSolver<memoryType>::rhs1, NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rhs1, 1.0, -1.0);
	
	NavierStokesSolver<memoryType>::logger.stopTimer("assembleRHS1");
}

//##############################################################################
//                           LINEAR SOLVES
//##############################################################################

template <typename memoryType>
void SLL0Solver<memoryType>::solveIntermediateVelocity()
{
	NavierStokesSolver<memoryType>::logger.startTimer("solveIntermediateVel");
	
	// Solve for qTilde ========================================================
	
	// [1] A.qTilde = r1 
	//     r1 here includes the gradient of the pressure at the previous time step
	
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	
	int  maxIters = db["velocitySolve"]["maxIterations"].get<int>();
	real relTol   = db["velocitySolve"]["tolerance"].get<real>();

	cusp::default_monitor<real> sys1Mon(NavierStokesSolver<memoryType>::rhs1, maxIters, relTol);
	cusp::krylov::cg(NavierStokesSolver<memoryType>::A, SuLaiLinSolver<memoryType>::qTilde, NavierStokesSolver<memoryType>::rhs1, sys1Mon, *NavierStokesSolver<memoryType>::PC1);
	
	NavierStokesSolver<memoryType>::iterationCount1 = sys1Mon.iteration_count();
	if (!sys1Mon.converged())
	{
		std::cout << "ERROR: Solve for q~ failed at time step " << NavierStokesSolver<memoryType>::timeStep << std::endl;
		std::cout << "Iterations   : " << NavierStokesSolver<memoryType>::iterationCount1 << std::endl;          
		std::cout << "Residual norm: " << sys1Mon.residual_norm() << std::endl;
		std::cout << "Tolerance    : " << sys1Mon.tolerance() << std::endl;
		std::exit(-1);
	}
	
	// Solve for f =============================================================
	
	// [2] (E.BN.ET)f = E.qTilde - uBn+1
	
	maxIters = 10000;
	relTol   = 1e-5;
	
	SuLaiLinSolver<memoryType>::assembleRHS3();  // assemble rhs3 to solve for f
	
	cusp::default_monitor<real> sys3Mon(SuLaiLinSolver<memoryType>::rhs3, maxIters, relTol);
	//cusp::krylov::bicgstab(F, f, rhs3, sys3Mon, *PC3);
	cusp::krylov::cg(SuLaiLinSolver<memoryType>::F, SuLaiLinSolver<memoryType>::f, SuLaiLinSolver<memoryType>::rhs3, sys3Mon);//, *PC3);
	int iterationCount3 = sys3Mon.iteration_count();
	if (!sys3Mon.converged())
	{
		std::cout << "ERROR: Solve for f failed at time step " << NavierStokesSolver<memoryType>::timeStep << std::endl;
		std::cout << "Iterations   : " << iterationCount3 << std::endl;          
		std::cout << "Residual norm: " << sys3Mon.residual_norm() << std::endl;
		std::cout << "Tolerance    : " << sys3Mon.tolerance() << std::endl;
		std::exit(-1);
	}
	
	// Obtain qDagger ==========================================================
	
	// [3] qDagger = qTilde - BN.ET.f
	
	cusp::multiply(SuLaiLinSolver<memoryType>::ET, SuLaiLinSolver<memoryType>::f, NavierStokesSolver<memoryType>::temp1);
	cusp::multiply(NavierStokesSolver<memoryType>::BN, NavierStokesSolver<memoryType>::temp1, SuLaiLinSolver<memoryType>::qDagger);
	cusp::blas::axpby(SuLaiLinSolver<memoryType>::qTilde, SuLaiLinSolver<memoryType>::qDagger, SuLaiLinSolver<memoryType>::qDagger, 1.0, -1.0);
	
	// Obtain q* ===============================================================
	
	// [3] q* = qDagger + BN.G.Pn
	
	cusp::multiply(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::lambda, NavierStokesSolver<memoryType>::temp1);
	cusp::multiply(NavierStokesSolver<memoryType>::BN, NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::qStar);
	cusp::blas::axpby(SuLaiLinSolver<memoryType>::qDagger, NavierStokesSolver<memoryType>::qStar, NavierStokesSolver<memoryType>::qStar, 1.0, 1.0);

	NavierStokesSolver<memoryType>::logger.stopTimer("solveIntermediateVel");
}

template <typename memoryType>
SLL0Solver<memoryType>::SLL0Solver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template class SLL0Solver<host_memory>;
template class SLL0Solver<device_memory>;
