#include "SLL1Solver.h"
#include <sys/stat.h>

//##############################################################################
//                           LINEAR SOLVES
//##############################################################################

template <typename memoryType>
void SLL1Solver<memoryType>::solveIntermediateVelocity()
{
	NavierStokesSolver<memoryType>::logger.startTimer("solveIntermediateVel");
	
	// Solve for qTilde ========================================================
	
	// [1] A.qTilde = r1
	
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
	
	// Obtain q* ===============================================================
	
	// [3] q* = qTilde - BN.ET.f
	
	cusp::multiply(SuLaiLinSolver<memoryType>::ET, SuLaiLinSolver<memoryType>::f, NavierStokesSolver<memoryType>::temp1);
	cusp::multiply(NavierStokesSolver<memoryType>::BN, NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::qStar);
	cusp::blas::axpby(SuLaiLinSolver<memoryType>::qTilde, NavierStokesSolver<memoryType>::qStar, NavierStokesSolver<memoryType>::qStar, 1.0, -1.0);

	NavierStokesSolver<memoryType>::logger.stopTimer("solveIntermediateVel");
}

template <typename memoryType>
SLL1Solver<memoryType>::SLL1Solver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template class SLL1Solver<host_memory>;
template class SLL1Solver<device_memory>;
