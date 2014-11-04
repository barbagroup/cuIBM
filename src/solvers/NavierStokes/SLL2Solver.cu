#include "SLL2Solver.h"
#include <sys/stat.h>

//##############################################################################
//                           LINEAR SOLVES
//##############################################################################

template <typename memoryType>
void SLL2Solver<memoryType>::projectionStep()
{
	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");
	
	// Calculate qTilde ========================================================
	
	cusp::multiply(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::lambda, NavierStokesSolver<memoryType>::temp1);
	cusp::multiply(NavierStokesSolver<memoryType>::BN, NavierStokesSolver<memoryType>::temp1, SuLaiLinSolver<memoryType>::qTilde);
	cusp::blas::axpby(NavierStokesSolver<memoryType>::qStar, SuLaiLinSolver<memoryType>::qTilde, SuLaiLinSolver<memoryType>::qTilde, 1.0, -1.0);
	
	// Solve for f =============================================================
	
	SuLaiLinSolver<memoryType>::assembleRHS3();  // assemble rhs3 to solve for f
	
	int  maxIters = 10000;
	real relTol = 1e-5;
	
	cusp::default_monitor<real> sys3Mon(SuLaiLinSolver<memoryType>::rhs3, maxIters, relTol);
	//cusp::krylov::cg(F, f, rhs3, sys3Mon, *PC3);
	cusp::krylov::bicgstab(SuLaiLinSolver<memoryType>::F, SuLaiLinSolver<memoryType>::f, SuLaiLinSolver<memoryType>::rhs3, sys3Mon);//, *PC3);
	int iterationCount3 = sys3Mon.iteration_count();
	if (!sys3Mon.converged())
	{
		std::cout << "ERROR: Solve for f failed at time step " << NavierStokesSolver<memoryType>::timeStep << std::endl;
		std::cout << "Iterations   : " << iterationCount3 << std::endl;          
		std::cout << "Residual norm: " << sys3Mon.residual_norm() << std::endl;
		std::cout << "Tolerance    : " << sys3Mon.tolerance() << std::endl;
		std::exit(-1);
	}
	
	// Obtain q^n+1 ===============================================================
	
	cusp::multiply(SuLaiLinSolver<memoryType>::ET, SuLaiLinSolver<memoryType>::f, NavierStokesSolver<memoryType>::temp1);
	cusp::multiply(NavierStokesSolver<memoryType>::BN, SuLaiLinSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::q);
	cusp::blas::axpby(SuLaiLinSolver<memoryType>::qTilde, NavierStokesSolver<memoryType>::q, NavierStokesSolver<memoryType>::q, 1.0, -1.0);

	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
	
}

template <typename memoryType>
SLL2Solver<memoryType>::SLL2Solver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template class SLL2Solver<host_memory>;
template class SLL2Solver<device_memory>;
