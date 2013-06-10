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

#include <solvers/NavierStokes/SLL1Solver.h>
#include <sys/stat.h>
#include <cusp/io/matrix_market.h>

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
	relTol = 1e-5;
	
	SuLaiLinSolver<memoryType>::assembleRHS3();  // assemble rhs3 to solve for f
	cusp::io::write_matrix_market_file(SuLaiLinSolver<memoryType>::F, "F.mtx");
	cusp::io::write_matrix_market_file(SuLaiLinSolver<memoryType>::rhs3, "rhs3.mtx");
	
	cusp::default_monitor<real> sys3Mon(SuLaiLinSolver<memoryType>::rhs3, maxIters, relTol);
	//cusp::krylov::bicgstab(F, f, rhs3, sys3Mon, *PC3);
	cusp::krylov::cg(SuLaiLinSolver<memoryType>::F, SuLaiLinSolver<memoryType>::f, SuLaiLinSolver<memoryType>::rhs3, sys3Mon);//, *PC3);
	cusp::io::write_matrix_market_file(SuLaiLinSolver<memoryType>::f, "f.mtx");
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
	
	cusp::wrapped::multiply(SuLaiLinSolver<memoryType>::ET, SuLaiLinSolver<memoryType>::f, SuLaiLinSolver<memoryType>::temp3);
	cusp::wrapped::multiply(NavierStokesSolver<memoryType>::BN, SuLaiLinSolver<memoryType>::temp3, NavierStokesSolver<memoryType>::qStar);
	cusp::blas::axpby(SuLaiLinSolver<memoryType>::qTilde, NavierStokesSolver<memoryType>::qStar, NavierStokesSolver<memoryType>::qStar, 1.0, -1.0);

	NavierStokesSolver<memoryType>::logger.stopTimer("solveIntermediateVel");
}

template class SLL1Solver<host_memory>;
template class SLL1Solver<device_memory>;
