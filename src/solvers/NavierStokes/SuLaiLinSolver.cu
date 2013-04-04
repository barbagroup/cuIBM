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

#include <solvers/NavierStokes/SuLaiLinSolver.h>
#include <sys/stat.h>
#include <cusp/io/matrix_market.h>

template <typename memoryType>
void SuLaiLinSolver<memoryType>::initialise()
{	
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	
	NSWithBody<memoryType>::initialiseBodies();
	int numB  = NSWithBody<memoryType>::B.totalPoints; 
	
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	qTilde.resize(numUV);
	if(numB > 0)
	{
		f.resize(2*numB);
		velB.resize(2*numB);
		rhs3.resize(2*numB);
		temp3.resize(2*numB);
		E.resize(2*numB, numUV, 24*numB);
	}
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	
	// Assemble matrices
	NavierStokesSolver<memoryType>::assembleMatrices();
	generateE();
	generateF();
	cusp::io::write_matrix_market_file(F, "F.mtx");
	cusp::io::write_matrix_market_file(E, "E.mtx");
	NavierStokesSolver<memoryType>::logger.startTimer("preconditioner3");
	PC3 = new preconditioner< cusp::coo_matrix<int, real, memoryType> >(F, DIAGONAL);
	NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner3");
}

template <typename memoryType>
void SuLaiLinSolver<memoryType>::updateSolverState()
{
	if (NSWithBody<memoryType>::B.bodiesMove)
	{
		NSWithBody<memoryType>::updateBodies();
		generateE();
		generateF();
		
		NavierStokesSolver<memoryType>::logger.startTimer("preconditioner3");
		PC3->update(F);
		NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner3");
	}
}

//##############################################################################
//                          ASSEMBLE MATRICES
//##############################################################################

template <>
void SuLaiLinSolver<device_memory>::generateF()
{
	logger.startTimer("generateF");
	
	cooD temp; // Should this temp matrix be created each time step?
	cusp::wrapped::multiply(E, BN, temp);
	cusp::wrapped::multiply(temp, ET, F);
	
	logger.stopTimer("generateF");
}

template <>
void SuLaiLinSolver<host_memory>::generateF()
{
	logger.startTimer("generateF");
	
	cooH temp;
	cusp::wrapped::multiply(E, BN, temp);
	cusp::wrapped::multiply(temp, ET, F);
	F.sort_by_row_and_column();
	
	logger.stopTimer("generateF");
}

//##############################################################################
//                          GENERATE VECTORS
//##############################################################################

template <typename memoryType>
void SuLaiLinSolver<memoryType>::assembleRHS3()
{
	NavierStokesSolver<memoryType>::logger.startTimer("assembleRHS3");
	
	cusp::wrapped::multiply(E, qTilde, temp3);
	cusp::blas::axpby(temp3, velB, rhs3, 1.0, -1.0);
	
	NavierStokesSolver<memoryType>::logger.stopTimer("assembleRHS3");
}

//##############################################################################
//                           LINEAR SOLVES
//##############################################################################

//NavierStokesSolver<memoryType>::

/// MODIFY 

template <typename memoryType>
void SuLaiLinSolver<memoryType>::solveIntermediateVelocity()
{
	NavierStokesSolver<memoryType>::logger.startTimer("solveIntermediateVel");
	
	// Solve for qTilde ========================================================
	
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	
	int  maxIters = db["velocitySolve"]["maxIterations"].get<int>();
	real relTol = db["velocitySolve"]["tolerance"].get<real>();

	cusp::default_monitor<real> sys1Mon(NavierStokesSolver<memoryType>::rhs1, maxIters, relTol);
	cusp::krylov::cg(NavierStokesSolver<memoryType>::A, qTilde, NavierStokesSolver<memoryType>::rhs1, sys1Mon, *NavierStokesSolver<memoryType>::PC1);
	
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
	
	assembleRHS3();
	//cusp::print(rhs3);
	
	maxIters = 10000;
	relTol = 1e-5;
	
	cusp::io::write_matrix_market_file(rhs3, "rhs3.mtx");
	
	cusp::default_monitor<real> sys3Mon(rhs3, maxIters, relTol);
	cusp::krylov::cg(F, f, rhs3, sys3Mon, *PC3);
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
	
	cusp::wrapped::multiply(ET, f, temp3);
	cusp::wrapped::multiply(NavierStokesSolver<memoryType>::BN, temp3, NavierStokesSolver<memoryType>::qStar);
	cusp::blas::axpby(qTilde, NavierStokesSolver<memoryType>::qStar, NavierStokesSolver<memoryType>::qStar, 1.0, -1.0);

	NavierStokesSolver<memoryType>::logger.stopTimer("solveIntermediateVel");
}

/*template <typename memoryType>
void SuLaiLinSolver<memoryType>::intermediateProjectionStep()
{
	NavierStokesSolver<memoryType>::logger.startTimer("intermediateProjectionStep");

	cusp::wrapped::multiply(ET, f, temp3);
	cusp::wrapped::multiply(NavierStokesSolver<memoryType>::BN, temp3, NavierStokesSolver<memoryType>::qStar);
	cusp::blas::axpby(qTilde, NavierStokesSolver<memoryType>::qStar, NavierStokesSolver<memoryType>::qStar, 1.0, 1.0);

	NavierStokesSolver<memoryType>::logger.stopTimer("intermediateProjectionStep");
}*/

//##############################################################################
//                               OUTPUT
//##############################################################################

template <typename memoryType>
void SuLaiLinSolver<memoryType>::writeData()
{
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	NavierStokesSolver<memoryType>::writeCommon();

	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	real dt = db["simulation"]["dt"].get<real>();

	// Print forces calculated using both the T&C method and the CV approach
	calculateForce();
	NSWithBody<memoryType>::forceFile << NavierStokesSolver<memoryType>::timeStep*dt << '\t' << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << '\t';
	NSWithBody<memoryType>::calculateForce();
	NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << std::endl;
	
	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}

#include "SuLaiLin/generateE.inl"
#include "SuLaiLin/generateVelB.inl"
#include "SuLaiLin/calculateForce.inl"

template class SuLaiLinSolver<host_memory>;
template class SuLaiLinSolver<device_memory>;
