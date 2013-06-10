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

/**
* @file  SuLaiLinSolver.h
* @brief Solves the flow using the IB method described by Taira and Colonius (2007)
*
* <b>The immersed boundary method: a projection approach</b> \n
* Taira, K and Colonius, T \n
* Journal of Computational Physics \n
* Volume 225 Number 2 \n
* 2007
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
	if(numB > 0) // if bodies are present in the flow, create the following ararys
	{
		f.resize(2*numB);
		velB.resize(2*numB);
		rhs3.resize(2*numB);
		temp3.resize(2*numB);
		E.resize(2*numB, numUV, 24*numB);
	}

	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	
	NavierStokesSolver<memoryType>::assembleMatrices();  // Assemble matrices
	generateE();
	generateF();
	generateVelB();
	cusp::io::write_matrix_market_file(E, "E.mtx");
	cusp::io::write_matrix_market_file(F, "F.mtx");
	NavierStokesSolver<memoryType>::logger.startTimer("preconditioner3");
	PC3 = new preconditioner< cusp::coo_matrix<int, real, memoryType> >(F, DIAGONAL);
	NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner3");
}

//##############################################################################
//                            TIME STEPPING
//##############################################################################

template <typename memoryType>
void SuLaiLinSolver<memoryType>::updateSolverState()
{
	if (NSWithBody<memoryType>::B.bodiesMove)
	{
		NSWithBody<memoryType>::updateBodies();  // move the bodies
		generateE();                             // generate the E matrix (normalised interpolation matrix)
		generateF();                             // perform the matrix-matrix-matrix multiplication: F = E.BN.ET
		generateVelB();                          // generate the vector containing the body velocity at the boundary points
	
		NavierStokesSolver<memoryType>::logger.startTimer("preconditioner3");
	
		PC3->update(F);  // update the preconditioner for matrix F
	
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
	
	// rhs3 = (Eq~ - uBn+1)
	cusp::wrapped::multiply(E, qTilde, temp3);
	cusp::blas::axpby(temp3, velB, rhs3, 1.0, -1.0);
	
	NavierStokesSolver<memoryType>::logger.stopTimer("assembleRHS3");
}

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

	// Print forces calculated using both the SLL method and the CV approach
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
