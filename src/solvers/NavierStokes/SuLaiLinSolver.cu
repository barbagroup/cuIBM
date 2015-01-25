/***************************************************************************//**
 * \file SuLaiLinSolver.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c SuLaiLinSolver.
 */


#include "SuLaiLinSolver.h"
#include <sys/stat.h>

/**
 * \brief To be documented.
 */
template <typename memoryType>
void SuLaiLinSolver<memoryType>::initialise()
{	
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	
	NSWithBody<memoryType>::initialiseBodies();
	int totalPoints  = NSWithBody<memoryType>::B.totalPoints; 
	
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");

	qTilde.resize(numUV);
	qDagger.resize(numUV);
	if(totalPoints > 0) // if bodies are present in the flow, create the following ararys
	{
		f.resize(2*totalPoints);
		velB.resize(2*totalPoints);
		rhs3.resize(2*totalPoints);
		temp3.resize(2*totalPoints);
		E.resize(2*totalPoints, numUV, 24*totalPoints);
	}

	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	
	NavierStokesSolver<memoryType>::assembleMatrices();  // Assemble matrices
	generateE();
	generateF();
	generateVelB();
	NavierStokesSolver<memoryType>::logger.startTimer("preconditioner3");
	PC3 = new preconditioner< cusp::coo_matrix<int, real, memoryType> >(F, DIAGONAL);
	NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner3");
}

//##############################################################################
//                            TIME STEPPING
//##############################################################################

/**
 * \brief To be documented.
 */
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

/**
 * \brief To be documented.
 */
template <>
void SuLaiLinSolver<device_memory>::generateF()
{
	logger.startTimer("generateF");
	
	cooD temp; // Should this temp matrix be created each time step?
	cusp::multiply(E, BN, temp);
	cusp::multiply(temp, ET, F);
	
	logger.stopTimer("generateF");
}

/**
 * \brief To be documented.
 */
template <>
void SuLaiLinSolver<host_memory>::generateF()
{
	logger.startTimer("generateF");
	
	cooH temp;
	cusp::multiply(E, BN, temp);
	cusp::multiply(temp, ET, F);
	F.sort_by_row_and_column();
	
	logger.stopTimer("generateF");
}

//##############################################################################
//                          GENERATE VECTORS
//##############################################################################

/**
 * \brief To be documented.
 */
template <typename memoryType>
void SuLaiLinSolver<memoryType>::assembleRHS3()
{
	NavierStokesSolver<memoryType>::logger.startTimer("assembleRHS3");
	
	// rhs3 = (Eq~ - uBn+1)
	cusp::multiply(E, qTilde, temp3);
	cusp::blas::axpby(temp3, velB, rhs3, 1.0, -1.0);
	
	NavierStokesSolver<memoryType>::logger.stopTimer("assembleRHS3");
}

//##############################################################################
//                               OUTPUT
//##############################################################################

/**
 * \brief To be documented.
 */
template <typename memoryType>
void SuLaiLinSolver<memoryType>::writeData()
{
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	NSWithBody<memoryType>::writeCommon();

	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	real dt = db["simulation"]["dt"].get<real>();

	// Print forces calculated using both the SLL method and the CV approach
	calculateForce();
	NSWithBody<memoryType>::forceFile << NavierStokesSolver<memoryType>::timeStep*dt << '\t' << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << '\t';
	NSWithBody<memoryType>::calculateForce();
	NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << std::endl;
	
	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}

// include inline files located in "SuLaiLin"
#include "SuLaiLin/generateE.inl"
#include "SuLaiLin/generateVelB.inl"
#include "SuLaiLin/calculateForce.inl"

// specialization of the class SuLaiLinSolver
template class SuLaiLinSolver<host_memory>;
template class SuLaiLinSolver<device_memory>;
