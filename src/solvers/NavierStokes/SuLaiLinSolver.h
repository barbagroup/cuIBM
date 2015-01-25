/***************************************************************************//**
 * \file  SuLaiLinSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c SuLaiLinSolver.
 */


#pragma once

#include "NSWithBody.h"


/**
 * \class SuLaiLinSolver
 * \brief Immersed boundary method from Su, Lai and Lin.
 */
template <typename memoryType>
class SuLaiLinSolver : public NSWithBody<memoryType>
{
protected:
	cusp::coo_matrix<int, real, memoryType>
	     E, ET, F;
	
	cusp::array1d<real, memoryType>
	     qTilde, qDagger, rhs3, f, temp3, velB; // remove qDagger later and make it a part of SLL0Solver
	
	preconditioner< cusp::coo_matrix<int, real, memoryType> > *PC3;
	
	void generateE();
	void generateF();
	
	virtual void updateSolverState();
	
	void generateVelB();
	void assembleRHS3();
	
	virtual void calculateForce();
	
public:
	virtual void initialise();
	virtual void writeData();

	/**
	 * \brief Returns the name of the solver as a string.
	 */
	virtual std::string name()
	{
		return "Su, Lai & Lin";
	}
};
