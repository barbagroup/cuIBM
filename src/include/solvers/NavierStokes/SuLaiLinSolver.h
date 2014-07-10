/***************************************************************************//**
* \file  SuLaiLinSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Solves the flow using the IB method described by Taira and Colonius (2007)
*
* <b>The immersed boundary method: a projection approach</b> \n
* Taira, K and Colonius, T \n
* Journal of Computational Physics \n
* Volume 225 Number 2 \n
* 2007
*
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>

/********************//**
* \class SuLaiLinSolver
* \brief Immersed boundary method described by Taira and Colonius (2007)
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
	
	void updateSolverState();
	
	void generateVelB();
	void assembleRHS3();
	
	void calculateForce();
	
public:
	void initialise();
	void writeData();
	std::string name()
	{
		return "Su, Lai & Lin";
	}
};
