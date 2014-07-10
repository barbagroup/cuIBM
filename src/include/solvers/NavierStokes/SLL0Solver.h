/***************************************************************************//**
* \file SSL0Solver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief
*/

#pragma once

#include <solvers/NavierStokes/SuLaiLinSolver.h>

/********************//**
* \class SLL0Solver
* \brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class SLL0Solver : public SuLaiLinSolver<memoryType>
{
protected:
	/**
	* \brief
	*/
	void assembleRHS1();

	/**
	* \brief
	*/
	void solveIntermediateVelocity();
	
public:
	/**
	* \brief Constructor of the class \c SLL0Solver
	*/
	SLL0Solver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	/**
	* \brief Return the name of the solver
	*/
	std::string name()
	{
		return "SLL0";
	}
};
