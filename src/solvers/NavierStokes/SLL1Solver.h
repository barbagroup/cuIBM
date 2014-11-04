/***************************************************************************//**
* \file  SLL1Solver.h
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

#include "SuLaiLinSolver.h"

/********************//**
* \class SLL1Solver
* \brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class SLL1Solver : public SuLaiLinSolver<memoryType>
{
protected:
	/**
	* \brief
	*/
	virtual void solveIntermediateVelocity();
	
public:
	/**
	* \brief Constructor of the class \c SLL1Solver
	*/
	SLL1Solver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	/**
	* \brief Return the name of the solver
	*/
	virtual std::string name()
	{
		return "SLL1";
	}
};
