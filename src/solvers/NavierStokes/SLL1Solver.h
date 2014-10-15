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

#pragma once

#include "SuLaiLinSolver.h"

/**
* @brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class SLL1Solver : public SuLaiLinSolver<memoryType>
{
protected:
	virtual void solveIntermediateVelocity();
	
public:
	SLL1Solver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "SLL1";
	}
};
