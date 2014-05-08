/**
* @file  SLL2Solver.h
* @brief Solves the flow using the IB method described by Taira and Colonius (2007)
*
* <b>The immersed boundary method: a projection approach</b> \n
* Taira, K and Colonius, T \n
* Journal of Computational Physics \n
* Volume 225 Number 2 \n
* 2007
*/

#pragma once

#include <solvers/NavierStokes/SuLaiLinSolver.h>

/**
* @brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class SLL2Solver : public SuLaiLinSolver<memoryType>
{
protected:
	void projectionStep();

public:
	SLL2Solver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	std::string name()
	{
		return "SLL2";
	}
};
