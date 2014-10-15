#pragma once

#include "SuLaiLinSolver.h"

/**
* @brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class SLL0Solver : public SuLaiLinSolver<memoryType>
{
protected:
	
	virtual void assembleRHS1();
	virtual void solveIntermediateVelocity();
	
public:
	SLL0Solver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "SLL0";
	}
};
