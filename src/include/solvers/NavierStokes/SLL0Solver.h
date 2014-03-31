#pragma once

#include <solvers/NavierStokes/SuLaiLinSolver.h>

/**
* @brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class SLL0Solver : public SuLaiLinSolver<memoryType>
{
protected:
	
	void assembleRHS1();
	void solveIntermediateVelocity();
	
public:
	std::string name()
	{
		return "SLL0";
	}
};
