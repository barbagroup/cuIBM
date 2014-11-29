/**
* @file  DiffusionSolver.h
* @brief Solves diffusion with an immersed boundary.
*/

#pragma once

#include "DirectForcingSolver.h"

template <typename memoryType>
class DiffusionSolver : public DirectForcingSolver<memoryType>
{
private:
	virtual void solvePoisson();
	virtual void projectionStep();
	
public:
	virtual void initialise();
	DiffusionSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "Diffusion";
	}
};
