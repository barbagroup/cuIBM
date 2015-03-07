/***************************************************************************//**
 * \file  DiffusionSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c DiffusionSolver. 
 */

#pragma once

#include "DirectForcingSolver.h"

/**
 * \class DiffusionSolver
 * \brief Solves unsteady diffusion equation with an immersed boundary 
 *        using the direct forcing method.
 *
 * Inherits from the class DirectForcingSolver. Sets the convection term to 
 * zero, and does not solve any Poisson equation. The intermediate velocity is 
 * the same as the velocity at the next time step.
 * Each component of velocity is solved independently, and do not influence each other.
 * To solve the heat equation for only one variable, set boundary conditions 
 * and look at the solution only for one component of velocity.
 */
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
