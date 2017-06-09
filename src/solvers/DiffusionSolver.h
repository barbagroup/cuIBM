/**
 * \file  DiffusionSolver.h
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
  // do nothing -- no Poisson system needs to be solved
	virtual void solvePoisson();

  // project velocity onto divergence-free space
	virtual void projectionStep();
	
public:
  // constructor -- get simulation parameters and grid
  DiffusionSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);

  // initialize solver
	virtual void initialise();

	/**
   * \brief Return name of solver.
   */
	virtual std::string name()
	{
		return "Diffusion";
	}
}; // DiffusionSolver
