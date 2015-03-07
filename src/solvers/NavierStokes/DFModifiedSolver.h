/**
* @file  DFModifiedSolver.h
* @brief Modified direct forcing immersed boundary method.
*/

#pragma once

#include "DirectForcingSolver.h"
/**
 * \class DirectForcingSolver
 * \brief A fully discrete formulation of the direct forcing method, with an 
 *        explicit pressure term in the intermediate velocity solve.
 *
 * It inherits from the base class DirectForcingSolver. The difference here 
 * is that the pressure at time step \a n is included in the equation to solve 
 * the intermediate velocity. The quantity \a phi that is solved in the Poisson 
 * equation is a pressure correction that is projected on to the intermediate 
 * velocity to satisfy the divergence-free constraint. The pressure at time 
 * step \a n+1 is calculated by adding \a phi to the pressure at time step \a n.
 */
template <typename memoryType>
class DFModifiedSolver : public DirectForcingSolver<memoryType>
{
protected:
	virtual void calculateExplicitLambdaTerms();
	virtual void projectionStep();
	
public:
	DFModifiedSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "DFModified";
	}
};