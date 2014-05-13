/**
* @file  NSWithBody.h
* @brief Generic Navier-Stokes solver in the presence of immersed boundaries
*/

#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>

/**
* @brief Generic Navier-Stokes solver in the presence of immersed boundaries
*/
template <typename memoryType>
class NSWithBody : public NavierStokesSolver<memoryType>
{
protected:
	bodies<memoryType> B; ///< Bodies in the flow

	real forceX, forceY;
	
	std::ofstream forceFile;

	virtual void calculateForce();
	void initialiseBodies();
	void updateBodies();

public:
	void writeCommon();
	void shutDown();
};
