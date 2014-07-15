/***************************************************************************//**
* \file NSWithBody.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Generic Navier-Stokes solver in the presence of immersed boundaries
*/

#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>

/**
* \class NSWithBody
* \brief Generic Navier-Stokes solver in the presence of immersed boundaries.
*
* The class inherits from the class \c NavierStokesSolver.
*
*/
template <typename memoryType>
class NSWithBody : public NavierStokesSolver<memoryType>
{
protected:
	bodies<memoryType> B;		///< bodies in the flow

	real forceX,				///< force acting on each body in the x-direction
		 forceY;				///< force acting on each body in the y-direction
	
	std::ofstream forceFile;	///< file to write the forces

	/**
	* \brief To be documented
	*/
	virtual void calculateForce();
	
	/**
	* \brief To be documented
	*/
	void initialiseBodies();
	
	/**
	* \brief To be documented
	*/
	void updateBodies();

public:
	/**
	* \brief To be documented
	*/
	void writeCommon();
	
	/**
	* \brief To be documented
	*/
	void shutDown();
};
