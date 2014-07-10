/***************************************************************************//**
* \file  NSWithBody.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Generic Navier-Stokes solver in the presence of immersed boundaries
*/

#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>

/********************//**
* \class NSWithBody
* \brief Generic Navier-Stokes solver in the presence of immersed boundaries
*/
template <typename memoryType>
class NSWithBody : public NavierStokesSolver<memoryType>
{
protected:
	bodies<memoryType> B;		///< bodies in the flow

	real forceX,				///< force acting on the body in the x-direction
		 forceY;				///< force acting on the body in the y-direction
	
	std::ofstream forceFile;	///< file to write the forces

	/**
	* \brief
	*/
	virtual void calculateForce();
	
	/**
	* \brief
	*/
	void initialiseBodies();
	
	/**
	* \brief
	*/
	void updateBodies();

public:
	/**
	* \brief
	*/
	void writeCommon();
	
	/**
	* \brief
	*/
	void shutDown();
};
