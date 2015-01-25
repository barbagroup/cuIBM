/***************************************************************************//**
 * \file NSWithBody.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c NSWithBody.
 */


#pragma once

#include "NavierStokesSolver.h"
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

	// calculate forces acting on an immersed body
	virtual void calculateForce();
	
	// store parameters of simulation and initialize immersed bodies
	void initialiseBodies();
	
	// update location and motion of each immersed body at current time
	void updateBodies();

public:
	// write flow variables and position of body points into files
	virtual void writeCommon();
	
	// close iteration file and force file
	virtual void shutDown();
};
