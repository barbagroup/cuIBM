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

	/**
	 * \brief Calculates drag and lift for each immersed body.
	 */
	virtual void calculateForce();
	
	/**
	 * \brief Initializes immersed bodies.
	 */
	void initialiseBodies();
	
	/**
	 * \brief Updates location and velocity of immersed bodies at current time.
	 */
	void updateBodies();

public:
	/**
	 * \brief Writes flow variables and location of body points into files.
	 */
	virtual void writeCommon();
	
	/**
	 * \brief Closes iteration file and force file.
	 */
	virtual void shutDown();
};
