/***************************************************************************//**
 * \file domain.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c domain.
 */


#pragma once

#include "types.h"


/**
 * \class domain
 * \brief Stores information about the computational grid.
 */
class domain
{
public:
	int   nx, ///< number of cells in the x-direction
	      ny; ///< number of cells in the y-direction
	
	vecH  x,  ///< x-coordinates of the nodes
	      y,  ///< y-coordinates of the nodes
	      dx, ///< cell-widths in the x-direction
	      dy; ///< cell-widths in the y-direction
	
	vecD  xD,  ///< x-coordinates of the nodes stored on the device
	      yD,  ///< y-coordinates of the nodes stored on the device
	      dxD, ///< x- cell widths stored on the device
	      dyD; ///< y- cell widths stored on the device
	
	vecH  xu,  ///< x-coordinates where the x-components of velocity are evaluated
	      yu,  ///< y-coordinates where the x-components of velocity are evaluated
	      xv,  ///< x-coordinates where the y-components of velocity are evaluated
	      yv;  ///< y-coordinates where the y-components of velocity are evaluated
};
