/***************************************************************************//**
* \file domain.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the class \c domain
*/


#pragma once

#include "types.h"


/**
* \class domain
* \brief Store the mesh grid information
*/
class domain
{
public:
	int   nx, ///< number of cells in the x-direction
	      ny; ///< number of cells in the y-direction
	
	vecH  x,  ///< x-coordinates of the nodes
	      y,  ///< y-coordinates of the nodes
	      dx, ///< cell widths in the x-direction
	      dy; ///< cell widths in the y-direction
	
	vecD  xD,  ///< x-coordinates of the nodes stored on the device
	      yD,  ///< y-coordinates of the nodes stored on the device
	      dxD, ///< x- cell widths stored on the device
	      dyD; ///< y- cell widths stored on the device
	
	vecH  xu,  ///< x-coordinates of the locations at which the x-component of velocity is evaluated
	      yu,  ///< y-coordinates of the locations at which the x-component of velocity is evaluated
	      xv,  ///< x-coordinates of the locations at which the y-component of velocity is evaluated
	      yv;  ///< y-coordinates of the locations at which the y-component of velocity is evaluated
};