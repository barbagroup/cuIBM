/**
* \file domain.h
* \brief Stores the grid information
*/
#pragma once

#include "types.h"
/**
* \brief Stores the grid information
*/
class domain
{
public:
	int   nx, ///< Number of cells in the x-direction
	      ny; ///< Number of cells in the y-direction
	
	vecH  x,  ///< x-coordinates of the nodes
	      y,  ///< y-coordinates of the nodes
	      dx, ///< Cell widths in the x-direction
	      dy; ///< Cell widths in the y-direction
	
	vecD  xD,  ///< x-coordinates of the nodes stored on the device
	      yD,  ///< y-coordinates of the nodes stored on the device
	      dxD, ///< x- cell widths stored on the device
	      dyD; ///< y- cell widths stored on the device
	
	vecH  xu,  ///< x-coordinates of the locations at which the x-component of velocity is evaluated
	      yu,  ///< y-coordinates of the locations at which the x-component of velocity is evaluated
	      xv,  ///< x-coordinates of the locations at which the y-component of velocity is evaluated
	      yv;  ///< y-coordinates of the locations at which the y-component of velocity is evaluated
};
