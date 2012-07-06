/**
* \file domain.h
* \brief Stores the grid information
*/
#pragma once

#include <types.h>
/**
* \brief Stores the grid information
*/
class domain
{
public:
	int   nx, ny;
	
	vecH  x, y,   ///< coordinates of the nodes
	      dx, dy; ///< cell widths
	
	vecD  xD, yD,
	      dxD, dyD;
	
	vecH  xu, yu, xv, yv;
};
