/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

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
