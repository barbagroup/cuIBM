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
* @file  bodies.h
* @brief Stores information about the body points in arrays
*/

#pragma once

#include <domain.h>
#include <parameterDB.h>
#include <body.h>

template <typename memoryType>
class bodies
{
public:
	int  numBodies,   ///< number of bodies
	     totalPoints; ///< total number of boundary points (all bodies)

	bool bodiesMove;  ///< Tells whether the body is moving or not

	cusp::array1d<int, memoryType>
		numPoints,    ///< number of points in a body
		offsets,      ///< array index of the first point of a body
		I,            ///< x-index of the cell in which a body point is present
		J;            ///< y-index of the cell in which a body point is present
	
	cusp::array1d<int, memoryType>
		startI,       ///< Starting cell index of the bounding box of a body
		startJ,       ///< Starting cell index of the bounding box of a body
		numCellsX,    ///< Number of cells in the x-direction in the bounding box of a body
		numCellsY;    ///< Number of cells in the y-direction in the bounding box of a body
		
	cusp::array1d<real, memoryType>
		xmin,  ///< Lowest x-coordinate for the bounding box of a body
		xmax,  ///< Highest x-coordinate for the bounding box of a body
		ymin,  ///< Lowest y-coordinate for the bounding box of a body
		ymax;  ///< Highest y-coordinate for the bounding box of a body

	cusp::array1d<real, memoryType>
		X,     ///< reference x-coordinates of the boundary points
		Y,     ///< reference y-coordinates of the boundary points
		ds,    ///< vector containing the lengths of all the boundary segments
		ones,  ///< vector of size \link totalPoints \endlink with all elements 1
		x,     ///< actual x-coordinate of the boundary points
		y,     ///< actual y-coordinate of the boundary points
		uB,    ///< x-velocity of the boundary points
		vB;    ///< y-velocity of the boundary points

	/**
	* @brief Initialise the arrays in the class with information from the
	*        @link body @endlink instances.
	*
	* Information regarding the coordinates of the body points and the motion 
	* of the bodies is stored on the host as an array of instances of the 
	* class @link body @endlink. This function transfers that information to
	* arrays on the device, where they are stored as a structure of arrays.  
	* This makea computation more efficient.
	*
	* @param db Database that contains all the simulation parameters
	* @param D  Information about the computaional grid.
	*/
	void initialise(parameterDB &db, domain &D);
	
	/**
	* @brief Calculates the indices of the cells in which the boundary points 
	*        are present.
	* 
	* This information is useful when transferring data between the 
	* boundary points and the computational grid.
	*
	* @param D Information about the computational grid.
	*/
	void calculateCellIndices(domain &D);
	
	/**
	* @brief Calculates the bounding boxes for each body.
	* @param D Information about the computational grid.
	*/
	void calculateBoundingBoxes(domain &D);
	
	/**
	* @brief Update the locations of the body points
	* This is done using the formulae:
	*
	* \f$x_{i,m} = X^c_m + (X_{i,m} - X^0_m) \cos\theta - (Y_{i,m} - Y^0_m) \sin\theta\f$ and
	*
	* \f$y_{i,m} = Y^c_m + (X_{i,m} - X^0_m) \sin\theta + (Y_{i,m} - Y^0_m) \cos\theta\f$
	*/
	void update(parameterDB &db, domain &D, real Time);
};
