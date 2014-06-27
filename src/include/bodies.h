/***************************************************************************//**
* \file  bodies.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the class \c bodies.
*
* Stores information about the body points in arrays.
*
*/

#pragma once

#include <domain.h>
#include <parameterDB.h>
#include <body.h>

/***************************************************************************//**
* \class bodies
* \brief Stores information about the body
*/
template <typename memoryType>
class bodies
{
public:
	int  numBodies,   ///< number of bodies
	     totalPoints; ///< total number of boundary points (all bodies)

	bool bodiesMove;  ///< tells whether the body is moving or not

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
		forceX,
		forceY;

	cusp::array1d<real, memoryType>
		X,     ///< reference x-coordinates of the boundary points
		Y,     ///< reference y-coordinates of the boundary points
		ds,    ///< vector containing the lengths of all the boundary segments
		ones,  ///< vector of size \link totalPoints \endlink with all elements 1
		x,     ///< actual x-coordinate of the boundary points
		y,     ///< actual y-coordinate of the boundary points
		uB,    ///< x-velocity of the boundary points
		vB;    ///< y-velocity of the boundary points

	/********************//**
	* \brief Initializes the arrays in the class with information from \c body instances
	*/
	void initialise(parameterDB &db, domain &D);
	
	/********************//**
	* \brief Calculates the indices of the cells in which the boundary points are present
	*/
	void calculateCellIndices(domain &D);
	
	/********************//**
	* \brief Calculates the bounding boxes for each body
	*/
	void calculateBoundingBoxes(parameterDB &db, domain &D);
	
	/********************//**
	* \brief Updates the locations of the body points
	*/
	void update(parameterDB &db, domain &D, real Time);

	void writeToFile(std::string &caseFolder, int timeStep);

	void writeToFile(real *bx, real *by, std::string &caseFolder, int timeStep);
};
