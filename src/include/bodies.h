/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
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

#pragma once

#include <domain.h>
#include <parameterDB.h>
#include <body.h>

template <typename memoryType>
class bodies
{
public:
	int  numBodies,
	     totalPoints;

	bool bodiesMove;

	cusp::array1d<int, memoryType>
		numPoints, offsets,
		I, J;
	
	cusp::array1d<int, memoryType>
		startI, startJ, numCellsX, numCellsY;
		
	cusp::array1d<real, memoryType>
		xmin, xmax, ymin, ymax;

	cusp::array1d<real, memoryType>
		X, Y, ds,
		x, y, uB, vB,
		X0_x, X0_y, Xc_x, Xc_y, Theta0;

	void initialise(parameterDB &db, domain &D);
	void calculateCellIndices(domain &D);
	void calculateBoundingBoxes(domain &D);
	
	void update(parameterDB &db);
};
