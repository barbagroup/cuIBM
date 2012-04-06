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

	cusp::array1d<real, memoryType>
		X, Y, ds,
		x, y, uB, vB,
		X0_x, X0_y, Xc_x, Xc_y, Theta0;

	void initialise(parameterDB &db, domain &D);
	void calculateCellIndices(domain &D);
	void update();
};
