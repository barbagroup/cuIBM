#pragma once

#include <flowDescription.h>
#include <domain.h>

template <typename memoryType>
class bodies
{
public:
	int  numBodies,
	     totalPoints;
	
	bool bodiesMove;

	cusp::array1d<real, memoryType>
		X, Y, ds,
		x, y, uB, vB,
		X0_x, X0_y, Xc_x, Xc_y, Theta0;

	cusp::array1d<int, memoryType>
		numPoints, offsets,
		I, J;
	
	void initialise(flowDescription &F, domain &D);
	void calculateCellIndices(domain &D);
	void update();
};