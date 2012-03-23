#pragma once

#include <flowDescription.h>
#include <domain.h>

class bodies
{
public:
	int  numBodies;
	int  *numPoints, totalPoints;
	int  *offsets;
	real *X, *Y,	///< Reference location of boundary points
	     *ds;
	real *x, *y;
	int  *I, *J;
	real *uB, *vB;
	bool bodiesMove;
	
	real *X0_x, *X0_y,	///< Reference centre of rotation
	     *Xc_x, *Xc_y,	///< Actual centre of rotation
	     *Theta0;		///< Reference angle
	
	void initialise(flowDescription &F, domain &D);
	void calculateCellIndices(domain &D);
	void update();
};