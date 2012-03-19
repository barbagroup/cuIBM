#pragma once

class bodies
{
public:
	int  numBodies;
	int  *offsets;
	real *x, *y;
	bool bodiesMove;
	
	real X0_x, X0_y,	///< Reference centre of rotation
	     Xc_x, Xc_y,	///< Actual centre of rotation
	     Theta0;		///< Reference angle

	real *X, *Y,	///< Reference location of boundary points
	     *ds;	///< segment lengths
	
	void update();
};