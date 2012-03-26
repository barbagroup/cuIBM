/**
* \file body.h
* \brief
*/
#pragma once
/**
* Class description
*/
class body
{
	public:
		int  numPoints; ///< number of boundary points

		real X0_x, X0_y, ///< Reference centre of rotation
		     Theta0; ///< Reference angle

		real *X, *Y; ///< Reference location of boundary points
		
		bool tFlag, rFlag;

		real VelX, VelY, Omg;
		real AmpX, OmegaX, PhiX,
		     AmpY, OmegaY, PhiY,
		     AmpT, OmegaT, PhiT;
};