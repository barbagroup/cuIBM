/***************************************************************************//**
 * \file helpers.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the discrete delta function.
 */


#include "helpers.h"


/**
 * \brief Discrete delta function from Roma et al. (1999).
 *
 * A.M. Roma, C.S. Peskin, M.J. Berger.
 * An adaptative version of the immersed boundary method.
 * J. Comput. Phys. 153, 509-534 (1999).
 *
 * \param x coordinate of the grid point
 * \param h grid-spacing in the uniform region
 *
 * \return the value of the discrete delta function
 */
real dhRoma(real x, real h)
{
	real r = fabs(x)/h;
	
	if(r>1.5)
		return 0.0;
	else if(r>0.5 && r<=1.5)
		return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
	else
		return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
}

/**
 * \brief Two-dimensional discrete delta function.
 *
 * \param x x-coordinate of the grid point
 * \param y y-coordinate of the grid point
 * \param h grid-spacing in the uniform region
 *
 * \return the value of the discrete Delta function in 2D
 */
real delta(real x, real y, real h)
{
	return dhRoma(x, h) * dhRoma(y, h);
}
