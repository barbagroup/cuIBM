/***************************************************************************//**
 * \file generateBC2.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the kernels to generate elements of the right hand-side
 *        of the Poisson solver.
 */


#pragma once

#include <types.h>


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

// compute inhomogeneous terms of the discrete divergence operator
// from the bottom and top boundaries at the v-velocity locations
__global__
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny);

// compute inhomogeneous terms of the discrete divergence operator
// from the left and right boundaries at the u-velocity locations
__global__
void fillBC2_u(real *bc2, real *xminus, real *xplus, real *dy, int nx, int ny);

// compute inhomogeneous terms of the discrete divergence operator
// from the no-slip constraint at the body-point locations
__global__
void fillBC2_uvB(real *bc2, real *uB, real *vB, int totalPoints, int nx, int ny);

} // end of namespace kernels
