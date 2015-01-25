/***************************************************************************//**
 * \file updateRHS1.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the kernels to update the right hand-side 
 *        of the intermediate velocity flux solver.
 *        It replaces the right-hand side elements at the tagged points, with
 *        values obtained from the interpolation relations at those points.
 */


#pragma once

#include <types.h>


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

__global__
void updateRHS1(real *rhs1, int numUV, int *tags);

__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *cf, real *uv);

__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *cf, real *uv);

__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *cf, real *cf2, real *uv);

__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *cf, real *cf2, real *uv);

} // end of namespace kernels
