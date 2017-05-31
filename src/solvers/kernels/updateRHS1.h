/**
 * \file updateRHS1.h
 * \brief Declaration of the kernels to update the right hand-side 
 *        of the intermediate velocity flux solver.
 *
 * It replaces the right-hand side elements at the tagged points, with
 * values obtained from the interpolation relations at those points.
 */


#pragma once

#include "utilities/types.h"


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

// update the RHS of the velocity system with zeros at forcing points
__global__
void updateRHS1(real *rhs1, int numUV, int *tags);

// 1D linear interpolation along gridline in x-direction
// update the RHS of the velocity system for the x-components of the velocity
__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *cf, real *uv);

// 1D linear interpolation along gridline in y-direction
// update the RHS of the velocity system for the y-components of the velocity
__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *cf, real *uv);

// quadratic interpolation along gridline in x-direction
// update the RHS of the velocity system for the x-components of the velocity
__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *cf, real *cf2, real *uv);

// quadratic interpolation along gridline in x-direction
// update the RHS of the velocity system for the y-components of the velocity
__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *cf, real *cf2, real *uv);

} // End of namespace kernels
