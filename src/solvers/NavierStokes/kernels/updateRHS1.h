/***************************************************************************//**
* \file updateRHS1.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the CUDA kernels required to update vector \c rhs1
*/

#pragma once

#include <types.h>

/**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
__global__
void updateRHS1(real *rhs1, int numUV, int *tags);

/**
* \brief To be documented
*/
__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *cf, real *uv);

/**
* \brief To be documented
*/
__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *cf, real *uv);

} // end of namespace kernels
