/***************************************************************************//**
* \file generateRN.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate vector \c rn
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
* \brief CUDA kernel to compute element of the vector \c rn and \c H,
*		that corresponds to a x-velocity value
*/
__global__
void convectionTermU(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

/**
* \brief CUDA kernel to compute element of the vector \c rn and \c H,
*		that corresponds to a y-velocity value
*/
__global__
void convectionTermV(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

/**
* \brief CUDA kernel to compute element of the vector \c rn and \c H,
*		that corresponds to a y-velocity value
*		on the bottom or top boundaries
*/
__global__
void convectionTermVBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop);

/**
* \brief CUDA kernel to compute element of the vector \c rn and \c H,
*		that corresponds to a x-velocity value
* 		on the bottom or top boundaries
*/
__global__
void convectionTermUBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

/**
* \brief CUDA kernel to compute element of the vector \c rn and \c H,
*		that corresponds to a x-velocity value
*		on the left or right boundaries
*/
__global__
void convectionTermULeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcLeft, real *bcRight);

/**
* \brief CUDA kernel to compute element of the vector \c rn and \c H,
*		that corresponds to a y-velocity value
*		on the left or right boundaries
*/
__global__
void convectionTermVLeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

} // end of namespace kernels
