/***************************************************************************//**
* \file generateRN.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate vector rn
*/

#pragma once

#include <types.h>

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
__global__
void convectionTermU(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

/**
* \brief To be documented
*/
__global__
void convectionTermV(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

/**
* \brief To be documented
*/
__global__
void convectionTermVBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop);

/**
* \brief To be documented
*/
__global__
void convectionTermUBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

/**
* \brief To be documented
*/
__global__
void convectionTermULeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcLeft, real *bcRight);

/**
* \brief To be documented
*/
__global__
void convectionTermVLeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

} // end of namespace kernels
