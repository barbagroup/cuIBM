/***************************************************************************//**
 * \file generateRN.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the kernels to generate the explicit terms of the momentum equation.
 */


#pragma once

#include <types.h>


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

// compute u-component of the explicit terms of the discretized momentum equation
// and element of the explicit convective terms
__global__
void convectionTermU(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

// compute v-component of the explicit terms of the discretized momentum equation
// and element of the explicit convective terms
__global__
void convectionTermV(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

// compute v-component of the explicit terms of the discretized momentum equation
// and element of the explicit convective terms
// top or bottom boundaries
__global__
void convectionTermVBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop);

// compute u-component of the explicit terms of the discretized momentum equation
// and element of the explicit convective terms
// top or bottom boundaries
__global__
void convectionTermUBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

// compute u-component of the explicit terms of the discretized momentum equation
// and element of the explicit convective terms
// left or right boundaries
__global__
void convectionTermULeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcLeft, real *bcRight);

// compute v-component of the explicit terms of the discretized momentum equation
// and element of the explicit convective terms
// left or right boundaries
__global__
void convectionTermVLeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

} // end of namespace kernels
