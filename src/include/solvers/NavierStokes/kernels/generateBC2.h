/***************************************************************************//**
* \file generateBC2.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate elements of vector \c bc2
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
* \brief Kernel to generate Dirichlet boundary condition term
*		coming from the discrete divergence operator
*		at bottom and top v-velocity boundary points
*/
__global__
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny);

/**
* \brief Kernel to generate Dirichlet boundary condition term
*		coming from the discrete divergence operator
*		at left and right u-velocity boundary points
*/
__global__
void fillBC2_u(real *bc2, real *xminus, real *xplus, real *dy, int nx, int ny);

/**
* \brief Kernel to generate no-slip condition terms (x- and y- components)
*		at a body-point location
*/
__global__
void fillBC2_uvB(real *bc2, real *uB, real *vB, int totalPoints, int nx, int ny);

} // end of namespace kernels
