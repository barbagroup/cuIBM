/***************************************************************************//**
* \file generateBC2.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate vector bc2
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
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny);

/**
* \brief To be documented
*/
__global__
void fillBC2_u(real *bc2, real *xminus, real *xplus, real *dy, int nx, int ny);

/**
* \brief To be documented
*/
__global__
void fillBC2_uvB(real *bc2, real *uB, real *vB, int totalPoints, int nx, int ny);

} // end of namespace kernels
