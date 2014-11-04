/***************************************************************************//**
* \file generateM.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate the mass matrix \c M and its inverse
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
* \brief Kernel to compute element of the mass matrix \c M, and its inverse,
*		that corresponds to a x-velocity value
*		on the device
*/
__global__
void fillM_u(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

/**
* \brief Kernel to compute element of the mass matrix \c M, and its inverse,
*		that corresponds to a y-velocity value
*		on the device
*/
__global__
void fillM_v(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

} // end of namespace kernels
