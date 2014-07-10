/***************************************************************************//**
* \file generateM.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate matrix M
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
* \brief Assemble the part of M that corresponds to the x-velocity values
*/
__global__
void fillM_u(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

/**
* \brief Assemble the part of M that corresponds to the y-velocity values
*/
__global__
void fillM_v(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

} // end of namespace kernels
