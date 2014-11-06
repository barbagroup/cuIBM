/***************************************************************************//**
 * \file generateM.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the kernels to generate the mass matrix and its inverse.
 */


#pragma once

#include <types.h>


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

// compute element of the mass matrix and its inverse for a x-velocity node
__global__
void fillM_u(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

// compute element of the mass matrix and its inverse for a y-velocity node
__global__
void fillM_v(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

} // end of namespace kernels
