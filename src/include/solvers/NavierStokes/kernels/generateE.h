/***************************************************************************//**
* \file generateE.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate matrix E
*/

#pragma once

#include <types.h>
#include <helpers.h>

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
void generateEHost(int *ERows,  int *ECols,  real *EVals,
                   int nx, int ny, real *x, real *y, real *dx,
                   int totalPoints, real *xB, real *yB, int *I, int *J);

/**
* \brief To be documented
*/
__global__ \
void generateE(int *ERows,  int *ECols,  real *EVals,
               int nx, int ny, real *x, real *y, real *dx,
               int totalPoints, real *xB, real *yB, int *I, int *J);

} // end of namespace kernels
