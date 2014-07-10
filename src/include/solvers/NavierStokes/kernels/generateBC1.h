/***************************************************************************//**
* \file generateBC1.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate vector bc1
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
void bc1DirichletU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc);

/**
* \brief To be documented
*/
__global__
void bc1DirichletV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dy, real C, real *bc, int numUbc);

/**
* \brief To be documented
*/
__global__
void bc1ConvectiveU(real *bc1, int N, int nx, int offset, int stride, real *dx, real *dy, real C, real *bc, real *q, real alpha);

/**
* \brief To be documented
*/
__global__
void bc1ConvectiveV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dx, real *dy, real C, real *bc, int numUbc, real *q, real alpha);

/**
* \brief To be documented
*/
__global__
void bc1SpecialU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc, real time);
	
} // end of namespace kernels
