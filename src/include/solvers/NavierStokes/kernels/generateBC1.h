/***************************************************************************//**
* \file generateBC1.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the kernels required to generate elements of vector \c bc1
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
*		coming from the discrete Laplacian operator
*		applied to the u-velocity at a given boundary
*/
__global__
void bc1DirichletU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc);

/**
* \brief Kernel to generate Dirichlet boundary condition term
*		coming from the discrete Laplacian operator
*		applied to the v-velocity at a given boundary
*/
__global__
void bc1DirichletV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dy, real C, real *bc, int numUbc);

/**
* \brief Kernel to generate convective boundary condition term
*		coming from the discrete Laplacian operator
*		applied to the u-velocity at a given boundary
*/
__global__
void bc1ConvectiveU(real *bc1, int N, int nx, int offset, int stride, real *dx, real *dy, real C, real *bc, real *q, real alpha);

/**
* \brief Kernel to generate convective boundary condition term
*		coming from the discrete Laplacian operator
*		applied to the v-velocity at a given boundary
*/
__global__
void bc1ConvectiveV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dx, real *dy, real C, real *bc, int numUbc, real *q, real alpha);

/**
* \brief Kernel to generate "special" boundary condition term
*		coming from the discrete Laplacian operator
*		applied to the u-velocity at a given boundary
*/
__global__
void bc1SpecialU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc, real time);
	
} // end of namespace kernels
