/***************************************************************************//**
 * \file generateBC1.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the kernels to generate right hand-side terms of the
 *        intermediate velocity flux solver.
 */


#pragma once

#include <types.h>


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

// compute inhomogeneous term from the discrete Laplacian operator
// for the u-velocity at a given boundary with a Dirichlet-type condition
__global__
void bc1DirichletU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc);

// compute inhomogeneous term from the discrete Laplacian operator
// for the v-velocity at a given boundary with a Dirichlet-type condition
__global__
void bc1DirichletV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dy, real C, real *bc, int numUbc);

// compute inhomogeneous term from the discrete Laplacian operator
// for the u-velocity at a given boundary with a convective-type condition
__global__
void bc1ConvectiveU(real *bc1, int N, int nx, int offset, int stride, real *dx, real *dy, real C, real *bc, real *q, real alpha);

// compute inhomogeneous term from the discrete Laplacian operator
// for the v-velocity at a given boundary with a convective-type condition
__global__
void bc1ConvectiveV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dx, real *dy, real C, real *bc, int numUbc, real *q, real alpha);

// compute inhomogeneous term from the discrete Laplacian operator
// for the u-velocity at a given boundary with a special-type condition
__global__
void bc1SpecialU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc, real time);
	
} // end of namespace kernels
