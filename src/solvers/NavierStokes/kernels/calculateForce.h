/***************************************************************************//**
 * \file  calculateForce.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the kernels to calculate the forces acting on a body
 *        The method is described in Lai & Peskin (2000).
 */


#pragma once

#include <types.h>


/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */
namespace kernels
{

__global__
void dragLeftRight(real *FxX, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void dragBottomTop(real *FxY, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void dragUnsteady(real *FxU, real *q, real *qOld, real *dx, real *dy, real dt,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void liftLeftRight(real *FyX, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void liftBottomTop(real *FyY, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void liftUnsteady(real *FyU, real *q, real *qOld, real *dx, real *dy, real dt,
                  int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void forceX(real *f, real *q, real *rn, int *tags,
            int nx, int ny, real *dx, real *dy, 
            real dt, real alpha, real nu);

// doing nothing
__global__
void forceY();
                   
} // end of namespace kernels
