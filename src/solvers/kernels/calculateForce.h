/**
 * \file calculateForce.h
 * \brief Declaration of the kernels to calculate the forces acting on a body
 *        The method is described in Lai & Peskin (2000).
 */


#pragma once

#include "utilities/types.h"


/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */
namespace kernels
{

// calculate drag using a control-volume approach (left-right).
__global__
void dragLeftRight(real *FxX, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

// calculate drag using a control-volume approach (bottom-top).
__global__
void dragBottomTop(real *FxY, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

// calculate drag using a control-volume approach (unsteady).
__global__
void dragUnsteady(real *FxU, real *q, real *qOld, real *dx, real *dy, real dt,
                   int nx, int ny, int I, int J, int ncx, int ncy);

// calculate lift using a control-volume approach (left-right).
__global__
void liftLeftRight(real *FyX, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

// calculate lift using a control-volume approach (bottom-top).
__global__
void liftBottomTop(real *FyY, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

// calculate lift using a control-volume approach (unsteady).
__global__
void liftUnsteady(real *FyU, real *q, real *qOld, real *dx, real *dy, real dt,
                  int nx, int ny, int I, int J, int ncx, int ncy);

// kernel not usable
__global__
void forceX(real *f, real *q, real *rn, int *tags,
            int nx, int ny, real *dx, real *dy, 
            real dt, real alpha, real nu);

// kernel not usable
__global__
void forceY();
                   
} // End of namespace kernels
