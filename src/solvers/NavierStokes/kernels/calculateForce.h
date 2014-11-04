/***************************************************************************//**
* \file  calculateForce.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Calculate forces acting on the body using a control volume approach
*/

#pragma once

#include <types.h>

/**
* \namespace kernels
* \brief Contain all the custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief Calculate drag using a control-volume approach (left-right)
*/
__global__
void dragLeftRight(real *FxX, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

/**
* \brief Calculate drag using a control-volume approach (bottom-top)
*/
__global__
void dragBottomTop(real *FxY, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

/**
* \brief Calculate drag using a control-volume approach (unsteady)
*/                   
__global__
void dragUnsteady(real *FxU, real *q, real *qOld, real *dx, real *dy, real dt,
                   int nx, int ny, int I, int J, int ncx, int ncy);

/**
* \brief Calculate lift using a control-volume approach (left-right)
*/
__global__
void liftLeftRight(real *FyX, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

/**
* \brief Calculate lift using a control-volume approach (bottom-right)
*/
__global__
void liftBottomTop(real *FyY, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

/**
* \brief Calculate lift using a control-volume approach (unsteady)
*/
__global__
void liftUnsteady(real *FyU, real *q, real *qOld, real *dx, real *dy, real dt,
                  int nx, int ny, int I, int J, int ncx, int ncy);

/**
* \brief To be documented
*/
__global__
void forceX(real *f, real *q, real *rn, int *tags,
            int nx, int ny, real *dx, real *dy, 
            real dt, real alpha, real nu);

/**
* \brief Doing nothing !
*/
__global__
void forceY();
                   
} // end of namespace kernels
