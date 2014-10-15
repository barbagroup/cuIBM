/**
* @file  calculateForce.h
* @brief Calculates force on the body using a control volume approach
*/

#pragma once

#include <types.h>

/**
* @namespace kernels
* @brief     Contains all the custom-written CUDA kernels.
*/
namespace kernels
{

/**
* @brief Used to calculate drag using the control volume approach.
*
* This function evaluates the contribution of the left and right parts of the 
* control surface.
*
* @param q   Raw pointer to the vector storing all the fluxes
* @param dx  Raw pointer to the vector storing the cell widths in the x-direction
* @param dy  Raw pointer to the vector storing the cell widths in the y-direction
* @param nx  Number of cells in the x-direction
* @param ny  Number of cells in the y-direction
* @param I   x-index of bottom-left corner cell of the control surface
* @param J   y-index of top-right corner cell of the control surface
* @param ncx Number of cells in the x-direction in the control volume
* @param ncy Number of cells in the y-direction in the control volume
*/
__global__
void dragLeftRight(real *FxX, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

/**
* @brief Used to calculate drag using the control volume approach. 
*
* This function evaluates the contribution of the bottom and top parts of the 
* control surface.
*
* @param q   Raw pointer to the vector storing all the fluxes
* @param dx  Raw pointer to the vector storing the cell widths in the x-direction
* @param dy  Raw pointer to the vector storing the cell widths in the y-direction
* @param nx  Number of cells in the x-direction
* @param ny  Number of cells in the y-direction
* @param I   x-index of bottom-left corner cell of the control surface
* @param J   y-index of top-right corner cell of the control surface
* @param ncx Number of cells in the x-direction in the control volume
* @param ncy Number of cells in the y-direction in the control volume
*/
__global__
void dragBottomTop(real *FxY, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

/**
* @brief Used to calculate drag using the control volume approach. 
*
* This function evaluates the unsteady contribution of the control volume.
*
* @param q    Raw pointer to the vector storing all the fluxes
* @param qOld Raw pointer to the flux vector at the previous time step
* @param dx   Raw pointer to the vector storing the cell widths in the x-direction
* @param dy   Raw pointer to the vector storing the cell widths in the y-direction
* @param dt   Time increment
* @param nx   Number of cells in the x-direction
* @param ny   Number of cells in the y-direction
* @param I    x-index of bottom-left corner cell of the control surface
* @param J    y-index of top-right corner cell of the control surface
* @param ncx  Number of cells in the x-direction in the control volume
* @param ncy  Number of cells in the y-direction in the control volume
*/                   
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

__global__
void forceY();
                   
} // end of namespace kernels
