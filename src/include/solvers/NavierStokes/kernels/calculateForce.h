/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

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
void forceX(real *f, real *q, real *rn, int *tagsX, int *tagsY,
            int nx, int ny, real *dx, real *dy, 
            real dt, real alpha, real nu);

__global__
void forceY();
                   
} // end of namespace kernels
