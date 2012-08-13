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

#pragma once

#include <types.h>

namespace kernels
{

/**
* @brief Assemble the part of M that correspond to the x-velocity values
* @param MRows    Raw pointer to the array storing the row indices of the elements of M
* @param MCols    Raw pointer to the array storing the column indices of the elements of M
* @param MVals    Raw pointer to the array storing the values of the elements of M
* @param MinvRows Raw pointer to the array storing the row indices of the elements of M-inverse
* @param MinvCols Raw pointer to the array storing the column indices of the elements of M-inverse
* @param MinvVals Raw pointer to the array storing the values of the elements of M-inverse
* @param nx       Number of cells in the domain in the x-direction
* @param ny       Number of cells in the domain in the y-direction
* @param dx       Cell widhts in the x-direction
* @param dy       Cell widths in the y-direction
* @param dt       Time increment
*/
__global__
void fillM_u(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

/**
* @brief Assemble the part of M that correspond to the y-velocity values
* @param MRows    Raw pointer to the array storing the row indices of the elements of M
* @param MCols    Raw pointer to the array storing the column indices of the elements of M
* @param MVals    Raw pointer to the array storing the values of the elements of M
* @param MinvRows Raw pointer to the array storing the row indices of the elements of M-inverse
* @param MinvCols Raw pointer to the array storing the column indices of the elements of M-inverse
* @param MinvVals Raw pointer to the array storing the values of the elements of M-inverse
* @param nx       Number of cells in the domain in the x-direction
* @param ny       Number of cells in the domain in the y-direction
* @param dx       Cell widhts in the x-direction
* @param dy       Cell widths in the y-direction
* @param dt       Time increment
*/
__global__
void fillM_v(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

} // end of namespace kernels
