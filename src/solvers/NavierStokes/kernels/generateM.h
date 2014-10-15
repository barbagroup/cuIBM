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
