/***************************************************************************//**
* \file generateVelB.h
* \author Krishnan, A. (anush@bu.edu
* \brief Declaration of the CUDA kernels required to generate vector \c velB
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
* \brief CUDA kernel to store element of x- and y- body velocity arrays
*		into one single array
*/
__global__
void fill_velB(real *velB, real *uB, real *vB, int totalPoints);

} // end of namespace kernels
