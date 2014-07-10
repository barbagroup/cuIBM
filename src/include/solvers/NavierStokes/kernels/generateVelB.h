/***************************************************************************//**
* \file generateVelB.h
* \author Krishnan, A. (anush@bu.edu
* \brief Declaration of the kernels required to generate vector VelB
*/

#pragma once

#include <types.h>

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
__global__
void fill_velB(real *velB, real *uB, real *vB, int totalPoints);

} // end of namespace kernels
