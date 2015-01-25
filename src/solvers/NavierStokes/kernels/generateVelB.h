/***************************************************************************//**
 * \file generateVelB.h
 * \author Anush Krishnan (anush@bu.edu
 * \brief Declaration of the kernels required to generate body-velocities.
 */


#pragma once

#include <types.h>


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

// store an element of the u- and v- body-velocities into one single array
__global__
void fill_velB(real *velB, real *uB, real *vB, int totalPoints);

} // end of namespace kernels
