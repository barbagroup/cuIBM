/**
* @file  generateVelB.h
* @brief Contains kernals required for generation of vector VelB.
*/

#pragma once

#include <types.h>

namespace kernels
{

__global__
void fill_velB(real *velB, real *uB, real *vB, int totalPoints);

}
