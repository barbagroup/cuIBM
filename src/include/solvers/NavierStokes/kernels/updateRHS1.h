#pragma once

#include <types.h>

namespace kernels
{
__global__
void updateRHS1(real *rhs1, int numUV, int *tags);

__global__
void updateRHS1(real *rhs1, int numUV, int *tagsX, int *tagsY);

__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tagsX, int *tagsY, real *cfX, real *cfY, real *uvX, real *uvY);

__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tagsX, int *tagsY, real *cfX, real *cfY, real *uvX, real *uvY);

} // end of namespace kernels
