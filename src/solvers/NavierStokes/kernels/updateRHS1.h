#pragma once

#include <types.h>

namespace kernels
{
__global__
void updateRHS1(real *rhs1, int numUV, int *tags);

__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *cf, real *uv);

__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *cf, real *uv);

__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *cf, real *cf2, real *uv);

__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *cf, real *cf2, real *uv);

} // end of namespace kernels
