#pragma once

#include <types.h>

namespace kernels
{

__global__
void convectionTermU(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha);

__global__
void convectionTermV(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha);

__global__
void convectionTermVBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, \
                              real *bcBottom, real *bcTop);

__global__
void convectionTermUBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

__global__
void convectionTermULeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, \
                              real *bcLeft, real *bcRight);

__global__
void convectionTermVLeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

__global__
void updateRN(real *rn, int numUV, int *tags);

__global__
void updateRN(real *rn, int numUV, int *tagsX, int *tagsY);

} // end of namespace kernels
