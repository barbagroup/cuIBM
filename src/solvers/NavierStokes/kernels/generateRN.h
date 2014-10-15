#pragma once

#include <types.h>

namespace kernels
{

__global__
void convectionTermU(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

__global__
void convectionTermV(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real gamma, real zeta, real alpha, real nu);

__global__
void convectionTermVBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop);

__global__
void convectionTermUBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

__global__
void convectionTermULeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcLeft, real *bcRight);

__global__
void convectionTermVLeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight);

} // end of namespace kernels
