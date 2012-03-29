#pragma once

#include <types.h>
namespace kernels
{

__global__
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny);

__global__
void fillBC2_u(real *bc2, real *xminus, real *xplus, real *dy, int nx, int ny);

__global__
void fillBC2_uvB(real *bc2, real *uB, real *vB, int totalPoints, int nx, int ny);

}