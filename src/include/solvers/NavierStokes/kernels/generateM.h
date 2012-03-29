#pragma once

#include <types.h>

namespace kernels
{

__global__
void fillM_u(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

__global__
void fillM_v(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt);

} // end of namespace kernel