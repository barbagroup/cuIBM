#pragma once

#include <types.h>

namespace kernels
{

__global__
void updateQFadlun(int *QRows, int *QCols, real *QVals, int QSize, int *tags);

}
