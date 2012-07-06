#pragma once

#include <types.h>

namespace kernels
{

__global__
void updateQFadlun(int *QRows, int *QCols, real *QVals, int QSize, int *tags);

__global__
void updateQFadlun(int *QRows, int *QCols, real *QVals, int QSize, int *tagsX, int *tagsY);

}
