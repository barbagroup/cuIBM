/**
* @file  generateA.h
* @breif Contains kernels required for the generation of matrix A.
*/
#pragma once

#include <types.h>

namespace kernels
{

__global__
void generateA(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha);

__global__
void generateAFadlun(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha, int *tagsX, int *tagsY);

} // end of namespace kernels
