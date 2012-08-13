/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

#pragma once

#include <types.h>

namespace kernels
{

__global__
void bc1DirichletU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc);

__global__
void bc1DirichletV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dy, real C, real *bc, int numUbc);
	
__global__
void bc1ConvectiveU(real *bc1, int N, int nx, int offset, int stride, real *dx, real *dy, real C, real *bc, real *q, real alpha);

__global__
void bc1ConvectiveV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dx, real *dy, real C, real *bc, int numUbc, real *q, real alpha);
	
} // end of namespace kernels
