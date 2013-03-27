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

#include <solvers/NavierStokes/kernels/generateA.h>

namespace kernels
{

__global__
void generateA(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha)
{
	for (int I=threadIdx.x + blockIdx.x*blockDim.x; I<ASize; I += blockDim.x*gridDim.x)
	{
		ARows[I] = LRows[I];
		ACols[I] = LCols[I];
		AVals[I] = -alpha*LVals[I] + (LRows[I]==LCols[I])*MVals[LRows[I]];
	}
}

__global__
void generateAFadlun(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha, int *tagsX, int *tagsY)
{
	for(int I=threadIdx.x + blockIdx.x*blockDim.x; I<ASize; I += blockDim.x*gridDim.x)
	{
		ARows[I] = LRows[I];
		ACols[I] = LCols[I];
		AVals[I] = - (tagsX[LRows[I]] == -1 && tagsY[LRows[I]] == -1)*alpha*LVals[I] // if the current location is untagged, add -alpha*L
		           - (tagsX[LRows[I]] != -1 || tagsY[LRows[I]] != -1)*LVals[I] // if the current location is tagged, add -L
		           + (LRows[I]==LCols[I])*MVals[LRows[I]]; // if it is a diagonal, add M
	}
}
	
} // end of namespace kernels
