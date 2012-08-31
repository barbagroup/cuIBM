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

template <>
void FadlunEtAlSolver<device_memory>::generateA(int alpha)
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	int  ASize = 5*( (nx-1)*ny + nx*(ny-1) ) - 4*(nx+ny) + 4,
	     numUV = (nx-1)*ny + nx*(ny-1);

	A.resize(numUV, numUV, ASize);

	int  *LRows = thrust::raw_pointer_cast(&(L.row_indices[0])),
	     *LCols = thrust::raw_pointer_cast(&(L.column_indices[0])),
	     *ARows = thrust::raw_pointer_cast(&(A.row_indices[0])),
	     *ACols = thrust::raw_pointer_cast(&(A.column_indices[0])),
	     *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));

	real *MVals = thrust::raw_pointer_cast(&(M.values[0])),
	     *LVals = thrust::raw_pointer_cast(&(L.values[0])),
	     *AVals = thrust::raw_pointer_cast(&(A.values[0]));

	const int blockSize = 256;
	dim3 dimGrid( int((ASize-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);

	kernels::generateAFadlun <<<dimGrid, dimBlock>>> (ARows, ACols, AVals, MVals, LRows, LCols, LVals, ASize, alpha, tags_r);
}

template <>
void FadlunEtAlSolver<host_memory>::generateA(int alpha)
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	int  ASize = 5*( (nx-1)*ny + nx*(ny-1) ) - 4*(nx+ny) + 4,
	     numUV = (nx-1)*ny + nx*(ny-1);

	A.resize(numUV, numUV, ASize);

	for(int i=0; i<ASize; i++)
	{
		A.row_indices[i] = L.row_indices[i];
		A.column_indices[i] = L.column_indices[i];
		if(tags[A.row_indices[i]]==-1)
			A.values[i] = -alpha*L.values[i];
		else
			A.values[i] = -L.values[i];
		if(A.row_indices[i] == A.column_indices[i])
			A.values[i] += M.values[A.row_indices[i]];
	}
}
