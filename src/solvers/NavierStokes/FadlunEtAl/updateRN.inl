/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
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

#include <solvers/NavierStokes/kernels/generateRN.h>

template <>
void FadlunEtAlSolver<device_memory>::updateRN()
{
	const int blocksize = 256;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numUV = (nx-1)*ny + nx*(ny-1);
	
	real *rn_r   = thrust::raw_pointer_cast(&(rn[0]));
//	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));
	int  *tagsX_r = thrust::raw_pointer_cast(&(tagsXD[0])),\
	     *tagsY_r = thrust::raw_pointer_cast(&(tagsYD[0]));
	
	dim3 dimGrid( int((numUV-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
//	kernels::updateRN <<<dimGrid, dimBlock>>> (rn_r, numUV, tags_r);
	kernels::updateRN <<<dimGrid, dimBlock>>> (rn_r, numUV, tagsX_r, tagsY_r);
}

template <>
void FadlunEtAlSolver<host_memory>::updateRN()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numUV = (nx-1)*ny + nx*(ny-1);
	
	for(int i=0; i<numUV; i++)
	{
//		if(tags[i]!=-1)
		if(tagsX[i]!=-1 || tagsY[i]!=-1)
			rn[i] = 0.0;
	}
}
