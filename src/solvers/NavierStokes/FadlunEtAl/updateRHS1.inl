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

#include <solvers/NavierStokes/kernels/updateRHS1.h>

template <>
void FadlunEtAlSolver<device_memory>::updateRHS1()
{
	const int blocksize = 256;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numUV = (nx-1)*ny + nx*(ny-1);
	
	real *rhs1_r   = thrust::raw_pointer_cast(&(rhs1[0]));
//	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));
	int  *tagsX_r = thrust::raw_pointer_cast(&(tagsXD[0])),\
	     *tagsY_r = thrust::raw_pointer_cast(&(tagsYD[0]));
	
	dim3 dimGrid( int((numUV-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
	// 1-d interpolation
//	kernels::updateRHS1 <<<dimGrid, dimBlock>>> (rhs1_r, numUV, tags_r);

	// 2-d interpolation
	kernels::updateRHS1 <<<dimGrid, dimBlock>>> (rhs1_r, numUV, tagsX_r, tagsY_r);
}

template <>
void FadlunEtAlSolver<host_memory>::updateRHS1()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny,
		 numUV = (nx-1)*ny + nx*(ny-1);
	
	for(in j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			//if(tags[I]!=-1) // this is for 1-d interpolation
		
			// MODIFY THIS FOR MOVING BODIES
			
			// 2-d interpolation
			if(tagsX[I]!=-1 || tagsY[I]!=-1)
				rhs1[I] = 0.0;
		
			// MULTIPLY WITH SCALE AND DT AND DIVIDE BY FACE WIDTH
			//rhs1[I] = (tagsX[I]==-1 && tagsY[I]==-1)*rhs1[I] + (tagsX[I]!=-1)*(1.0-coeffsX[I])*(1.0-coeffsY[I])*uX[I] + (tagsY[I]!=-1)*coeffsX[I]*(1.0-coeffsY[I])*uY[I];
		}
	}
	
	for(in j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = numU + j*nx + i;
			//if(tags[i]!=-1) // this is for 1-d interpolation
		
			// MODIFY THIS FOR MOVING BODIES
			// 2-d interpolation
			if(tagsX[i]!=-1 || tagsY[i]!=-1)
				rhs1[i] = 0.0;
		
			// MULTIPLY WITH SCALE AND DT AND DIVIDE BY FACE WIDTH
			//rhs1[i] = (tagsX[i]==-1 && tagsY[i]==-1)*rhs[i] + (tagsX[i]!=-1)*(1.0-coeffsX[i])*(1.0-coeffsY[i])*vX[i] + (tagsY[i]!=-1)*coeffsX[i]*(1.0-coeffsY[i])*vY[i];
		}
	}
}
