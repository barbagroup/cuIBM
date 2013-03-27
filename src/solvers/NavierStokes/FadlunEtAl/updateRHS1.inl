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
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	real *rhs1_r   = thrust::raw_pointer_cast(&(rhs1[0]));
//	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));

	int  *tagsX_r = thrust::raw_pointer_cast(&(tagsXD[0])),
	     *tagsY_r = thrust::raw_pointer_cast(&(tagsYD[0]));

	real *dx_r = thrust::raw_pointer_cast( &(domInfo->dxD[0]) ),
	     *dy_r = thrust::raw_pointer_cast( &(domInfo->dyD[0]) ),
	     *coeffsX_r = thrust::raw_pointer_cast( &(coeffsXD[0]) ),
	     *coeffsY_r = thrust::raw_pointer_cast( &(coeffsYD[0]) ),
	     *uvX_r   = thrust::raw_pointer_cast( &(uvXD[0]) ),
	     *uvY_r   = thrust::raw_pointer_cast( &(uvYD[0]) );

	const int blocksize = 256;
	dim3 dimGrid( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
	// 1-d interpolation
//	kernels::updateRHS1 <<<dimGrid, dimBlock>>> (rhs1_r, numUV, tags_r);

	// 2-d interpolation
	//kernels::updateRHS1 <<<dimGrid, dimBlock>>> (rhs1_r, (nx-1)*ny+nx*(ny-1), tagsX_r, tagsY_r);
	
	kernels::updateRHS1X <<<dimGrid, dimBlock>>>(rhs1_r, nx, ny, dt, dx_r, tagsX_r, tagsY_r, coeffsX_r, coeffsY_r, uvX_r, uvY_r);
	kernels::updateRHS1Y <<<dimGrid, dimBlock>>>(rhs1_r, nx, ny, dt, dy_r, tagsX_r, tagsY_r, coeffsX_r, coeffsY_r, uvX_r, uvY_r);
}

template <>
void FadlunEtAlSolver<host_memory>::updateRHS1()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny,
	     I = 0;
	
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			//if(tags[I]!=-1) // this is for 1-d interpolation
		
			// MODIFY THIS FOR MOVING BODIES
			
			// 2-d interpolation
			//if(tagsX[I]!=-1 || tagsY[I]!=-1)
			//	rhs1[I] = 0.0;

			rhs1[I] = (tagsX[I]==-1 && tagsY[I]==-1)*rhs1[I] + ( (tagsX[I]!=-1)*(1.0-coeffsX[I])*(1.0-coeffsY[I])*uvX[I] + (tagsY[I]!=-1)*coeffsX[I]*(1.0-coeffsY[I])*uvY[I] ) * 0.5*(domInfo->dx[i+1]+domInfo->dx[i]) / dt;
		}
	}
	
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = numU + j*nx + i;
			//if(tags[i]!=-1) // this is for 1-d interpolation
		
			// MODIFY THIS FOR MOVING BODIES
			// 2-d interpolation
			//if(tagsX[i]!=-1 || tagsY[i]!=-1)
			//	rhs1[i] = 0.0;
		
			rhs1[I] = (tagsX[I]==-1 && tagsY[I]==-1)*rhs1[I] + ( (tagsX[I]!=-1)*(1.0-coeffsX[I])*(1.0-coeffsY[I])*uvX[I] + (tagsY[I]!=-1)*coeffsX[I]*(1.0-coeffsY[I])*uvY[I] ) * 0.5*(domInfo->dy[j+1]+domInfo->dy[j]) / dt;
		}
	}
}
