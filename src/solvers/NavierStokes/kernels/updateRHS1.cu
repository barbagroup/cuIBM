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

#define BSZ 16

namespace kernels
{

// 1-d interpolation
/*__global__
void updateRHS1(real *rhs1, int numUV, int *tags)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I>=numUV)
		return;
	
	rhs1[I] = rhs1[I]*(tags[I]==-1);
}*/

// 2-d interpolation
__global__
void updateRHS1(real *rhs1, int numUV, int *tagsX, int *tagsY)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I>=numUV)
		return;
	
	rhs1[I] = rhs1[I]*(tagsX[I]==-1 && tagsY[I]==-1);
}

__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tagsX, int *tagsY, real *coeffsX, real *coeffsY, real *uvX, real *uvY)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	int i = I % (nx-1);
	
	if( I < (nx-1)*ny )
	{
		rhs1[I] = 
		rhs1[I] = (tagsX[I]==-1 && tagsY[I]==-1)*rhs1[I] + ( (tagsX[I]!=-1)*(1.0-coeffsX[I])*(1.0-coeffsY[I])*uvX[I] + (tagsY[I]!=-1)*coeffsX[I]*(1.0-coeffsY[I])*uvY[I] ) * 0.5*(dx[i+1]+dx[i]) / dt;
	}
}

__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tagsX, int *tagsY, real *coeffsX, real *coeffsY, real *uvX, real *uvY)
{
	int numU = (nx-1)*ny;
	int	I = blockIdx.x*blockDim.x + threadIdx.x + numU;
	int j = (I-numU) / nx;
	
	if( I < numU + nx*(ny-1) )
	{
		rhs1[I] = (tagsX[I]==-1 && tagsY[I]==-1)*rhs1[I] + ( (tagsX[I]!=-1)*(1.0-coeffsX[I])*(1.0-coeffsY[I])*uvX[I] + (tagsY[I]!=-1)*coeffsX[I]*(1.0-coeffsY[I])*uvY[I] ) * 0.5*(dy[j+1]+dy[j]) / dt;
	}
}

} // end of namespace kernels
