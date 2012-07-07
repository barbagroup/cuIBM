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

#include <solvers/NavierStokes/kernels/generateBC2.h>
namespace kernels
{

__global__
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=nx)
		return;
	bc2[i] -= yminus[i+nx-1]*dx[i];
	bc2[(ny-1)*nx + i] += yplus[i+nx-1]*dx[i];
}

__global__
void fillBC2_u(real *bc2, real *xminus, real *xplus, real *dy, int nx, int ny)
{
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	if(j>=ny)
		return;
	bc2[j*nx] -= xminus[j]*dy[j];
	bc2[j*nx+nx-1] += xplus[j]*dy[j];
}

__global__
void fillBC2_uvB(real *bc2, real *uB, real *vB, int totalPoints, int nx, int ny)
{
	int k = threadIdx.x + blockIdx.x*blockDim.x;
	if(k>=totalPoints)
		return;
	bc2[nx*ny + k] = uB[k];
	bc2[nx*ny + k + totalPoints] = vB[k];
}

}