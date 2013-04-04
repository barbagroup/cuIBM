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

#include <solvers/NavierStokes/kernels/generateBC2.h>

template <>
void NavierStokesSolver<host_memory>::generateBC2()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	//int  numP  = nx*ny;
	
	real *dx = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	/// rhs2 of size np+2nb
	
	/// boundary conditions from the continuity equations
		
	cusp::blas::fill(bc2, 0.0);
	
	for(int i=0; i<nx; i++)
	{
		bc2[i] -= bc[YMINUS][i+nx-1]*dx[i]; // v[0][i+1]*dx;
		bc2[(ny-1)*nx + i] += bc[YPLUS][i+nx-1]*dx[i]; // v[ny][i+1]*dx;
	}

	for(int j=0; j<ny; j++)
	{
		// left
		bc2[j*nx] -= bc[XMINUS][j]*dy[j]; // u[j+1][0]*dy;
		// right
		bc2[j*nx+nx-1] += bc[XPLUS][j]*dy[j]; // u[j+1][nx]*dy;
	}
}

template<>
void NavierStokesSolver<device_memory>::generateBC2()
{
	// raw pointers for cup arrays
	real *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	real *bc2_r = thrust::raw_pointer_cast(&(bc2[0])),
	     *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0]));
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	cusp::blas::fill(bc2, 0.0);
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx+ny-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	
	kernels::fillBC2_v <<<dimGrid, dimBlock>>> (bc2_r, yminus, yplus, dxD, nx, ny);
	kernels::fillBC2_u <<<dimGrid, dimBlock>>> (bc2_r, xminus, xplus, dyD, nx, ny);
}
