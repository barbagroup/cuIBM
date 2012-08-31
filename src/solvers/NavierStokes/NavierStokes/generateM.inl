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

#include <solvers/NavierStokes/kernels/generateM.h>

template <>
void NavierStokesSolver<host_memory>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I;
	real value;
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	M.resize(numUV, numUV, numUV);
	Minv.resize(numUV, numUV, numUV);
	
	for (int j=0; j < ny; j++)
	{
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;
			value = 0.5*(domInfo->dx[i+1]+domInfo->dx[i]) / domInfo->dy[j] / dt;
			
			M.row_indices[I] = I;
			M.column_indices[I] = I;
			M.values[I] = value;
			
			Minv.row_indices[I] = I;
			Minv.column_indices[I] = I;
			Minv.values[I] = 1.0/value;
		}
	}
	for (int j=0; j < ny-1; j++)
	{
		for (int i=0; i < nx; i++)
		{
			I = j*nx + i + numU;
			
			value  = 0.5*(domInfo->dy[j+1]+domInfo->dy[j]) / domInfo->dx[i] / dt;
			
			M.row_indices[I] = I;
			M.column_indices[I] = I;
			M.values[I] = value;
			
			Minv.row_indices[I] = I;
			Minv.column_indices[I] = I;
			Minv.values[I] = 1.0/value;
		}
	}
}

template<>
void NavierStokesSolver<device_memory>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	real *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	M.resize(numUV, numUV, numUV);
	Minv.resize(numUV, numUV, numUV);
	
	int  *MRows = thrust::raw_pointer_cast(&(M.row_indices[0])),
	     *MCols = thrust::raw_pointer_cast(&(M.column_indices[0])),
	     *MinvRows = thrust::raw_pointer_cast(&(Minv.row_indices[0])),
	     *MinvCols = thrust::raw_pointer_cast(&(Minv.column_indices[0]));
	
	real *MVals = thrust::raw_pointer_cast(&(M.values[0])),
	     *MinvVals = thrust::raw_pointer_cast(&(Minv.values[0]));
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx*ny-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	kernels::fillM_u <<<dimGrid, dimBlock>>> (MRows, MCols, MVals, MinvRows, MinvCols, MinvVals, nx, ny, dxD, dyD, dt);
	kernels::fillM_v <<<dimGrid, dimBlock>>> (MRows, MCols, MVals, MinvRows, MinvCols, MinvVals, nx, ny, dxD, dyD, dt);
}
