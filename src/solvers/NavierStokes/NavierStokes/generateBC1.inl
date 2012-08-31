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

#include <solvers/NavierStokes/kernels/generateBC1.h>

template <>
void NavierStokesSolver<device_memory>::generateBC1Full(real alpha)
{
	// raw pointers from cusp arrays
	real *bc1_r = thrust::raw_pointer_cast(&bc1[0]),
	     *q_r   = thrust::raw_pointer_cast(&q[0]),
	     *dx  = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy  = thrust::raw_pointer_cast(&(domInfo->dy[0])),
	     *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	real *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0]));

	real dx0, dx1, dy0, dy1,
	     nu = (*paramDB)["flow"]["nu"].get<real>(),
         dt = (*paramDB)["simulation"]["dt"].get<real>();

	boundaryCondition **bcInfo = (*paramDB)["flow"]["boundaryConditions"].get<boundaryCondition **>();
	
	const int blocksize = 256;
	
	real C, beta;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	
	// zero the bc1 vector
	cusp::blas::fill(bc1, 0.0);
	
	dim3 dimGridx( int((nx - 0.5)/blocksize) + 1, 1),
	     dimGridy( int((ny - 0.5)/blocksize) + 1, 1);
			
	dim3 dimBlock(blocksize, 1);
	
	/// bottom
	if(bcInfo[YMINUS][0].type == DIRICHLET)
	{
		/// u
		dy0	= 0.5*dy[0];
		dy1	= 0.5*(dy[0] + dy[1]);
		/// multiply by 0.5 for the Crank-Nicolson scheme and 2.0 for the non-uniform central difference
		C	= alpha * 2.0 * nu / (dy0 * (dy0+dy1));
		kernels::bc1DirichletU <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, 0, 1, dxD, C, yminus);
		
		/// v
		C	= alpha * 2.0 * nu / (dy[0] * (dy[0]+dy[1]));
		kernels::bc1DirichletV <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, 0, 1, dyD, C, yminus, nx-1);
	}
	//else if(F.nbc.bottom_type==BC_CONVECTIVE)
	
	/// top
	if(bcInfo[YPLUS][0].type == DIRICHLET)
	{
		/// u
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= alpha * 2.0 * nu / (dy1 * (dy0+dy1));
		kernels::bc1DirichletU <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, (nx-1)*(ny-1), 1, dxD, C, yplus);
		
		/// v
		C	= alpha * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		kernels::bc1DirichletV <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, nx*(ny-2), 1, dyD, C, yplus, nx-1);
	}
	else if(bcInfo[YPLUS][0].type == CONVECTIVE)
	{
		/// u
		beta = bcInfo[YPLUS][1].value * dt / (0.5 * dy[ny-1]);
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= alpha * 2.0 * nu / (dy1 * (dy0+dy1));
		kernels::bc1ConvectiveU <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, (nx-1)*(ny-1), 1, dxD, dyD, C, yplus, q_r, beta);
		
		/// v
		beta = bcInfo[YPLUS][1].value * dt / dy[ny-1];
		C	= alpha * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		kernels::bc1ConvectiveV <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, nx*(ny-2), 1, dxD, dyD, C, yplus, nx-1, q_r, beta);
	}
	
	/// left
	if(bcInfo[XMINUS][0].type == DIRICHLET)
	{
		/// u
		C = alpha * 2.0 * nu / ( dx[0] * (dx[0]+dx[1]) );
		kernels::bc1DirichletU <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, 0, (nx-1), dxD, C, xminus);
		
		/// v
		dx0	= 0.5*dx[0];
		dx1	= 0.5*(dx[0] + dx[1]);
		C = alpha * 2.0 * nu / (dx0 * (dx0+dx1));
		kernels::bc1DirichletV <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, 0, nx, dyD, C, xminus, ny);
	}
	//else if(F.nbc.left_type==BC_CONVECTIVE)

	/// right
	if(bcInfo[XPLUS][0].type == DIRICHLET)
	{
		/// u
		C = alpha * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		kernels::bc1DirichletU <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, nx-2, (nx-1), dxD, C, xplus);
		
		/// v
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= alpha * 2.0 * nu / (dx1 * (dx0+dx1));
		kernels::bc1DirichletV <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, nx-1, nx, dyD, C, xplus, ny);
	}
	else if(bcInfo[XPLUS][0].type == CONVECTIVE)
	{
		/// u
		beta = bcInfo[XPLUS][0].value * dt / dx[nx-1];
		C = alpha * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		kernels::bc1ConvectiveU <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, nx-2, nx-1, dxD, dyD, C, xplus, q_r, beta);
		
		/// v
		beta = bcInfo[XPLUS][0].value * dt / (0.5*dx[nx-1]);
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= alpha * 2.0 * nu / (dx1 * (dx0+dx1));
		kernels::bc1ConvectiveV <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, nx-1, nx, dxD, dyD, C, xplus, ny, q_r, beta);
	}
}

template <>
void NavierStokesSolver<host_memory>::generateBC1Full(real alpha)
{
}
