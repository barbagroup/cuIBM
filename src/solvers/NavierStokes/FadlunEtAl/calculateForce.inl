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

#include <solvers/NavierStokes/kernels/calculateForce.h>

#define BSZ 16

template <>
void FadlunEtAlSolver<host_memory>::calculateForceF()
{
}

template <>
void FadlunEtAlSolver<device_memory>::calculateForceF()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	     
	real alpha = intgSchm.alphaImplicit[subStep],
	     nu = (*paramDB)["flow"]["nu"].get<real>();
	     
	cusp::array1d<real, device_memory>
	     f((nx-1)*ny + nx*(ny-1), 0.0),
	     temp((nx-1)*ny + nx*(ny-1));
	     
	// raw pointers for cup arrays
	real *f_r  = thrust::raw_pointer_cast(&f[0]),
	     *q_r  = thrust::raw_pointer_cast(&q[0]),
	     *rn_r = thrust::raw_pointer_cast(&rn[0]),
	     *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	int  *tagsX_r = thrust::raw_pointer_cast(&(tagsXD[0])),\
	     *tagsY_r = thrust::raw_pointer_cast(&(tagsYD[0]));
	
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	real h = 1;//domInfo->dx[ B.I[0] ];
	     
	dim3 dimGridx( int( (nx-1-0.5)/(BSZ-2) ) + 1, int( (ny-0.5)/(BSZ-2) ) + 1 ),
	     dimGridy( int( (nx-0.5)/(BSZ-2) ) + 1, int( (ny-1-0.5)/(BSZ-2) ) + 1 );
	dim3 dimBlock(BSZ, BSZ);
	
	// call the kernel
	kernels::forceX <<<dimGridx, dimBlock>>> (f_r, q_r, rn_r, tagsX_r, tagsY_r, nx, ny, dxD, dyD, dt, alpha, nu);
	
	//cusp::multiply(Q, lambda, temp);
	//cusp::blas::axpy(temp, f, 1.0);
	
	forceX = (h*h)*thrust::reduce( f.begin(), f.begin()+(nx-1)*ny );
	forceY = (h*h)*thrust::reduce( f.begin()+(nx-1)*ny, f.end() );
}
