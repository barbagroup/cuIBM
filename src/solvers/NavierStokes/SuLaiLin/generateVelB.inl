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

#include <solvers/NavierStokes/kernels/generateVelB.h>

template <>
void SuLaiLinSolver<host_memory>::generateVelB()
{
	/// no-slip condition on the body surface
	for(int k=0; k<B.totalPoints; k++)
	{
		velB[k] = B.uB[k];
		velB[k+B.totalPoints] = B.vB[k];
	}
}

template<>
void SuLaiLinSolver<device_memory>::generateVelB()
{
	// raw pointers for cup arrays
	real *velB_r  = thrust::raw_pointer_cast(&(velB[0])),
	     *uB     = thrust::raw_pointer_cast(&(B.uB[0])),
	     *vB     = thrust::raw_pointer_cast(&(B.vB[0]));

	cusp::blas::fill(velB, 0.0);
	
	const int blockSize = 256;
	dim3 dimGrid( int((B.totalPoints-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	
	kernels::fill_velB <<<dimGrid, dimBlock>>> (velB_r, uB, vB, B.totalPoints);
}
