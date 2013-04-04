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

#include <solvers/NavierStokes/kernels/generateE.h>

template <>
void SuLaiLinSolver<host_memory>::generateE()
{
}

template <>
void SuLaiLinSolver<device_memory>::generateE()
{
	logger.startTimer("generateE");
	
	int  *ERows = thrust::raw_pointer_cast(&(E.row_indices[0])),
	     *ECols = thrust::raw_pointer_cast(&(E.column_indices[0]));
	real *EVals = thrust::raw_pointer_cast(&(E.values[0]));
	
	real *x  = thrust::raw_pointer_cast(&(domInfo->xD[0])),
	     *y  = thrust::raw_pointer_cast(&(domInfo->yD[0])),
	     *dx = thrust::raw_pointer_cast(&(domInfo->dxD[0]));
	
	real *xB = thrust::raw_pointer_cast(&(B.x[0])),
	     *yB = thrust::raw_pointer_cast(&(B.y[0]));
	
	int *I = thrust::raw_pointer_cast(&(B.I[0])),
	    *J = thrust::raw_pointer_cast(&(B.J[0]));
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	     
	const int blockSize = 256;
	     
	dim3 dimGrid(int((B.totalPoints-0.5)/blockSize)+1, 1);
	dim3 dimBlock(blockSize, 1);

	kernels::generateE <<<dimGrid, dimBlock>>> (ERows,  ECols,  EVals, nx, ny, x, y, dx, B.totalPoints, xB, yB, I, J);
    
    logger.stopTimer("generateE");
	
	logger.startTimer("transposeE");
	cusp::transpose(E, ET);
	logger.stopTimer("transposeE");
}
