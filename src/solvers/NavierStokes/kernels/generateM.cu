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

#include <solvers/NavierStokes/kernels/generateM.h>

namespace kernels
{

__global__
void fillM_u(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= (nx-1)*ny) return;
	
	int  i = I % (nx-1);
	int  j = I / (nx-1);
	real value = 0.5*(dx[i]+dx[i+1])/dy[j]/dt;
	
	MRows[I] = I;
	MCols[I] = I;
	MVals[I] = value;
	
	MinvRows[I] = I;
	MinvCols[I] = I;
	MinvVals[I] = 1.0/value;
}

__global__
void fillM_v(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= nx*(ny-1)) return;
	
	int  numU = (nx-1)*ny;
	int  i = I % nx;
	int  j = I / nx;
	real value = 0.5*(dy[j]+dx[j+1])/dx[i]/dt;
	
	MRows[I+numU] = I+numU;
	MCols[I+numU] = I+numU;
	MVals[I+numU] = value;
	
	MinvRows[I+numU] = I+numU;
	MinvCols[I+numU] = I+numU;
	MinvVals[I+numU] = 1.0/value;
}

} // end of namespace kernels