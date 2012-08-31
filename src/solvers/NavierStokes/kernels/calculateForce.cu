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

namespace kernels
{

__global__
void dragLeftRight(real *FxY, real *q, real *lambda, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx >= ncy)
		return;
	int  Ip = (J+idx)*nx + I,
	     Iu = (J+idx)*(nx-1) + (I-1);
	FxY[idx] = -(
	              // multiply the pressure with the surface area to get p dy
	              (lambda[Ip+ncx]-lambda[Ip-1])*dy[J+idx]
	              +
	              // divide q^2 by dy, so that just u^2 dy is obtained
	              (
	                0.25*(q[Iu+ncx+1] + q[Iu+ncx])*(q[Iu+ncx+1] + q[Iu+ncx]) -
	                0.25*(q[Iu] + q[Iu-1])*(q[Iu] + q[Iu-1])
	              )/dy[J+idx]
	              -
	              // no multiplication or division since du/dx dy = dq/dx
	              (
	                2*(q[Iu+ncx+1] - q[Iu+ncx])/dx[I+ncx] -
	                2*(q[Iu] - q[Iu-1])/dx[I-1]
	              )
	            );	
}

__global__
void dragBottomTop(real *FxX, real *q, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx >= ncx)
		return;
	int  Iu = J*(nx-1) + (I+idx-1),
	     Iv = (nx-1)*ny + (J-1)*nx + I+idx;
	FxX[idx] = -(
	              // multiply by dS
	              (
	                0.25 * (q[Iu+ncy*(nx-1)]/dy[J+ncy] + q[Iu+(ncy-1)*(nx-1)]/dy[J+ncy-1])
	                     * (q[Iv+ncy*nx]/dx[I+idx] + q[Iv+ncy*nx-1]/dx[I+idx-1])
	                -
	                0.25 * (q[Iu]/dy[J] + q[Iu-(nx-1)]/dy[J-1]) * (q[Iv]/dx[I+idx] + q[Iv-1]/dx[I+idx-1])
	              )
	              -
	              // multiply by dS (cannot use the leftRight trick in this case)
	              (
	                (
	                  (q[Iu+ncy*(nx-1)]/dy[J+ncy] - q[Iu+(ncy-1)*(nx-1)]/dy[J+ncy-1])/2.0/(dy[J+ncy]+dy[J+ncy-1]) +
	                  (q[Iv+ncy*nx]/dx[I+idx] - q[Iv+ncy*nx-1]/dx[I+idx-1])/2.0/(dx[I+idx]+dx[I+idx-1])
	                ) -
	                (
	                  (q[Iu]/dy[J] - q[Iu-(nx-1)]/dy[J-1])/2.0/(dy[J]+dy[J-1]) +
	                  (q[Iv]/dx[I+idx] - q[Iv-1]/dx[I+idx-1])/2.0/(dx[I+idx]+dx[I+idx-1])
	                )
	              )
	            )*0.5*(dx[I+idx]+dx[I+idx-1]);	

}

__global__
void dragUnsteady(real *FxU, real *q, real *qOld, real *dx, real *dy, real dt,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	     
	if(idx >= (ncx+1)*ncy)
		return;
	
	int i = idx%(ncx+1),
	    j = idx/(ncx+1);
	    
	int Iu = (J+j)*(nx-1) + (I+i-1);
	
	FxU[idx] = - ( q[Iu] - qOld[Iu])/dt * 0.5*(dx[I+i]+dx[I+i-1]);
}

} // end of namespace kernels
