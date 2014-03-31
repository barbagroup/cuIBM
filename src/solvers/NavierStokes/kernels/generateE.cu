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

__device__ \
real dhRomaDeviceE(real x, real h)
{
	real r = fabs(x)/h;
	
	if(r>1.5)
		return 0.0;
	else if(r>0.5 && r<=1.5)
		return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
	else
		return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
}

__device__ \
real deltaDeviceE(real x, real y, real h)
{
	return dhRomaDeviceE(x, h) * dhRomaDeviceE(y, h);
}

namespace kernels
{
	
void generateEHost(int *ERows,  int *ECols,  real *EVals,
                   int nx, int ny, real *x, real *y, real *dx,
                   int totalPoints, real *xB, real *yB, int *I, int *J)
{
	for(int bodyIdx=0; bodyIdx<totalPoints; bodyIdx++)
	{
		int  Ib=I[bodyIdx],
		     Jb=J[bodyIdx],
		     EIdx  = bodyIdx*12,
		     i, j;

		real Dx = dx[Ib];
	
		// uB = integral (u * delta * dxdy)
		// E = E_hat * R^-1 => divide E_hat by Dx
	
		// populate x-components
		for(j=Jb-1; j<=Jb+1; j++)
		{
			for(i=Ib-2; i<=Ib+1; i++)
			{
				ERows[EIdx] = bodyIdx;
				ECols[EIdx] = j*(nx-1) + i;
				EVals[EIdx] = Dx*delta(x[i+1]-xB[bodyIdx], 0.5*(y[j]+y[j+1])-yB[bodyIdx], Dx);
				EIdx++;
			}
		}

		// populate y-components
		for(j=Jb-2; j<=Jb+1; j++)
		{
			for(i=Ib-1; i<=Ib+1; i++)
			{
				ERows[EIdx+12*totalPoints-12] = bodyIdx + totalPoints;
				ECols[EIdx+12*totalPoints-12] = j*nx + i + (nx-1)*ny;
				EVals[EIdx+12*totalPoints-12] = Dx*delta(0.5*(x[i]+x[i+1])-xB[bodyIdx], y[j+1]-yB[bodyIdx], Dx);
				EIdx++;
			}
		}
	}
}

__global__ \
void generateE(int *ERows,  int *ECols,  real *EVals,
               int nx, int ny, real *x, real *y, real *dx,
               int totalPoints, real *xB, real *yB, int *I, int *J)
{
	int bodyIdx = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(bodyIdx < totalPoints)
	{
		int  Ib=I[bodyIdx],
			 Jb=J[bodyIdx],
			 EIdx  = bodyIdx*12,
			 i, j;

		real Dx = dx[Ib];
	
		// uB = integral u * delta * dxdy = Ehat * u
		// E = Ehat * R^-1 => divide by Dx
		// E = Dx * delta
	
		// populate x-components
		for(j=Jb-1; j<=Jb+1; j++)
		{
			for(i=Ib-2; i<=Ib+1; i++)
			{
				ERows[EIdx] = bodyIdx;
				ECols[EIdx] = j*(nx-1) + i;
				EVals[EIdx] = Dx*deltaDeviceE(x[i+1]-xB[bodyIdx], 0.5*(y[j]+y[j+1])-yB[bodyIdx], Dx);
				EIdx++;
			}
		}

		// populate y-components
		for(j=Jb-2; j<=Jb+1; j++)
		{
			for(i=Ib-1; i<=Ib+1; i++)
			{
				ERows[EIdx+12*totalPoints-12] = bodyIdx + totalPoints;
				ECols[EIdx+12*totalPoints-12] = j*nx + i + (nx-1)*ny;
				EVals[EIdx+12*totalPoints-12] = Dx*deltaDeviceE(0.5*(x[i]+x[i+1])-xB[bodyIdx], y[j+1]-yB[bodyIdx], Dx);
				EIdx++;
			}
		}
	}
}

} // end of namespace kernels
