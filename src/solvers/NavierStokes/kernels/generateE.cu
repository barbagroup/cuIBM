/***************************************************************************//**
* \file generateE.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the kernels required to generate elements of the interpolation matrix
*/

#include <solvers/NavierStokes/kernels/generateE.h>

/**
* \brief Discrete delta function defined by Roma et al. (1999)
*
* \param x x- or y- component of the vector defined between a body point and a grid point
* \param h grid-spacing near the body point
*
* \return the value of the discrete Delta function
*/
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

/**
* \brief Compute the discrete Delta function in 2D 
*		using the definition by Roma et al. (1999)
*
* \param x x-component of the vector defined between a body point and a grid point
* \param y y-component of the vector defined between a body point and a grid point
* \param h grid-spacing near the body point
*
* \return the discrete Delta function in 2D
*/
__device__ \
real deltaDeviceE(real x, real y, real h)
{
	return dhRomaDeviceE(x, h) * dhRomaDeviceE(y, h);
}

/**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{
	
/**
* \brief Compute elements of the interpolation matrix on the host
*
* \param ERows rows of the COO matrix \c E
* \param ECols columns of the COO matrix \c E
* \param EVals values of the COO matrix \c E
* \param nx number of cells in the x-direction
* \param ny number of cells in the y-direction
* \param x x-component of the Eulerian grid points
* \param y y-component of the Eulerian grid points
* \param dx cell-widths in the x-direction
* \param totalPoints total number of body points (all bodies included)
* \param xB x-coordinate of the body points
* \param yB y-coordinate of the body points
* \param I x-index of the cell in which the body point is located
* \param J y-index of the cell in which the body point is located
*/
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

/**
* \brief Kernel to compute an element of the interpolation matrix on the device
*
* \param ERows rows of the COO matrix \c E
* \param ECols columns of the COO matrix \c E
* \param EVals values of the COO matrix \c E
* \param nx number of cells in the x-direction
* \param ny number of cells in the y-direction
* \param x x-component of the Eulerian grid points
* \param y y-component of the Eulerian grid points
* \param dx cell-widths in the x-direction
* \param totalPoints total number of body points (all bodies included)
* \param xB x-coordinate of the body points
* \param yB y-coordinate of the body points
* \param I x-index of the cell in which the body point is located
* \param J y-index of the cell in which the body point is located
*/
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
