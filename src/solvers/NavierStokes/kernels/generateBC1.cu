/***************************************************************************//**
* \file generateBC1.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the kernels required to generate vector bc1
*/

#include <solvers/NavierStokes/kernels/generateBC1.h>

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
__global__
void bc1DirichletU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
	
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I	= (offset + idx*stride),
		i	= I % (nx-1);
	bc1[I] += bc[idx] * C * 0.5*(dx[i] + dx[i+1]);
}

/**
* \brief To be documented
*/
__global__
void bc1DirichletV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dy, real C, real *bc, int numUbc)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
		
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I	= (offset + idx*stride),
		j	= I / nx;
	bc1[numU + I] += bc[idx+numUbc] * C * 0.5*(dy[j] + dy[j+1]);
}

/**
* \brief To be documented
*/
__global__
void bc1ConvectiveU(real *bc1, int N, int nx, int offset, int stride, real *dx, real *dy, real C, real *bc, real *q, real beta)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
	
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I = (offset + idx*stride),
		i = I % (nx-1),
		j = I / (nx-1);
	
	bc[idx] = (1.0-beta)*bc[idx] + beta*q[I]/dy[j];
	
	bc1[I] += bc[idx] * C * 0.5*(dx[i] + dx[i+1]);
}

/**
* \brief To be documented
*/
__global__
void bc1ConvectiveV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dx, real *dy, real C, real *bc, int numUbc, real *q, real beta)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
		
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I = (offset + idx*stride),
		i = I % nx,
		j = I / nx;
	
	bc[idx+numUbc] = (1.0-beta)*bc[idx+numUbc] + beta*q[numU + I]/dx[i];
	
	bc1[numU + I] += bc[idx+numUbc] * C * 0.5*(dy[j] + dy[j+1]);
}

/**
* \brief To be documented
*/
__global__
void bc1SpecialU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc, real time)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
	
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I = (offset + idx*stride),
		i = I % (nx-1);
		
	const real T = 10.0;
	
	bc[idx] = sin(M_PI*time/T);
	
	bc1[I] += bc[idx] * C * 0.5*(dx[i] + dx[i+1]);
}

} // end of namespace kernels
