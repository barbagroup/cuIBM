#include <solvers/NavierStokes/kernels/generateBC1.h>

namespace kernels
{

__global__
void bc1DirichletU(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
	
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I	= (offset + idx*stride),
		i	= I % (nx-1);
	bc1[I] += bc[idx] * C * 0.5*(dx[i] + dx[i+1]);
}

__global__
void bc1DirichletV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dy, real C, real *bc, int numUbc)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
		
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I	= (offset + idx*stride),
		j	= I / nx;
	bc1[numU + I] += bc[idx+numUbc] * C * 0.5*(dy[j] + dy[j+1]);
}

__global__
void bc1ConvectiveU(real *bc1, int N, int nx, int offset, int stride, real *dx, real *dy, real C, real *bc, real *q, real alpha)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
	
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I = (offset + idx*stride),
		i = I % (nx-1),
		j = I / (nx-1);
	
	bc[idx] = (1.0-alpha)*bc[idx] + alpha*q[I]/dy[j];
	
	bc1[I] += bc[idx] * C * 0.5*(dx[i] + dx[i+1]);
}

__global__
void bc1ConvectiveV(real *bc1, int N, int nx, int numU, int offset, int stride, real *dx, real *dy, real C, real *bc, int numUbc, real *q, real alpha)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
		
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I = (offset + idx*stride),
		i = I % nx,
		j = I / nx;
	
	bc[idx+numUbc] = (1.0-alpha)*bc[idx+numUbc] + alpha*q[numU + I]/dx[i];
	
	bc1[numU + I] += bc[idx+numUbc] * C * 0.5*(dy[j] + dy[j+1]);
}

} // end of namespace kernels