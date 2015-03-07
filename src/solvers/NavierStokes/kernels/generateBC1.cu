/***************************************************************************//**
 * \file generateBC1.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the kernels to generate right hand-side terms of the
 *        intermediate velocity flux solver.
 */


#include "generateBC1.h"


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

/**
 * \brief Computes inhomogeneous term from the discrete Laplacian operator 
 *        for the u-velocity at a given boundary with a Dirichlet-type condition.
 *
 * Each element is also normalized with the corresponding 
 * diagonal element of the matrix \f$ \hat{M} \f$.
 *
 * \param bc1 array that contains the boundary conditions
 * \param N number of u-velocity points on the boundary
 * \param nx number of cells in the x-direction
 * \param offset index in vector \c bc1 of the first element u-velocity point on the boundary
 * \param stride index increment
 * \param dx cell-widths in the x-direction
 * \param C coefficient from the Laplacian discretization at the boundary
 * \param bc boundary velocity
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
 * \brief Computes inhomogeneous term from the discrete Laplacian operator 
 *        for the v-velocity at a given boundary with a Dirichlet-type condition.
 *
 * Each element is also normalized with the corresponding 
 * diagonal element of the matrix \f$ \hat{M} \f$.
 *
 * \param bc1 array that contains the boundary conditions
 * \param N number of v-velocity points on the boundary
 * \param nx number of cells in the x-direction
 * \param numU number of u-velocity points in the domain
 * \param offset index in vector \c bc1 of the first element v-velocity point on the boundary
 * \param stride index increment
 * \param dy cell-widths in the y-direction
 * \param C coefficient from the Laplacian discretization at the boundary
 * \param bc boundary velocity
 * \param numUbc number of u-velocity points on the boundary
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
 * \brief Computes inhomogeneous term from the discrete Laplacian operator 
 *        for the u-velocity at a given boundary with a convective-type condition.
 *
 * Each element is also normalized with the corresponding 
 * diagonal element of the matrix \f$ \hat{M} \f$.
 *
 * \param bc1 array that contains the boundary conditions
 * \param N number of u-velocity points on the boundary
 * \param nx number of cells in the x-direction
 * \param offset index in vector \c bc1 of the first element u-velocity point on the boundary
 * \param stride index increment
 * \param dx cell-widths in the x-direction
 * \param dy cell-widths in the y-direction
 * \param C coefficient from the Laplacian discretization at the boundary
 * \param bc boundary velocity
 * \param q flux vector
 * \param beta linear interpolation coefficient (U*dt/dx)
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
 * \brief Computes inhomogeneous term from the discrete Laplacian operator 
 *        for the v-velocity at a given boundary with a convective-type condition.
 *
 * Each element is also normalized with the corresponding 
 * diagonal element of the matrix \f$ \hat{M} \f$.
 *
 * \param bc1 array that contains the boundary conditions
 * \param N number of v-velocity points on the boundary
 * \param nx number of cells in the x-direction
 * \param numU number of u-velocity points in the domain
 * \param offset index in vector \c bc1 of the first element v-velocity point on the boundary
 * \param stride index increment
 * \param dx cell-widths in the x-direction
 * \param dy cell-widths in the y-direction
 * \param C coefficient from the Laplacian discretization at the boundary
 * \param bc boundary velocity
 * \param numUbc number of u-velocity points on the boundary
 * \param q flux vector
 * \param beta linear interpolation coefficient (U*dt/dx)
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
 * \brief Computes inhomogeneous term from the discrete Laplacian operator 
 *        for the u-velocity at a given boundary with a special-type condition.
 *
 * It computes a sinusoidal displacement of the boundary in the x-direction.
 * Each element is also normalized with the corresponding 
 * diagonal element of the matrix \f$ \hat{M} \f$.
 *
 * \param bc1 array that contains the boundary conditions
 * \param N number of v-velocity points on the boundary
 * \param nx number of cells in the x-direction
 * \param offset index in vector \c bc1 of the first element v-velocity point on the boundary
 * \param stride index increment
 * \param dx cell-widths in the x-direction
 * \param C coefficient from the Laplacian discretization at the boundary
 * \param bc boundary velocity
 * \param time the current time
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
