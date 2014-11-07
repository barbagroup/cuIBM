/***************************************************************************//**
 * \file generateBC2.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the kernels to generate elements of the right hand-side
 *        of the Poisson solver.
 */


#include "generateBC2.h"


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

/**
 * \brief Computes inhomogeneous terms of the discrete divergence operator 
 *        from the bottom and top boundaries at the v-velocity locations.
 *
 * \param bc2 array that contains boundary conditions
 * \param yminus bottom-boundary velocities
 * \param yplus top-boundary velocities
 * \param dx cell-widths in the x-direction
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 */
__global__
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=nx)
		return;
	bc2[i] -= yminus[i+nx-1]*dx[i];
	bc2[(ny-1)*nx + i] += yplus[i+nx-1]*dx[i];
}

/**
 * \brief Computes inhomogeneous terms of the discrete divergence operator 
 *        from the left and right boundaries at the u-velocity locations.
 *
 * \param bc2 array that contains boundary conditions
 * \param xminus left-boundary velocities
 * \param xplus right-boundary velocities
 * \param dy cell-widths in the x-direction
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 */
__global__
void fillBC2_u(real *bc2, real *xminus, real *xplus, real *dy, int nx, int ny)
{
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	if(j>=ny)
		return;
	bc2[j*nx] -= xminus[j]*dy[j];
	bc2[j*nx+nx-1] += xplus[j]*dy[j];
}

/**
 * \brief Computes inhomogeneous terms of the discrete divergence operator 
 *        from the no-slip constraint at the body-point locations.
 *
 * \param bc2 array that contains boundary conditions
 * \param uB x-component of the body-velocity
 * \param vB y-component of the body-velcoity
 * \param totalPoints number of body-points (all bodies included)
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 */
__global__
void fillBC2_uvB(real *bc2, real *uB, real *vB, int totalPoints, int nx, int ny)
{
	int k = threadIdx.x + blockIdx.x*blockDim.x;
	if(k>=totalPoints)
		return;
	bc2[nx*ny + k] = uB[k];
	bc2[nx*ny + k + totalPoints] = vB[k];
}

} // end of namespace kernels
