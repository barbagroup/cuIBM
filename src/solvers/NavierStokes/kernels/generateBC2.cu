/***************************************************************************//**
* \file generateBC2.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the kernels required to generate vector bc2
*/

#include <solvers/NavierStokes/kernels/generateBC2.h>

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
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=nx)
		return;
	bc2[i] -= yminus[i+nx-1]*dx[i];
	bc2[(ny-1)*nx + i] += yplus[i+nx-1]*dx[i];
}

/**
* \brief To be documented
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
* \brief To be documented
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
