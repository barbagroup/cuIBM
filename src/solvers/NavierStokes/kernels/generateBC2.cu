#include <solvers/NavierStokes/kernels/generateBC2.h>
namespace kernels
{

__global__
void fillBC2_v(real *bc2, real *yminus, real *yplus, real *dx, int nx, int ny)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=nx)
		return;
	bc2[i] -= yminus[i+nx-1]*dx[i];
	bc2[(ny-1)*nx + i] += yplus[i+nx-1]*dx[i];
}

__global__
void fillBC2_u(real *bc2, real *xminus, real *xplus, real *dy, int nx, int ny)
{
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	if(j>=ny)
		return;
	bc2[j*nx] -= xminus[j]*dy[j];
	bc2[j*nx+nx-1] += xplus[j]*dy[j];
}

__global__
void fillBC2_uvB(real *bc2, real *uB, real *vB, int totalPoints, int nx, int ny)
{
	int k = threadIdx.x + blockIdx.x*blockDim.x;
	if(k>=totalPoints)
		return;
	bc2[nx*ny + k] = uB[k];
	bc2[nx*ny + k + totalPoints] = vB[k];
}

}