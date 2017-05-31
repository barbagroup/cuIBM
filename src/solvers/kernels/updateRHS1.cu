/**
 * \file updateRHS1.cu
 * \brief Implementation of the kernels to update the right hand-side
 *        of the intermediate velocity flux solver.
 *        It replaces the right-hand side elements at the tagged points, with
 *        values obtained from the interpolation relations at those points.
 */


#include "updateRHS1.h"


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

/**
 * \brief Update the RHS vector of the velocity system at forcing nodes.
 *
 * Elements associated with a forcing node are set to 0.
 *
 * \param rhs1 vector to update
 * \param numUV number of velocity points in the domain
 * \param tags vector used to differentiate regular points from forcing points
 */
__global__
void updateRHS1(real *rhs1, int numUV, int *tags)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I>=numUV)
		return;
	
	rhs1[I] = rhs1[I]*(tags[I]==-1);
} // updateRHS1


/**
 * \brief Perform 1D linear interpolation at forcing points for the x-velocity.
 *
 * \param rhs1 vector to update
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param dt time-step size
 * \param dx grid-spacings along a gridline in the x-direction
 * \param tags vector used to differentiate regular points from forcing points
 * \param coeffs coefficients of interpolation
 * \param uv velocity vector
 */
__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *coeffs, real *uv)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	int i = I % (nx-1);
	
	if( I < (nx-1)*ny )
	{
		rhs1[I] = (tags[I]==-1)*rhs1[I] 
		          + ((tags[I]!=-1)*((1.0-coeffs[I])*uv[I])) * 0.5*(dx[i+1]+dx[i])/dt;
	}
} // updateRHS1X


/**
 * \brief Perform 1D linear interpolation at forcing points for the y-velocity.
 *
 * \param rhs1 vector to update
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param dt time-step size
 * \param dy grid-spacings along a gridline in the y-direction
 * \param tags vector used to differentiate regular points from forcing points
 * \param coeffs coefficients of interpolation
 * \param uv velocity vector
 */
__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *coeffs, real *uv)
{
	int numU = (nx-1)*ny;
	int	I = blockIdx.x*blockDim.x + threadIdx.x + numU;
	int j = (I-numU) / nx;
	
	if( I < numU + nx*(ny-1) )
	{
		rhs1[I] = (tags[I]==-1)*rhs1[I] 
		           + ((tags[I]!=-1)*((1.0-coeffs[I])*uv[I])) * 0.5*(dy[j+1]+dy[j])/dt;
	}
} // updateRHS1Y


/**
 * \brief Perform 1D quadratic interpolation at forcing points for the x-velocity.
 *
 * \param rhs1 vector to update
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param dt time-step size
 * \param dx grid-spacings along a gridline in the x-direction
 * \param tags vector used to differentiate regular points from forcing points
 * \param coeffs coefficients of interpolation
 * \param coeffs2 coefficients of interpolation
 * \param uv velocity vector
 */
__global__
void updateRHS1X(real *rhs1, int nx, int ny, real dt, real *dx, int *tags, real *coeffs, real *coeffs2, real *uv)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	int i = I % (nx-1);
	
	if( I < (nx-1)*ny )
	{
		rhs1[I] = (tags[I]==-1)*rhs1[I] 
		           + ((tags[I]!=-1)*((1.0-coeffs[I]-coeffs2[I])*uv[I])) * 0.5*(dx[i+1]+dx[i])/dt;
	}
} // updateRHS1X


/**
 * \brief Perform 1D quadratic interpolation at forcing points for the y-velocity.
 *
 * \param rhs1 vector to update
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param dt time-step size
 * \param dy grid-spacings along a gridline in the y-direction
 * \param tags vector used to differentiate regular points from forcing points
 * \param coeffs coefficients of interpolation
 * \param coeffs2 coefficients of interpolation
 * \param uv velocity vector
 */
__global__
void updateRHS1Y(real *rhs1, int nx, int ny, real dt, real *dy, int *tags, real *coeffs, real *coeffs2, real *uv)
{
	int numU = (nx-1)*ny;
	int	I = blockIdx.x*blockDim.x + threadIdx.x + numU;
	int j = (I-numU) / nx;
	
	if( I < numU + nx*(ny-1) )
	{
		rhs1[I] = (tags[I]==-1)*rhs1[I] 
		           + ((tags[I]!=-1)*((1.0-coeffs[I]-coeffs2[I])*uv[I])) * 0.5*(dy[j+1]+dy[j])/dt;
	}
} // updateRHS1Y

} // End of namespace kernels
