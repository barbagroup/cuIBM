/***************************************************************************//**
 * \file generateBC2.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c TairaColoniusSolver
 *        to generate the right hand-side of the Poisson system.
 */


#include <solvers/NavierStokes/kernels/generateBC2.h>


/**
 * \brief Generates the right hand-side of the Poisson system (on the host).
 *
 * It contains the inhomogeneous boundary conditions from the discrete divergence
 * operator, as well as the no-slip boundary condition at the body surface.
 *
 */
template <>
void TairaColoniusSolver<host_memory>::generateBC2()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	real *dx = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	// rhs2 is of size np+2nb	
	// boundary conditions are obtained from the continuity equations
		
	// initialise the values in bc2 to zero
	cusp::blas::fill(bc2, 0.0);
	
	for(int i=0; i<nx; i++)
	{
		bc2[i] -= bc[YMINUS][i+nx-1]*dx[i]; // v[0][i+1]*dx;
		bc2[(ny-1)*nx + i] += bc[YPLUS][i+nx-1]*dx[i]; // v[ny][i+1]*dx;
	}

	for(int j=0; j<ny; j++)
	{
		// left
		bc2[j*nx] -= bc[XMINUS][j]*dy[j]; // u[j+1][0]*dy;
		// right
		bc2[j*nx+nx-1] += bc[XPLUS][j]*dy[j]; // u[j+1][nx]*dy;
	}	

	int	row = nx*ny;
	/// no-slip condition on the body surface
	for(int k=0; k<B.totalPoints; k++)
	{
		bc2[row] = B.uB[k];
		row++;
	}
	for(int k=0; k<B.totalPoints; k++)
	{
		bc2[row] = B.vB[k];
		row++;
	}
}

/**
 * \brief Generates the right hand-side of the Poisson system (on the device).
 *
 * It contains the inhomogeneous boundary conditions from the discrete divergence
 * operator, as well as the no-slip boundary condition at the body surface.
 *
 */
template<>
void TairaColoniusSolver<device_memory>::generateBC2()
{
	// raw pointers for cup arrays
	real *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	real *bc2_r  = thrust::raw_pointer_cast(&(bc2[0])),
	     *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0])),
	     *uB     = thrust::raw_pointer_cast(&(B.uB[0])),
	     *vB     = thrust::raw_pointer_cast(&(B.vB[0]));
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	cusp::blas::fill(bc2, 0.0);
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx+ny-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	kernels::fillBC2_v   <<<dimGrid, dimBlock>>> (bc2_r, yminus, yplus, dxD, nx, ny);
	kernels::fillBC2_u   <<<dimGrid, dimBlock>>> (bc2_r, xminus, xplus, dyD, nx, ny);
	
	dim3 dimGridB( int((B.totalPoints-0.5)/blockSize) + 1, 1);
	kernels::fillBC2_uvB <<<dimGridB, dimBlock>>> (bc2_r, uB, vB, B.totalPoints, nx, ny);
}
