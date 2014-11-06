/***************************************************************************//**
 * \file generateBC2.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods to generate 
 *        inhomogeneous boundary conditions from the discrete continuity equation.
 */


#include <solvers/NavierStokes/kernels/generateBC2.h>


/**
 * \brief Generates inhomogeneous boundary conditions 
 *        from the discrete continuity equation (on the host).
 */
template <>
void NavierStokesSolver<host_memory>::generateBC2()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	//int  numP  = nx*ny;
	
	real *dx = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	// bc2 and rhs2 are of size np+2nb
	
	// the boundary conditions are obtained from the continuity equations
	
	cusp::blas::fill(bc2, 0.0);  // initialise the elements of bc2 to zero
	
	for(int i=0; i<nx; i++)
	{
		// bottom
		bc2[i] -= bc[YMINUS][i+nx-1]*dx[i]; // v[0][i+1]*dx;
		// top
		bc2[(ny-1)*nx + i] += bc[YPLUS][i+nx-1]*dx[i]; // v[ny][i+1]*dx;
	}

	for(int j=0; j<ny; j++)
	{
		// left
		bc2[j*nx] -= bc[XMINUS][j]*dy[j]; // u[j+1][0]*dy;
		// right
		bc2[j*nx+nx-1] += bc[XPLUS][j]*dy[j]; // u[j+1][nx]*dy;
	}
}

/**
 * \brief Generates inhomogeneous boundary conditions 
 *        from the discrete continuity equation (on the device).
 */
template<>
void NavierStokesSolver<device_memory>::generateBC2()
{
	// raw pointers for cup arrays
	real *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	real *bc2_r = thrust::raw_pointer_cast(&(bc2[0])),
	     *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0]));
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	// initialize array with zeros
	cusp::blas::fill(bc2, 0.0);
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx+ny-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	
	kernels::fillBC2_v <<<dimGrid, dimBlock>>> (bc2_r, yminus, yplus, dxD, nx, ny);
	kernels::fillBC2_u <<<dimGrid, dimBlock>>> (bc2_r, xminus, xplus, dyD, nx, ny);
}
