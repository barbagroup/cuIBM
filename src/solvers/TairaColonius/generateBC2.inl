/**
 * \file generateBC2.inl
 * \brief Implementation of the methods of the class \c TairaColoniusSolver
 *        to generate the right hand-side of the Poisson system.
 */


#include <solvers/kernels/generateBC2.h>


/**
 * \brief Generates the right hand-side of the Poisson system.
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
	
	int nx = domInfo->nx,
	    ny = domInfo->ny;
	
	cusp::blas::fill(bc2, 0.0);
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx+ny-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	kernels::fillBC2_v   <<<dimGrid, dimBlock>>> (bc2_r, yminus, yplus, dxD, nx, ny);
	kernels::fillBC2_u   <<<dimGrid, dimBlock>>> (bc2_r, xminus, xplus, dyD, nx, ny);
	
	dim3 dimGridB( int((B.totalPoints-0.5)/blockSize) + 1, 1);
	kernels::fillBC2_uvB <<<dimGridB, dimBlock>>> (bc2_r, uB, vB, B.totalPoints, nx, ny);
} // generateBC2
