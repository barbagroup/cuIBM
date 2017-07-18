/**
 * \file generateBC2.inl
 * \brief Implementation of the methods to generate 
 *        inhomogeneous boundary conditions from the discrete continuity equation.
 */


#include <solvers/kernels/generateBC2.h>


/**
 * \brief Generates inhomogeneous boundary conditions 
 *        from the discrete continuity equation.
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
	
	int nx = domInfo->nx,
	    ny = domInfo->ny;

	// initialize array with zeros
	cusp::blas::fill(bc2, 0.0);
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx+ny-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	
	kernels::fillBC2_v <<<dimGrid, dimBlock>>> (bc2_r, yminus, yplus, dxD, nx, ny);
	kernels::fillBC2_u <<<dimGrid, dimBlock>>> (bc2_r, xminus, xplus, dyD, nx, ny);
} // generateBC2
