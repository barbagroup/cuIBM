#include <solvers/NavierStokes/kernels/generateVelB.h>

template <>
void SuLaiLinSolver<host_memory>::generateVelB()
{
	/// no-slip condition on the body surface
	for(int k=0; k<B.totalPoints; k++)
	{
		velB[k] = B.uB[k];
		velB[k+B.totalPoints] = B.vB[k];
	}
}

template<>
void SuLaiLinSolver<device_memory>::generateVelB()
{
	// raw pointers for cup arrays
	real *velB_r  = thrust::raw_pointer_cast(&(velB[0])),
	     *uB     = thrust::raw_pointer_cast(&(B.uB[0])),
	     *vB     = thrust::raw_pointer_cast(&(B.vB[0]));

	cusp::blas::fill(velB, 0.0);
	
	const int blockSize = 256;
	dim3 dimGrid( int((B.totalPoints-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	
	kernels::fill_velB <<<dimGrid, dimBlock>>> (velB_r, uB, vB, B.totalPoints);
}
