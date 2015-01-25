/***************************************************************************//**
 * \file generateA.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to generate
 *        the matrix resulting from implicit terms in the discretized momentum equation.
 */


#include <solvers/NavierStokes/kernels/generateA.h>


/**
 * \brief Generates the matrix A on the device.
 *
 * \param alpha Coefficient of the implicit part of the diffusion term.
 *              1 for backward Euler, 0.5 for Crank-Nicolson and 0 for explicit Euler.
 */
template <>
void DirectForcingSolver  <device_memory>::generateA(real alpha)
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	int  ASize = 5*( (nx-1)*ny + nx*(ny-1) ) - 4*(nx+ny) + 4,
	     numUV = (nx-1)*ny + nx*(ny-1);

	A.resize(numUV, numUV, ASize);

	int  *LRows = thrust::raw_pointer_cast(&(L.row_indices[0])),
	     *LCols = thrust::raw_pointer_cast(&(L.column_indices[0])),
	     *ARows = thrust::raw_pointer_cast(&(A.row_indices[0])),
	     *ACols = thrust::raw_pointer_cast(&(A.column_indices[0])),
	     *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));

	real *MVals = thrust::raw_pointer_cast(&(M.values[0])),
	     *LVals = thrust::raw_pointer_cast(&(L.values[0])),
	     *AVals = thrust::raw_pointer_cast(&(A.values[0]));

	const int blockSize = 256,
	          gridSize  = 1000;
	dim3 dimGrid(gridSize, 1);
	dim3 dimBlock(blockSize, 1);

	kernels::generateADirectForcing <<<dimGrid, dimBlock>>> (ARows, ACols, AVals, MVals, LRows, LCols, LVals, ASize, alpha, tags_r);
}

/**
 * \brief Generates the matrix A on the host.
 *
 * \param alpha Coefficient of the implicit part of the diffusion term.
 *              1 for backward Euler, 0.5 for Crank-Nicolson and 0 for explicit Euler.
 */
template <>
void DirectForcingSolver<host_memory>::generateA(real alpha)
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	int  ASize = 5*( (nx-1)*ny + nx*(ny-1) ) - 4*(nx+ny) + 4,
	     numUV = (nx-1)*ny + nx*(ny-1);

	A.resize(numUV, numUV, ASize);

	for(int i=0; i<ASize; i++)
	{
		A.row_indices[i] = L.row_indices[i];
		A.column_indices[i] = L.column_indices[i];
		if(tags[A.row_indices[i]]==-1)
			A.values[i] = -alpha*L.values[i];
		else
			A.values[i] = -L.values[i];
		if(A.row_indices[i] == A.column_indices[i])
			A.values[i] += M.values[A.row_indices[i]];
	}
}
