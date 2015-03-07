/***************************************************************************//**
 * \file generateQT.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to generate
 *        the divergence matrix and the interpolation matrix.
 */


#include <solvers/NavierStokes/kernels/generateQT.h>


/**
 * \brief To be documented.
 */
template <>
void DirectForcingSolver<host_memory>::updateQ()
{
}

/**
 * \brief After the Q matrix has been set up for the entire grid, this function
 *        removes the non-zeros that correspond to the velocity nodes where the 
 *        interpolation is performed, i.e. it sets the fluxes there to zero.
 *        This routine is specific to DirectForcingSolver.
 */
template <>
void DirectForcingSolver<device_memory>::updateQ()
{
	const int blocksize = 256;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  QSize = 4*nx*ny-2*(nx+ny);
	
	int  *QRows = thrust::raw_pointer_cast(&(Q.row_indices[0])),
	     *QCols = thrust::raw_pointer_cast(&(Q.column_indices[0]));
	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));

	real *QVals = thrust::raw_pointer_cast(&(Q.values[0]));
	
	dim3 dimGrid( int((QSize-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
	kernels::updateQ <<<dimGrid, dimBlock>>> (QRows, QCols, QVals, QSize, tags_r);
}

/**
 * \brief After the Q matrix has been set up for the entire grid, this function
 *        removes the non-zeros that correspond to the velocity nodes where the 
 *        interpolation is performed, i.e. it sets the fluxes there to zero.
 *        This routine is specific to DirectForcingSolver.
 */
template <typename memoryType>
void DirectForcingSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	updateQ();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}
