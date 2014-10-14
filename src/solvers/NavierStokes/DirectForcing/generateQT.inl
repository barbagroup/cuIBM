#include <solvers/NavierStokes/kernels/generateQT.h>

template <>
void DirectForcingSolver<host_memory>::updateQ()
{
}

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

template <typename memoryType>
void DirectForcingSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	updateQ();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}
