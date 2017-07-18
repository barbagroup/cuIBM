/**
 * \file generateQT.inl
 * \brief Implementation of the methods to generate the discrete divergence matrix.
 */


#include <solvers/kernels/generateQT.h>


/**
 * \brief Generates the discrete divergence matrix.
 *
 * The gradient matrix is also computed as the transpose of the divergence one.
 */
template <>
void NavierStokesSolver<device_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int nx = domInfo->nx,
	    ny = domInfo->ny;
	
	int numU  = (nx-1)*ny;
	int numUV = numU + nx*(ny-1);
	int numP  = numU + ny;
	
	cooH QTHost(numP, numUV, 4*numP-2*(nx+ny));
	
	int *QTRows = thrust::raw_pointer_cast(&(QTHost.row_indices[0])),
	    *QTCols = thrust::raw_pointer_cast(&(QTHost.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QTHost.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);
	QT = QTHost;
	cusp::transpose(QT, Q);
	
	logger.stopTimer("generateQT");
} // generateQT
