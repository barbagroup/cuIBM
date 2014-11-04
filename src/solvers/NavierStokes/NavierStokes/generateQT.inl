/***************************************************************************//**
 * \file generateQT.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Generates the matrix containing the discrete divergence operator.
 */


#include <solvers/NavierStokes/kernels/generateQT.h>


/**
 * \brief Generates the matrix containing the discrete divergence operator on the host.
 */
template <>
void NavierStokesSolver<host_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	QT.resize(numP, numUV, 4*numP-2*(nx+ny));
	
	int  *QTRows = thrust::raw_pointer_cast(&(QT.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QT.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QT.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);
	cusp::transpose(QT, Q);
	
	logger.stopTimer("generateQT");
}

/**
 * \brief Generates the matrix containing the discrete divergence operator.
 */
template <>
void NavierStokesSolver<device_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	cooH QTHost(numP, numUV, 4*numP-2*(nx+ny));
	
	int  *QTRows = thrust::raw_pointer_cast(&(QTHost.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QTHost.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QTHost.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);
	QT = QTHost;
	cusp::transpose(QT, Q);
	
	logger.stopTimer("generateQT");
}
