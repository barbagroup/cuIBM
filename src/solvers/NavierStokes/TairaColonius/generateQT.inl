/***************************************************************************//**
 * \file generateQT.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c TairaColoniusSolver
 *        to generate the transposed gradient matrix and the interpolation matrix.
 */


#include <solvers/NavierStokes/kernels/generateQT.h>

#define BLOCKSIZE 256


/**
 * \brief Updates the interpolation matrix using the current locations of body points (on the host).
 *        Typically called after the body has moved.
 */
template <>
void TairaColoniusSolver<host_memory>::updateQT()
{
	logger.startTimer("updateQT");
	
	int  *QTRows = thrust::raw_pointer_cast(&(QT.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QT.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QT.values[0]));
	
	int  *ERows = thrust::raw_pointer_cast(&(E.row_indices[0])),
	     *ECols = thrust::raw_pointer_cast(&(E.column_indices[0]));
	real *EVals = thrust::raw_pointer_cast(&(E.values[0]));
	
	real *x  = &(domInfo->x[0]),
	     *y  = &(domInfo->y[0]),
	     *dx = &(domInfo->dx[0]);
	
	real *xB = thrust::raw_pointer_cast(&(B.x[0])),
	     *yB = thrust::raw_pointer_cast(&(B.y[0]));
	
	int *I = thrust::raw_pointer_cast(&(B.I[0])),
	    *J = thrust::raw_pointer_cast(&(B.J[0]));
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	     
	kernels::updateQTHost(QTRows, QTCols, QTVals,
	                      ERows,  ECols,  EVals,
	                      nx, ny, x, y, dx,
                          B.totalPoints, xB, yB, I, J);

	logger.stopTimer("updateQT");
	
	logger.startTimer("transposeQT");
	cusp::transpose(QT, Q);
	logger.stopTimer("transposeQT");
	
	logger.startTimer("transposeE");
	cusp::transpose(E, ET);
	logger.stopTimer("transposeE");
}

/**
 * \brief Updates the interpolation matrix using the current locations of body points (on the device).
 *        Typically called after the body has moved.
 */
template <>
void TairaColoniusSolver<device_memory>::updateQT()
{
	logger.startTimer("updateQT");
	
	int  *QTRows = thrust::raw_pointer_cast(&(QT.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QT.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QT.values[0]));
	
	int  *ERows = thrust::raw_pointer_cast(&(E.row_indices[0])),
	     *ECols = thrust::raw_pointer_cast(&(E.column_indices[0]));
	real *EVals = thrust::raw_pointer_cast(&(E.values[0]));
	
	real *x  = thrust::raw_pointer_cast(&(domInfo->xD[0])),
	     *y  = thrust::raw_pointer_cast(&(domInfo->yD[0])),
	     *dx = thrust::raw_pointer_cast(&(domInfo->dxD[0]));
	
	real *xB = thrust::raw_pointer_cast(&(B.x[0])),
	     *yB = thrust::raw_pointer_cast(&(B.y[0]));
	
	int *I = thrust::raw_pointer_cast(&(B.I[0])),
	    *J = thrust::raw_pointer_cast(&(B.J[0]));
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	     
	dim3 dimGrid(int((B.totalPoints-0.5)/BLOCKSIZE)+1, 1);
	dim3 dimBlock(BLOCKSIZE, 1);

	kernels::updateQT <<<dimGrid, dimBlock>>>
	         (QTRows, QTCols, QTVals,
	          ERows,  ECols,  EVals,
	          nx, ny, x, y, dx,
              B.totalPoints, xB, yB, I, J);
    
    logger.stopTimer("updateQT");
    
	logger.startTimer("transposeQT");
	cusp::transpose(QT, Q);
	logger.stopTimer("transposeQT");
	
	logger.startTimer("transposeE");
	cusp::transpose(E, ET);
	logger.stopTimer("transposeE");
}

/**
 * \brief Generates the transposed gradient matrix and interpolation matrix (on the host).
 */
template <>
void TairaColoniusSolver<host_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	/// QT is an (np + 2*nb) x nuv matrix
	QT.resize(numP + 2*B.totalPoints, numUV, 4*numP-2*(nx+ny) + 24*B.totalPoints);
	
	int  *QTRows = thrust::raw_pointer_cast(&(QT.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QT.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QT.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);

	logger.stopTimer("generateQT");
	
	updateQT();	
}

/**
 * \brief Generates the transposed gradient matrix and interpolation matrix (on the device).
 */
template <>
void TairaColoniusSolver<device_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	/// QT is an (np + 2*nb) x nuv matrix
	cooH QTHost(numP + 2*B.totalPoints, numUV, 4*numP-2*(nx+ny) + 24*B.totalPoints);
	
	int  *QTRows = thrust::raw_pointer_cast(&(QTHost.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QTHost.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QTHost.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);

	QT = QTHost;
	
	logger.stopTimer("generateQT");
	
	updateQT();
}
