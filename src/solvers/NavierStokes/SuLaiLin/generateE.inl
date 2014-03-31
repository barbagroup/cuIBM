#include <solvers/NavierStokes/kernels/generateE.h>

template <>
void SuLaiLinSolver<host_memory>::generateE()
{
	logger.startTimer("generateE");
	
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
	
	kernels::generateEHost(ERows, ECols, EVals,
	                       nx, ny, x, y, dx,
	                       B.totalPoints, xB, yB, I, J);
	
	logger.stopTimer("generateE");
	
	logger.startTimer("transposeE");
	cusp::transpose(E, ET);
	logger.stopTimer("transposeE");
}

template <>
void SuLaiLinSolver<device_memory>::generateE()
{
	logger.startTimer("generateE");
	
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
	     
	const int blockSize = 256;
	     
	dim3 dimGrid(int((B.totalPoints-0.5)/blockSize)+1, 1);
	dim3 dimBlock(blockSize, 1);

	kernels::generateE <<<dimGrid, dimBlock>>> (ERows, ECols, EVals,
	                                            nx, ny, x, y, dx, 
	                                            B.totalPoints, xB, yB, I, J);
    
    logger.stopTimer("generateE");
	
	logger.startTimer("transposeE");
	cusp::transpose(E, ET);
	logger.stopTimer("transposeE");
}
