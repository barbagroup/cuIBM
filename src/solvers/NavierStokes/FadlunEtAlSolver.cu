/***************************************************************************//**
 * \file FadlunEtAlSolver.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c FadlunEtAlSolver.
 */


#include "FadlunEtAlSolver.h"
#include <sys/stat.h>
#include <solvers/NavierStokes/kernels/generateQT.h>

template <typename memoryType>
FadlunEtAlSolver<memoryType>::FadlunEtAlSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}


/**
 * \brief Generates the matrix \c QT for FadlunEtAlSolver.
 */
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	G = NavierStokesSolver<memoryType>::Q;
	updateG();
}

template <>
void FadlunEtAlSolver<host_memory>::updateG()
{
}

template <>
void FadlunEtAlSolver<device_memory>::updateG()
{
	const int blocksize = 256;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  GSize = 4*nx*ny-2*(nx+ny);
	
	int  *GRows = thrust::raw_pointer_cast(&(G.row_indices[0])),
	     *GCols = thrust::raw_pointer_cast(&(G.column_indices[0]));
	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));

	real *GVals = thrust::raw_pointer_cast(&(G.values[0]));
	
	dim3 dimGrid( int((GSize-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
	kernels::updateQ <<<dimGrid, dimBlock>>> (GRows, GCols, GVals, GSize, tags_r);
}

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 */
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::calculateExplicitLambdaTerms()
{
	// temp1 = G.pressure
	cusp::multiply(G, DirectForcingSolver<memoryType>::pressure, NavierStokesSolver<memoryType>::temp1);
	// rn = rn - temp1
	cusp::blas::axpy(NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rn, -1.0);
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateC()
{
	NavierStokesSolver<memoryType>::generateC();
}

// specialization of the class
template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;

//----------------------------------------------------------------------------//

template <typename memoryType>
FEAModifiedSolver<memoryType>::FEAModifiedSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template <typename memoryType>
void FEAModifiedSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
}

template <typename memoryType>
void FEAModifiedSolver<memoryType>::generateC()
{
	NavierStokesSolver<memoryType>::generateC();
}

template class FEAModifiedSolver<host_memory>;
template class FEAModifiedSolver<device_memory>;
