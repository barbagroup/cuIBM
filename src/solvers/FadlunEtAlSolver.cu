/**
 * \file FadlunEtAlSolver.cu
 * \brief Implementation of the methods of the classes \c FadlunEtAlSolver and \c FEAModifiedSolver.
 */


#include <sys/stat.h>

#include "FadlunEtAlSolver.h"
#include "solvers/kernels/generateQT.h"


/**
 * \brief Constructor. Copies the database and information about the 
 *        computational grid.
 */
template <typename memoryType>
FadlunEtAlSolver<memoryType>::FadlunEtAlSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
} // FadlunEtAlSolver


/**
 * \brief Generates the matrix \c QT, Q and G for FadlunEtAlSolver.
 */
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	G = NavierStokesSolver<memoryType>::Q;
	updateG();
} // generateQT


/**
 * \brief Update the gradient operator.
 *
 * Zeros rows of G that correspond to the forcing nodes.
 * Calls the function updateQ, but passes the matrix G to it.
 */
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
} // updateG


/**
 * \brief Add the pressure gradient to the right hand side of the 
 *        momentum equation.
 */
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::calculateExplicitLambdaTerms()
{
	// temp_1 = G.pressure
	cusp::multiply(G, DirectForcingSolver<memoryType>::pressure, NavierStokesSolver<memoryType>::temp1);
	// r^n = r^n - temp_1
	cusp::blas::axpy(NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rn, -1.0);
} // calculateExplicitLambdaTerms


/**
 * \brief Generates the matrix for the Poisson equation.
 *
 * Calls the function from NavierStokesSolver because it does not need to set
 * any node inside the immersed boundary to zero.
 */
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateC()
{
	NavierStokesSolver<memoryType>::generateC();
} // generateC


// specialization of the class
template class FadlunEtAlSolver<device_memory>;


/**
 * \brief Constructor. Copies the database and information about the 
 *        computational grid.
 */
template <typename memoryType>
FEAModifiedSolver<memoryType>::FEAModifiedSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
} // FEAModifiedSolver


/**
 * \brief Generates the matrices QT and Q, which are the same as if an 
 *        immersed boundary is not present in the fluid.
 */
template <typename memoryType>
void FEAModifiedSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
} // generateQT


/**
 * \brief Generates the matrix for the Poisson equation.
 *
 * Calls the function from NavierStokesSolver because it does not need to set
 * any node inside the immersed boundary to zero.
 */
template <typename memoryType>
void FEAModifiedSolver<memoryType>::generateC()
{
	NavierStokesSolver<memoryType>::generateC();
} // generateC


// specialization of the class
template class FEAModifiedSolver<device_memory>;
