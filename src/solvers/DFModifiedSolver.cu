/**
 * \file DFModifiedSolver.cu
 * \brief Implementation of the methods of the class \c DFModifiedSolver.
 */


#include <sys/stat.h>

#include "DFModifiedSolver.h"


/**
 * \brief Constructor -- get simulation parameters and grid.
 */
template <typename memoryType>
DFModifiedSolver<memoryType>::DFModifiedSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
} // DFModifiedSolver


/**
 * \brief Compute and add the explicit pressure gradient to the RHS of velocity system.
 */
template <typename memoryType>
void DFModifiedSolver<memoryType>::calculateExplicitLambdaTerms()
{
	// temp_1 = Q.pressure
	cusp::multiply(NavierStokesSolver<memoryType>::Q, DirectForcingSolver<memoryType>::pressure, NavierStokesSolver<memoryType>::temp1);
	// r^n = r^n - temp1
	cusp::blas::axpy(NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rn, -1.0);
} // calculateExplicitLambdaTerms


/**
 * \brief Project velocity field onto divergence-free space and update pressure.
 */
template <typename memoryType>
void DFModifiedSolver<memoryType>::projectionStep()
{
	// u^{n+1} = u^{*} - B^N.Q.\phi 
	NavierStokesSolver<memoryType>::projectionStep();

	// p^{n+1} = p^{n} + \phi
	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");
	cusp::blas::axpy(NavierStokesSolver<memoryType>::lambda, DirectForcingSolver<memoryType>::pressure, 1.0);
	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
} // projectionStep


template class DFModifiedSolver<host_memory>;
template class DFModifiedSolver<device_memory>;
