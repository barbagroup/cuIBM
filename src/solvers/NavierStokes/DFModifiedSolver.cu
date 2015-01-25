#include "DFModifiedSolver.h"
#include <sys/stat.h>

template <typename memoryType>
DFModifiedSolver<memoryType>::DFModifiedSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

// rn = rn - Q.pressure
template <typename memoryType>
void DFModifiedSolver<memoryType>::calculateExplicitLambdaTerms()
{
	// temp1 = Q.pressure
	cusp::multiply(NavierStokesSolver<memoryType>::Q, DirectForcingSolver<memoryType>::pressure, NavierStokesSolver<memoryType>::temp1);
	// rn = rn - temp1
	cusp::blas::axpy(NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rn, -1.0);
}

template <typename memoryType>
void DFModifiedSolver<memoryType>::projectionStep()
{
	NavierStokesSolver<memoryType>::projectionStep();
	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");
	cusp::blas::axpy(NavierStokesSolver<memoryType>::lambda, DirectForcingSolver<memoryType>::pressure, 1.0);
	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
}

template class DFModifiedSolver<host_memory>;
template class DFModifiedSolver<device_memory>;