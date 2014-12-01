#include "FadlunEtAlSolver.h"
#include <sys/stat.h>

template <typename memoryType>
FadlunEtAlSolver<memoryType>::FadlunEtAlSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}

// rn = rn - Q.pressure
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::calculateExplicitLambdaTerms()
{
	// temp1 = Q.pressure
	cusp::multiply(NavierStokesSolver<memoryType>::Q, DirectForcingSolver<memoryType>::pressure, NavierStokesSolver<memoryType>::temp1);
	// rn = rn - temp1
	cusp::blas::axpy(NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rn, -1.0);
}

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
