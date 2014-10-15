#include "FadlunEtAlSolver.h"
#include <sys/stat.h>

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}

template <typename memoryType>
FadlunEtAlSolver<memoryType>::FadlunEtAlSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
