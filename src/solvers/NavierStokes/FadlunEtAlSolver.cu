#include <solvers/NavierStokes/FadlunEtAlSolver.h>

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateA()
{
}
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateA()
{
}
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialiseBodies()
{
}
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateBodies()
{
}
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialise()
{
	NavierStokesSolver<memoryType>::initialiseArrays();
	NavierStokesSolver<memoryType>::assembleMatrices();
	initialiseBodies();
}
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateSolverState()
{
	NavierStokesSolver<memoryType>::updateBoundaryConditions();
	if (B.bodiesMove) {
		updateBodies();
		updateA();
	}
}

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;