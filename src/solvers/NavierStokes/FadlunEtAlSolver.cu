#include <solvers/NavierStokes/FadlunEtAlSolver.h>

template <typename Matrix, typename Vector>
void FadlunEtAlSolver<Matrix, Vector>::generateA()
{
}
template <typename Matrix, typename Vector>
void FadlunEtAlSolver<Matrix, Vector>::updateA()
{
}
template <typename Matrix, typename Vector>
void FadlunEtAlSolver<Matrix, Vector>::initialiseBodies()
{
}
template <typename Matrix, typename Vector>
void FadlunEtAlSolver<Matrix, Vector>::updateBodies()
{
}
template <typename Matrix, typename Vector>
void FadlunEtAlSolver<Matrix, Vector>::initialise()
{
	NavierStokesSolver<Matrix, Vector>::initialiseArrays();
	NavierStokesSolver<Matrix, Vector>::assembleMatrices();
	initialiseBodies();
}
template <typename Matrix, typename Vector>
void FadlunEtAlSolver<Matrix, Vector>::updateSolverState()
{
	NavierStokesSolver<Matrix, Vector>::updateBoundaryConditions();
	if (B.bodiesMove) {
		updateBodies();
		updateA();
	}
}

template class FadlunEtAlSolver<cooH, vecH>;
template class FadlunEtAlSolver<cooD, vecD>;