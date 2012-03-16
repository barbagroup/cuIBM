#include <solvers/NavierStokes/FadlunEtAlSolver.h>

void FadlunEtAlSolver::generateA()
{
}
void FadlunEtAlSolver::updateA()
{
}
void FadlunEtAlSolver::initialiseBodies()
{
}
void FadlunEtAlSolver::updateBodies()
{
}
void FadlunEtAlSolver::initialise()
{
	initialiseArrays();
	assembleMatrices();
	initialiseBodies();
}
void FadlunEtAlSolver::updateSolverState()
{
	updateBoundaryConditions();
	if (B.bodiesMove) {
		updateBodies();
		updateA();
	}
}