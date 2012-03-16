#include <solvers/NavierStokes/NavierStokesSolver.h>

void NavierStokesSolver::assembleMatrices()
{
	generateA();
	generateQT();
	// QTBNQ
	generateC();
}

void NavierStokesSolver::initialise()
{
//	time_step = start_time_step;
	initialiseArrays();
	assembleMatrices();
}
void NavierStokesSolver::initialiseArrays()
{
}
void NavierStokesSolver::generateC()
{
}
void NavierStokesSolver::stepTime()
{
	for(int i=0; i<simPar->intSch.substeps; i++)
	{
		updateSolverState();
		
		generateRN();
		generateBC1();
		assembleRHS1();
	
		solveIntermediateVelocity();
	
		generateBC2();
		assembleRHS2();
	
		solvePoisson();
	
		projectionStep();
	}
	time_step++;
}
void NavierStokesSolver::assembleRHS1()
{
}
void NavierStokesSolver::solveIntermediateVelocity()
{
}
void NavierStokesSolver::assembleRHS2()
{
}
void NavierStokesSolver::solvePoisson()
{
}
void NavierStokesSolver::projectionStep()
{
}
void NavierStokesSolver::writeData()
{
	//if (time_step % nsave == 0)
	//{		
	//}
}
void NavierStokesSolver::updateBoundaryConditions()
{
}
void NavierStokesSolver::updateSolverState()
{
	updateBoundaryConditions();
}

bool NavierStokesSolver::finished()
{
	return (time_step < simPar->nt) ? false : true;
}