#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>

void NavierStokesSolver::assembleMatrices()
{
	generateA();
	generateQT();
	// QTBNQ
	generateC();
}
void NavierStokesSolver::initialise()
{
	timeStep = simPar->startStep;
	initialiseArrays();
	assembleMatrices();
}
void NavierStokesSolver::initialiseArrays()
{
}
void NavierStokesSolver::generateA()
{
}
void NavierStokesSolver::generateQT()
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
	timeStep++;
}
void NavierStokesSolver::generateRN()
{
}
void NavierStokesSolver::generateBC1()
{
}
void NavierStokesSolver::assembleRHS1()
{
}
void NavierStokesSolver::solveIntermediateVelocity()
{
}
void NavierStokesSolver::generateBC2()
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
	if (timeStep % simPar->nsave == 0)
	{		
	}
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
	return (timeStep < simPar->nt) ? false : true;
}

/**
* \brief Factory method to select the required IBM solver
* \param a Description
* \return
*/
NavierStokesSolver* NavierStokesSolver::createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
{
	NavierStokesSolver *solver;
	switch(sim_par.ibmSch)
	{
		case NAVIER_STOKES:
			solver = new NavierStokesSolver;
			break;
		case FADLUN_ET_AL:
			solver = new FadlunEtAlSolver;
			break;
	}
	solver->flowDesc = &flow_desc;
	solver->simPar = &sim_par;
	solver->domainInfo = &dom_info;
	return solver;
}