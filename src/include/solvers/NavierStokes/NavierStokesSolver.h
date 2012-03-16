#pragma once

#include <types.h>
#include <flowDescription.h>
#include <simulationParameters.h>
#include <domain.h>
#include <integrationScheme.h>

/**
* Navier-Stokes solver for a rectangular domain.
*/
class NavierStokesSolver
{
protected:
	flowDescription *flowDesc;
	simulationParameters *simPar;
	domain  *domainInfo;
	vector  bc[4];
	matrix  A, QT, Q, BN, QTBNQ;
	vector  q, q_star, lambda, rn, H;
	int     timeStep;
	
	void initialiseArrays();
	void assembleMatrices(); // contains subfunctions to calculate A, QT, BN, QTBNQ
	
	void generateA();
	void generateQT();
	void generateC();

	void generateRN();
	void generateBC1();
	void generateBC2();
	
	void assembleRHS1();
	void assembleRHS2();
	
	void solveIntermediateVelocity();
	void solvePoisson();
	void projectionStep();
	
	void updateBoundaryConditions();
	void updateSolverState();
	
public:
	void initialise();
	void stepTime();
	void writeData();
	bool finished();
	//friend NavierStokesSolver* allocator::createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info);
	static NavierStokesSolver* createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info);
};