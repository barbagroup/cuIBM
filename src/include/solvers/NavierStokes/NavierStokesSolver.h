#pragma once

#include <types.h>
#include <flowDescription.h>
#include <simulationParameters.h>
#include <domain.h>
#include <integrationScheme.h>

#include <cusp/print.h>
#include <cusp/elementwise.h>

/**
* Navier-Stokes solver for a rectangular domain.
*/
template <typename Matrix, typename Vector>
class NavierStokesSolver
{
protected:
	flowDescription *flowDesc;
	simulationParameters *simPar;
	domain  *domInfo;
	Vector  bcN[4], bcNP1[4];
	Matrix  M, Minv, L, A, QT, Q, BN, C;
	Vector  q, qStar, lambda, rn, H, rhs1, rhs2;
	int     timeStep;
	
	void initialiseArrays();
	void assembleMatrices(); // contains subfunctions to calculate A, QT, BN, QTBNQ
	
	void generateM();
	void generateL();
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
	static NavierStokesSolver<Matrix, Vector>* createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info);
};

/*template <>
NavierStokesSolver<coo_h, vec_h>* NavierStokesSolver<coo_h, vec_h>::createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
{
	NavierStokesSolver<coo_h, vec_h> *solver;
	switch(sim_par.ibmSch)
	{
		case NAVIER_STOKES:
			solver = new NavierStokesSolver<coo_h, vec_h>;
			break;
		case FADLUN_ET_AL:
			solver = new FadlunEtAlSolver<coo_h, vec_h>;
			break;
	}
	solver->flowDesc = &flow_desc;
	solver->simPar = &sim_par;
	solver->domainInfo = &dom_info;
	return solver;
}*/