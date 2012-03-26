#pragma once

#include <types.h>
#include <flowDescription.h>
#include <simulationParameters.h>
#include <domain.h>
#include <integrationScheme.h>
#include <io.h>

/**
* Navier-Stokes solver for a rectangular domain.
*/
template <typename Matrix, typename Vector>
class NavierStokesSolver
{
protected:
	options *opts;
	flowDescription *flowDesc;
	simulationParameters *simPar;
	domain  *domInfo;
	Vector  bc[4];//, bcNP1[4];
	vecD    bcHost[4];
	Matrix  M, Minv, L, A, QT, Q, BN, C;
	Vector  q, qStar, lambda, rn, H, rhs1, rhs2, bc1, bc2, temp2, temp1;
	vecD    bc2Host;
	int     timeStep;
	
	/**
	* Methods are defined as virtual when they are redefined with the same name in a derived class.
	*/	
	virtual void initialiseArrays();
	void initialiseFluxes();
	void initialiseBoundaryArrays();
	void assembleMatrices(); // contains subfunctions to calculate A, QT, BN, QTBNQ
	
	void generateM();
	void generateL();
	virtual void generateA();
	void generateBN();
	virtual void generateQT();
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
	virtual void initialise();
	void stepTime();
	void writeData();
	bool finished();
	static NavierStokesSolver<Matrix, Vector>* createSolver(options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info);
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
};