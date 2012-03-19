#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::initialise()
{
	timeStep = simPar->startStep;
	initialiseArrays();
	assembleMatrices();
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::assembleMatrices()
{
	generateA();
	generateBN();
	generateQT();
	generateC(); // QT*BN*Q
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::initialiseArrays()
{
	int nx = domInfo->nx,
		ny = domInfo->ny;
		
	int numU  = (nx-1)*ny;
	int numUV = numU + nx*(ny-1);
	int numP  = numU + ny;
//	int numB  = 0;
//	for(int k=0; k<flowDesc->numBodies; k++)
//		numB += flowDesc->B[k].numPoints;
	
	    q.resize(numUV);
	qStar.resize(numUV);
	   rn.resize(numUV);
	    H.resize(numUV);
	 rhs1.resize(numUV);
	
	//lambda.resize(numP+2*numB);
	//rhs2.resize(numP+2*numB);
	lambda.resize(numP);
	  rhs2.resize(numP);
	
	/// Initialise velocity fluxes
	int i;
	for(i=0; i < numU; i++)
	{
		q[i] = flowDesc->initialU * domInfo->dy[i/(nx-1)];
	}
	for(; i < numUV; i++)
	{
		q[i] = flowDesc->initialV * domInfo->dx[(i-numU)%nx];
	}
	qStar = q;
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateC()
{
	Matrix temp;
	cusp::multiply(QT, BN, temp);
	cusp::multiply(temp, Q, C);
	C.sort_by_row_and_column();
	C.values[0] += C.values[0];
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::stepTime()
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

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateBC1()
{
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::assembleRHS1()
{
	
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::solveIntermediateVelocity()
{/*
	cusp::default_monitor<real> sys1Mon(r1, 10000);
	cusp::krylov::cg(A, qStar, r1, sys1Mon, PC1);
*/}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateBC2()
{
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::assembleRHS2()
{
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::solvePoisson()
{/*
	cusp::default_monitor<real> sys2Mon(r2, 20000);
	cusp::krylov::cg(C, lambda, r2, sys2Mon, PC2);
*/}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::projectionStep()
{/*
	// QUESTION: Would you create a temp1 here at every time step or have a permanent temp1?
	cusp::multiply(Q, lambda, temp1);
	cusp::multiply(BN, temp1, q);
	cusp::blas::axpby(qStar, q, q, 1.0, -1.0 );
*/}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::writeData()
{
	if (timeStep % simPar->nsave == 0)
	{		
	}
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::updateBoundaryConditions()
{
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::updateSolverState()
{
	updateBoundaryConditions();
}

template <typename Matrix, typename Vector>
bool NavierStokesSolver<Matrix, Vector>::finished()
{
	return (timeStep < simPar->nt) ? false : true;
}

/**
* \brief Factory method to select the required IBM solver
* \param a Description
* \return
*/

template <typename Matrix, typename Vector>
NavierStokesSolver<Matrix, Vector>* NavierStokesSolver<Matrix, Vector>::createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
{
	NavierStokesSolver<Matrix, Vector> *solver = 0;
	switch(sim_par.ibmSch)
	{
		case SAIKI_BIRINGEN:
			break;
		case TAIRA_COLONIUS:
			break;
		case NAVIER_STOKES:
			solver = new NavierStokesSolver<Matrix, Vector>;
			break;
		case FADLUN_ET_AL:
			solver = new FadlunEtAlSolver<Matrix, Vector>;
			break;
	}
	solver->flowDesc = &flow_desc;
	solver->simPar = &sim_par;
	solver->domInfo = &dom_info;
	return solver;
}

#include "NavierStokes/generateA.inl"
#include "NavierStokes/generateQT.inl"
#include "NavierStokes/generateRN.inl"

template class NavierStokesSolver<coo_h, vec_h>;
template class NavierStokesSolver<coo_d, vec_d>;