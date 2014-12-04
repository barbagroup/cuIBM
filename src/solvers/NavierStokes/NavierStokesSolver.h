/**
* @file NavierStokesSolver.h
* @brief Solves the Navier-Stokes equations in a rectangular domain.
*/

#pragma once

#include <types.h>
#include <domain.h>
#include <integrationScheme.h>
#include <io/io.h>
#include <parameterDB.h>
#include <preconditioner.h>

/**
* Navier-Stokes solver for a rectangular domain.
*/
template <typename memoryType>
class NavierStokesSolver
{
protected:
	parameterDB *paramDB;
	domain      *domInfo;
	integrationScheme intgSchm;

	real QCoeff;
	
	cusp::coo_matrix<int, real, memoryType>
	     M,
	     Minv, 
	     L,
	     A, 
	     QT,
	     Q, 
	     BN, 
	     C;

	cusp::array1d<real, memoryType>
	     q, qStar, lambda, rn, H, rhs1, rhs2, bc1, bc2, temp2, temp1, bc[4], qOld;
	     
	preconditioner< cusp::coo_matrix<int, real, memoryType> > *PC1, *PC2;

	size_t  timeStep, subStep, iterationCount1, iterationCount2;
	
	Logger logger;
	
	std::ofstream iterationsFile;
	
	/**
	* Initialises stuff common to all IBM solvers
	*/
	void initialiseCommon();
	
	/**
	* @brief Initialises all required arrays
	* @param numQ Total number velocity variables (u and v)
	* @param numLambda Number of pressure variables + twice the number of body force variables
	*/
	void initialiseArrays(int numQ, int numLambda);
	virtual void initialiseFluxes();
	virtual void initialiseFluxes(real *q);
	void initialiseBoundaryArrays();
	void assembleMatrices(); // contains subfunctions to calculate A, QT, BN, QTBNQ

	// Methods are defined as virtual when they are redefined in a derived class with the same name.
	
	// functions to generate matrices
	void generateM();
	virtual void generateL();
	virtual void generateA(real alpha);
	void generateBN();	
	virtual void generateQT();
	void updateQ(real gamma);
	
	virtual void generateC();

	// generate explicit terms
	void generateRN();
	void calculateExplicitQTerms();
	virtual void calculateExplicitLambdaTerms();
	
	virtual void generateBC1();
	virtual void generateBC2();

	virtual void assembleRHS1();
	void assembleRHS2();

	virtual void solveIntermediateVelocity();
	virtual void solvePoisson();
	virtual void projectionStep();

	void updateBoundaryConditions();
	virtual void updateSolverState();

public:
	NavierStokesSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	/**
	* @brief Initialise stuff required for the simulation
	*/
	virtual void initialise();
	
	/**
	* @brief Calculate all the variables at the next time step
	*/
	void stepTime();
	
	/**
	* @brief Write the data to files
	*/
	virtual void writeCommon();
	virtual void writeData();
	
	/**
	* @brief Condition required to bring the simulation to a halt.
	* @return True if the simulation is over. False if it must continue.
	*/
	bool finished();
	
	/**
	* @brief Perform necessary actions to end the simulation
	*/
	virtual void shutDown();
	
	/**
	* @brief Give the name of the current solver 
	* @return String that describes the type of solver
	*/
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
};
