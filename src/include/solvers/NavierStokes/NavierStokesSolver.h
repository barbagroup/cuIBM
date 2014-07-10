/***************************************************************************//**
* \file NavierStokesSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Solves the Navier-Stokes equations in a rectangular domain
*/

#pragma once

#include <types.h>
#include <domain.h>
#include <integrationScheme.h>
#include <io/io.h>
#include <parameterDB.h>
#include <preconditioner.h>

/********************//**
* \class NavierStrokesSolver 
* \brief Navier-Stokes solver for a rectangular domain
*/
template <typename memoryType>
class NavierStokesSolver
{
protected:
	parameterDB *paramDB;		///< pointer to the database containing all simulation parameters
	domain      *domInfo;		///< pointer to the computational grid information
	integrationScheme intgSchm;	///< object of the class \c integrationScheme

	real QCoeff;
	
	cusp::coo_matrix<int, real, memoryType>
	     M,		///< diagonal mass matrix
	     Minv,	///< inverse of the mass matrix
	     L,		///< Discrete Laplacian 
	     A,		///< mass and Laplacian matrices
	     QT,	///< transposed matrix of Q
	     Q,		///< gathers the gradient matrix and the transpose of the extrapolation matrix
	     BN,	///< N-th order Taylor series expansion of  A^-1
	     C;		///< matrix product \f$ Q^T B^N Q \f$

	cusp::array1d<real, memoryType>
	     q, qStar, lambda, rn, H, rhs1, rhs2, bc1, bc2, temp2, temp1, bc[4], qOld;
	     
	preconditioner< cusp::coo_matrix<int, real, memoryType> > *PC1, *PC2;

	size_t  timeStep, subStep, iterationCount1, iterationCount2;
	
	Logger logger;
	
	std::ofstream iterationsFile;
	
	/**
	* \brief Initialize stuff common to all Immersed Boundary Method solvers
	*/
	void initialiseCommon();
	
	/**
	* \brief Initialize all required arrays
	*/
	void initialiseArrays(int numQ, int numLambda);
	
	/**
	* \brief Set the initial value of all the fluxes in the flow field
	*/
	virtual void initialiseFluxes();
	
	/**
	* \brief
	*/
	virtual void initialiseFluxes(real *q);
	
	/**
	* \brief Set the initial values of the boundary velocities
	*/
	void initialiseBoundaryArrays();
	
	/**
	* \brief Assemble all the required matrices
	*/
	void assembleMatrices(); // contains subfunctions to calculate A, QT, BN, QTBNQ

	// Methods are defined as virtual when they are redefined in a derived class with the same name.
	
	// functions to generate matrices
	/**
	* \brief Generate the mass matrix
	*/
	void generateM();

	/**
	* \brief Generate the Laplacian matrix
	*/
	virtual void generateL();

	/**
	* \brief Generate the matrix \f$ A = M + L \f$
	*/
	virtual void generateA(real alpha);

	/**
	* \brief Computes the Nth-order Taylor expansion of the matrix A
	*/
	void generateBN();	

	/**
	* \brief Generate the matrix containing the discrete divergence operator
	*/
	virtual void generateQT();

	/**
	* \brief
	*/
	void updateQ(real gamma);
	
	/**
	* \brief Compute the matrix product \f$ C = Q^T B^N Q \f$
	*/
	void generateC();

	// generate explicit terms
	void generateRN();
	void calculateExplicitQTerms();
	void calculateExplicitLambdaTerms();
	
	virtual void generateBC1();
	virtual void generateBC2();

	virtual void assembleRHS1();
	void assembleRHS2();

	virtual void solveIntermediateVelocity();
	void solvePoisson();
	virtual void projectionStep();

	void updateBoundaryConditions();
	virtual void updateSolverState();

public:
	/**
	* \brief Constructor of the class \c NavierStokesSolver
	*/
	NavierStokesSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);

	/**
	* \brief Initialize stuff required for the simulation
	*/
	virtual void initialise();
	
	/**
	* \brief Calculate all the variables at the next time step
	*/
	void stepTime();
	
	/**
	* \brief Write common data to files
	*/
	virtual void writeCommon();
	
	/**
	* \brief Write data to files
	*/
	virtual void writeData();
	
	/**
	* \brief Condition required to bring the simulation to a halt.
	*/
	bool finished();
	
	/**
	* \brief Perform necessary actions to end the simulation
	*/
	virtual void shutDown();
	
	// Factory methods are static (not entirely sure why)
	/**
	* \brief  Factory method to select the required IBM solver
	*/
	//	static NavierStokesSolver<memoryType>* createSolver(parameterDB &paramDB, domain &domInfo);
	
	/**
	* \brief Gives the name of the current solver 
	*/
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
};
