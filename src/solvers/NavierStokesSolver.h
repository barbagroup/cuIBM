/***************************************************************************//**
 * \file NavierStokesSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c NavierStokesSolver.
 */


#pragma once

#include "utilities/types.h"
#include "utilities/domain.h"
#include "utilities/integrationScheme.h"
#include "utilities/parameterDB.h"
#include "utilities/preconditioner.h"
#include "io/io.h"


/**
 * \class NavierStokesSolver
 * \brief Solves the Navier-Stokes equations in a rectangular domain.
 *
 * Methods are defined as virtual when they are redefined 
 * in a derived class with the same name.
 *
 */
template <typename memoryType>
class NavierStokesSolver
{
protected:
	parameterDB *paramDB;		///< database that contains all the simulation parameters
	domain      *domInfo;		///< computational grid information
	integrationScheme intgSchm;	///< instance of the class \c integrationScheme

	real QCoeff;
	
	cusp::coo_matrix<int, real, memoryType>
	    M,		///< diagonal mass matrix
	    Minv,	///< inverse of the mass matrix
	    L,		///< discrete Laplacian matrix
	    A,		///< combination of mass and Laplacian matrices
	    Q,		///< gradient matrix (and regularization matrix if immersed body)
	    QT,		///< transposed gradient matrix (and interpolation matrix if immersed body)
	    BN,		///< N-th order Taylor series expansion of the inverse of \c A
	    C;		///< matrix of the Poisson system

	cusp::array1d<real, memoryType>
	    q,			///< velocity flux vector
		qStar,		///< intermediate velocity flux vector
		qOld,		///< velocity flux vector at the previous time-step
		lambda,		///< pressure vector (and body forces if immersed body)
		rn,			///< explicit terms of the momentum equation
		bc1,		///< inhomogeneous boundary conditions from the discrete Laplacian operator
		bc2,		///< inhomogeneous boundary conditions from the discrete continuity equation
		rhs1,		///< right hand-side of the system for the intermediate velocity flux vector
		rhs2,		///< right hand-side of the Poisson system
		H,
		temp1,		///< temporary 1D Cusp array
		temp2,		///< temporary 1D Cusp array
		bc[4];		///< array that contains the boundary conditions of the rectangular domain
	     
	preconditioner< cusp::coo_matrix<int, real, memoryType> >
		*PC1,		///< preconditioner for the intermediate flux solver
		*PC2;		///< preconditioner for the Poisson solver

	size_t  timeStep,			///< iteration number
			subStep,			///< number of sub-iterations
			iterationCount1,	///< number of iteration to solve the intermediate fluxes
			iterationCount2;	///< number of iteration to solve the Poisson equation
	
	Logger logger;	///< instance of the class \c Logger to track time of different tasks
	
	std::ofstream iterationsFile;	///< file that contains the number of iterations
	
	// initialize parameters common to all Navier-Stokes solvers
	void initialiseCommon();
	
	// initialize all arrays required to solve the Navier-Stokes equations
	void initialiseArrays(int numQ, int numLambda);
	
	// initialize velocity flux vectors
	virtual void initialiseFluxes();
	
	// initialize velocity flux vectors
	virtual void initialiseFluxes(real *q);
	
	// initialize boundary velocity arrays with values stored in the database
	void initialiseBoundaryArrays();
	
	// assemble matrices of the intermediate flux solver and the Poisson solver
	void assembleMatrices();

	// generate the mass matrix and its inverse
	void generateM();

	// generate the discrete Laplacian matrix
	virtual void generateL();

	// generate the matrix resulting from the implicit flux terms
	virtual void generateA(real alpha);

	// generate approximate inverse of the matrix resulting from implicit velocity terms
	void generateBN();	

	// generate the discrete divergence matrix
	virtual void generateQT();

	// does nothing
	void updateQ(real gamma);
	
	// generate the matrix of the Poisson solver
	virtual void generateC();

	// generate explicit terms of the momemtum equation
	void generateRN();

	// generate explicit flux terms
	void calculateExplicitQTerms();

	// calculates the 
	virtual void calculateExplicitLambdaTerms();
	
	// generate inhomogeneous boundary conditions from the discrete Laplacian operator
	virtual void generateBC1();

	// generate inhomogeneous boundary conditions from the discrete continuity equation
	virtual void generateBC2();

	// assemble the right hand-side of the system for the intermediate flux
	virtual void assembleRHS1();

	// assemble the right hand-side of the Poisson system.
	void assembleRHS2();

	// solve for the intermediate flux velocity
	virtual void solveIntermediateVelocity();

	// solve the Poisson system for the pressure (and the body forces if immersed body)
	virtual void solvePoisson();

	// project the flux onto the divergence-free field
	virtual void projectionStep();

	// doing nothing
	void updateBoundaryConditions();

	// doing nothing
	virtual void updateSolverState();

public:
	// constructor -- copy the database and information about the computational grid
	NavierStokesSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);

	// initialize parameters, arrays and matrices required for the simulation
	virtual void initialise();
	
	// calculate the variables at the next time step
	void stepTime();
	
	// write numerical solution and number of iterations performed in each solver.
	virtual void writeCommon();
	
	// write data into files
	virtual void writeData();
	
	// evaluate the condition required to stop the simulation
	bool finished();
	
	// print timing information and clse the different files
	virtual void shutDown();
	
	/**
	 * \brief Returns the name of the current solver.
	 */
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
};
