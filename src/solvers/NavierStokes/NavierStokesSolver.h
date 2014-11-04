/***************************************************************************//**
* \file NavierStokesSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the class \c NavierStokesSolver
*/

#pragma once

#include <types.h>
#include <domain.h>
#include <integrationScheme.h>
#include <io/io.h>
#include <parameterDB.h>
#include <preconditioner.h>

/**
* \class NavierStokesSolver
* \brief Solve the Navier-Stokes equations in a rectangular domain
*/
template <typename memoryType>
class NavierStokesSolver
{
protected:
	parameterDB *paramDB;		///< pointer to the database that contains all the simulation parameters
	domain      *domInfo;		///< pointer to the computational grid information
	integrationScheme intgSchm;	///< instance of the class \c integrationScheme

	real QCoeff;
	
	cusp::coo_matrix<int, real, memoryType>
	    M,		///< diagonal mass matrix
	    Minv,	///< inverse of the mass matrix
	    L,		///< Discrete Laplacian matrix
	    A,		///< mass and Laplacian matrices
	    Q,		///< gradient matrix plus transposed extrapolation matrix (if immersed body)
	    QT,		///< transposed matrix of Q
	    BN,		///< N-th order Taylor series expansion of  \f$ A^{-1} \f$
	    C;		///< matrix product \f$ Q^T B^N Q \f$

	cusp::array1d<real, memoryType>
	    q,			///< fluxes
		qStar,		///< intermediate fluxes
		qOld,		///< fluxes at the previous time-step
		lambda,		///< pressures and body forces (if immersed body)
		rn,			///< explicit terms in the momentum equation
		bc1,		///< inhomogeneous boundary condition from the discrete Laplacian operator
		bc2,		///< inhomogeneous boundary condition from the discrete continuity equation
		rhs1,		///< \c rn + \c bc1
		rhs2,		///< \c QT * \c qStar - \c bc2
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
	
	Logger logger;	///< object of the class \c Logger to track time of different tasks
	
	std::ofstream iterationsFile;	///< file that contains the number of iterations
	
	// Methods are defined as virtual when they are redefined in a derived class with the same name.
	
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
	* \brief Set the initial value of all the fluxes in the flow field
	*/
	virtual void initialiseFluxes(real *q);
	
	/**
	* \brief Set the initial values of the boundary velocities
	*/
	void initialiseBoundaryArrays();
	
	/**
	* \brief Assemble all the required matrices
	*/
	void assembleMatrices();

	/**
	* \brief Generate the mass matrix and its inverse
	*/
	void generateM();

	/**
	* \brief Generate the discrete Laplacian matrix
	*/
	virtual void generateL();

	/**
	* \brief Add the mass matrix to the discrete Laplacian matrix
	*/
	virtual void generateA(real alpha);

	/**
	* \brief Compute the Nth-order Taylor expansion of the matrix \c A
	*/
	void generateBN();	

	/**
	* \brief Generate the matrix containing the discrete divergence operator
	*/
	virtual void generateQT();

	/**
	* \brief Doing nothing !
	*/
	void updateQ(real gamma);
	
	/**
	* \brief Compute the matrix product \f$ C = Q^T B^N Q \f$
	*/
	void generateC();

	/**
	* \brief Generate explicit terms in the momentum equation
	*/
	void generateRN();

	/**
	* \brief Calculate explicit flux terms
	*/
	void calculateExplicitQTerms();

	/**
	* \brief Calculate explicit lambda terms
	*/
	void calculateExplicitLambdaTerms();
	
	/**
	* \brief Generate the inhomogenous boundary conditions from the discrete Laplacian operator
	*/
	virtual void generateBC1();

	/**
	* \brief Generate the inhomogeneous boundary conditions from the discrete continuity equation
	*/
	virtual void generateBC2();

	/**
	* \brief Assemble the right hand-side of the system for the intermediate flux variables
	*/
	virtual void assembleRHS1();

	/**
	* \brief Assemble the right hand-side of the system for the pressure variables and body forces
	*/
	void assembleRHS2();

	/**
	* \brief Solve for the intermediate flux variables
	*/
	virtual void solveIntermediateVelocity();

	/**
	* \brief Solve for the pressure varaibles and body forces
	*/
	void solvePoisson();

	/**
	* \brief Projection of the intermediate fluxes onto the divergence-free field
	*/
	virtual void projectionStep();

	/**
	* \brief Doing nothing !
	*/
	void updateBoundaryConditions();

	/**
	* \brief Doing nothing !
	*/
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
	* \brief Perform necessary tasks to end the simulation
	*/
	virtual void shutDown();
	
	/**
	* \brief Give the name of the current solver 
	*/
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
};