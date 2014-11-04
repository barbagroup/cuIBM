/***************************************************************************//**
 * \file  TairaColoniusSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c TairaColoniusSolver.
 */


#pragma once

#include "NSWithBody.h"


/**
 * \class TairaColoniusSolver
 * \brief Solve the flow using the Immersed boundary method described by Taira and Colonius (2007).
 *
 * <b>The immersed boundary method: a projection approach</b> \n
 * Taira K. and Colonius T. \n
 * Journal of Computational Physics \n
 * Volume 225 Number 2 (2007).
 *
 */
template <typename memoryType>
class TairaColoniusSolver : public NSWithBody<memoryType>
{
private:
	cusp::coo_matrix<int, real, memoryType> 
		E,		///< Interpolation matrix from the Eulerian grid to the Lagrangian points
		ET;		///< Regularization matrix form the Lagrangian points to the Eulerian grid
	
	/**
	 * \brief Generates matrix containing discrete divergence operator and interpolation matrix.
	 */
	virtual void generateQT();

	/**
	 * \brief Updates the interpolation matrix with current location of body points.
	 */
	void updateQT();

	/**
	 * \brief Creates an array that contains 
	 *        the inhomogeneous boundary conditions from the continuity equation
	 *        and the no-slip boundary conditions.
	 */
	virtual void generateBC2();
	
	/**
	 * \brief Updates the location of the bodies and re-compute matrices.
	 */
	virtual void updateSolverState();
	
	/**
	 * \brief Calculates forces acting on each immersed body.
	 */
	virtual void calculateForce();
	
	/**
	 * \brief Cannot find the definition of this method !
	 */
	void generateE();
	
public:
	/**
	 * \brief Constructor of the class \c TairaColoniusSolver.
	 */
	TairaColoniusSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	/**
	 * \brief Initializes the solvers.
	 */
	virtual void initialise();
	
	/**
	 * \brief Calculates forces and write time and forces for each body in a file.
	 */
	virtual void writeData();
	
	/**
	 * \brief Returns the name of the solver.
	 */
	virtual std::string name()
	{
		return "Taira & Colonius";
	}
};
