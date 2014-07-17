/***************************************************************************//**
* \file  TairaColoniusSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the class \c TairaColoniusSolver
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>

/**
* \class TairaColoniusSolver
* \brief Solve the flow using the Immersed boundary method described by Taira and Colonius (2007).
*
* <b>The immersed boundary method: a projection approach</b> \n
* Taira, K and Colonius, T \n
* Journal of Computational Physics \n
* Volume 225 Number 2 \n
* 2007
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
	* \brief Generate matrix containing discrete divergence operator and interpolation matrix
	*/
	void generateQT();

	/**
	* \brief Update the interpolation matrix with current location of body points
	*/
	void updateQT();

	/**
	* \brief Create an array that contains 
	*		the inhomogeneous boundary conditions from the continuity equation
	*		and the no-slip boundary conditions
	*/
	void generateBC2();
	
	/**
	* \brief Update the location of the bodies and re-compute matrices
	*/
	void updateSolverState();
	
	/**
	* \brief Calculate forces acting on each immersed body
	*/
	void calculateForce();
	
	/**
	* \brief Cannot find the definition of this method !
	*/
	void generateE();
	
public:
	/**
	* \brief Constructor of the class \c TairaColoniusSolver
	*/
	TairaColoniusSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	/**
	* \brief Initialize the solver
	*/
	void initialise();
	
	/**
	* \brief Calculate forces and write time and forces for each body in a file
	*/
	void writeData();
	
	/**
	* \brief Return the name of the solver
	*/
	std::string name()
	{
		return "Taira & Colonius";
	}
};
