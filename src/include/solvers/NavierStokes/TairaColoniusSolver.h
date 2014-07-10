/***************************************************************************//**
* \file  TairaColoniusSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Solves the flow using the IB method described by Taira and Colonius (2007)
*
* <b>The immersed boundary method: a projection approach</b> \n
* Taira, K and Colonius, T \n
* Journal of Computational Physics \n
* Volume 225 Number 2 \n
* 2007
*
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>

/********************//**
* \class TairaColoniusSolver
* \brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class TairaColoniusSolver : public NSWithBody<memoryType>
{
private:
	cusp::coo_matrix<int, real, memoryType> E,		///< Interpolation matrix from the Eulerian grid to the Lagrangian points
											ET;		///< Regularization matrix form the Lagrangian points to the Eulerian grid
	

	/**
	* \brief Generate the sub-matrix containing the discrete divergence operator and the interpolation matrix
	*/
	void generateQT();

	void updateQT();

	void generateBC2();
	
	void updateSolverState();
	
	void calculateForce();
	
	void generateE();
	
public:
	/**
	* \brief Constructor of the class TairaColoniusSolver
	*/
	TairaColoniusSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	void initialise();
	
	void writeData();
	
	/**
	* \brief Return the name of the solver
	*/
	std::string name()
	{
		return "Taira & Colonius";
	}
};
