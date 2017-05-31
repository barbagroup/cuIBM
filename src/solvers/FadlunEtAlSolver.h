/**
 * \file  FadlunEtAlSolver.h
 * \brief Declaration of the classes \c FadlunEtAlSolver and \c FEAModifiedSolver.
 */


#pragma once

#include "DirectForcingSolver.h"
#include "DFModifiedSolver.h"


/**
 * \class FadlunEtAlSolver
 * \brief Direct forcing method proposed by Fadlun et al (2000).
 *
 * This is the same method that Fadlun et al. used in their 2000 paper.
 * An explicit pressure is used in the intermediate velocity step.
 * The quantity obtained in the Poisson equation is the pressure correction,
 * and has to be added to the pressure at the previous step
 * to obtain the pressure at the next time step.
 *
 * Inherits from DFModifiedSolver.
 *
 * The Poisson equations solved the standard Poisson system for the entire domain.
 * The divergence formula is the same in every cell.
 * But the equation for the intermediate velocity involves both 
 * the momentum equation as well as the interpolation relations.
 * Hence, both G and Q need to be defined separately, with G missing rows that
 * correspond to the forcing nodes.
 */
template <typename memoryType>
class FadlunEtAlSolver : public DFModifiedSolver<memoryType>
{
private:
	cusp::coo_matrix<int, real, memoryType> G;  ///< gradient operator

	// update the gradient operator
	virtual void updateG();

	// compute the gradient and divergence operators
	virtual void generateQT();

	// add explicit pressure gradient to RHS of velocity system
	virtual void calculateExplicitLambdaTerms();

	// compute the matrix Poisson system
	virtual void generateC();
	
public:
	// constructor -- copy database and grid
	FadlunEtAlSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	/**
	 * \brief Returns the name of the solver as a string.
	 */
	virtual std::string name()
	{
		return "Fadlun et al.";
	}
}; // FadlunEtAlSolver


/**
 * \class FEAModifiedSolver
 * \brief Direct forcing method proposed by Fadlun et al (2000), without 
 *        explicit pressure added when the intermediate velocity is solved.
 *
 * Detailed description is similar to that of FadlunEtAlSolver.
 * The difference is that no explicit pressure term is added when the
 * momentum equation is solved, and so it inherits from DirectForcingSolver.
 */
template <typename memoryType>
class FEAModifiedSolver : public DirectForcingSolver<memoryType>
{
private:
	// compute the gradient and divergence operators
	virtual void generateQT();
	// compute the matrix of the Poisson system
	virtual void generateC();
	
public:
	// constructor -- copy database and grid
	FEAModifiedSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	/**
	 * \brief Returns the name of the solver as a string.
	 */
	virtual std::string name()
	{
		return "FEAModified";
	}
}; // FEAModifiedSolver
