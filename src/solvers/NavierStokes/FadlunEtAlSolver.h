/***************************************************************************//**
 * \file  FadlunEtAlSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c FadlunEtAlSolver.
 */


#pragma once

#include "DirectForcingSolver.h"
#include "DFModifiedSolver.h"
/**
 * \class FadlunEtAlSolver
 * \brief Direct forcing method proposed by Fadlun et al (2000).
 *
 * This is the same method that Fadlen et al. used in their 2000 paper.
 * An explicit pressure is used in the intermediate velocity step.
 * The quantity obtained in the Poisson equation is the pressure correction,
 * and has to be added to the pressure at the previous step
 * to obtain the pressure at the next time step.
 */
template <typename memoryType>
class FadlunEtAlSolver : public DFModifiedSolver<memoryType>
{
private:
	cusp::coo_matrix<int, real, memoryType> G;

	virtual void updateG();
	virtual void generateQT();
	virtual void calculateExplicitLambdaTerms();
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
};

//----------------------------------------------------------------------------//

template <typename memoryType>
class FEAModifiedSolver : public DirectForcingSolver<memoryType>
{
private:
	virtual void generateQT();
	virtual void generateC();
	
public:
	FEAModifiedSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "FEAModified";
	}
};
