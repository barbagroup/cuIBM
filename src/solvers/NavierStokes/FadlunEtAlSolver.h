/***************************************************************************//**
 * \file  FadlunEtAlSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c FadlunEtAlSolver.
 */


#pragma once

#include "DirectForcingSolver.h"


/**
 * \class FadlunEtAlSolver
 * \brief Direct forcing method proposed by Fadlun et al (2000).
 *
 * with modifications by Kim et al (2001).
 *
 */
template <typename memoryType>
class FadlunEtAlSolver : public DirectForcingSolver<memoryType>
{
private:
	virtual void generateQT();
	
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
