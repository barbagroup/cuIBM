/**
* \file  DFImprovedSolver.h
* \brief Definition of the class \c DFImprovedSolver.
*/


#pragma once

#include "DirectForcingSolver.h"


/**
 * \class DFImprovedSolver
 * \brief Second-order fully-discrete direct forcing method.
 */
template <typename memoryType>
class DFImprovedSolver : public DirectForcingSolver<memoryType>
{
private:
	virtual void generateQT();
	
public:
	DFImprovedSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "Direct Forcing Improved";
	}
};
