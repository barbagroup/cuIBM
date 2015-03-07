/**
* @file  DFImprovedSolver.h
* @brief Second-order direct forcing method.
*/

#pragma once

#include "DirectForcingSolver.h"
/**
* Second-order fully-discrete direct forcing method
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
