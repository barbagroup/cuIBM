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
	cusp::coo_matrix<int, real, host_memory> CHost;

	//virtual void generateQT();
	//void updateQT();
	virtual void generateC();
	
public:
	DFImprovedSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "Direct Forcing Improved";
	}
};
