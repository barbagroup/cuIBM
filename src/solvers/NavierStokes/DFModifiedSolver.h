/**
* @file  DFModifiedSolver.h
* @brief Modified direct forcing immersed boundary method.
*/

#pragma once

#include "DirectForcingSolver.h"
/**
* The direct forcing method proposed by Fadlun et al (2000)
* with modifications by Kim et al (2001)
*/
template <typename memoryType>
class DFModifiedSolver : public DirectForcingSolver<memoryType>
{
protected:
	virtual void calculateExplicitLambdaTerms();
	virtual void projectionStep();
	
public:
	DFModifiedSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "DFModified";
	}
};