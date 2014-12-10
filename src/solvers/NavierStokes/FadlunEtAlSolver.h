/**
* @file  FadlunEtAlSolver.h
* @brief Modified direct forcing immersed boundary method.
*/

#pragma once

#include "DirectForcingSolver.h"
#include "DFModifiedSolver.h"
/**
* The direct forcing method proposed by Fadlun et al (2000)
* with modifications by Kim et al (2001)
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
	FadlunEtAlSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
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
