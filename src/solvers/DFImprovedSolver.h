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
  // generate the modified divergence operator
	virtual void generateQT();
	
public:
  // constructor -- get simulation parameters and grid
	DFImprovedSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);

  /**
   * \brief Return the name of solver.
   */
	virtual std::string name()
	{
		return "Direct Forcing Improved";
	}

}; // DFImprovedSolver
