/***************************************************************************//**
* \file boundaryCondition.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the class \c boundaryCondition
*/

#pragma once

#include <string>
#include <sstream>
#include "types.h"
#include "parameterDB.h"

/**
* \class boundaryCondition
* \brief Store the boundary conditions for a given system
*/
class boundaryCondition
{
public:
	bcType type; ///< type of boundary condition
	real  value; ///< numerical value associated with the boundary condition
	
	/**
	* \brief Constructor of the class \c boundaryCondition.
	*
	* Initialize with a Dirichlet-type boundary condition
	* with a value sets to zero.
	*
	*/
	boundaryCondition() : type(DIRICHLET), value(0) {};

	/**
	* \brief Other constructor of the class \c boundaryCondition.
	*
	* Initialize with a given boundary condition type
	* and a given value.
	*
	*/
	boundaryCondition(bcType _type, real _value) : type(_type), value(_value) {}; 

  /*const char *print()
  {
    std::stringstream ss; 
    ss << toString(this->type);
    ss << " : ";
    ss << this->value;
    std::string st = ss.str();
    //std::cout << st << std::endl;
    return ss.str().c_str();
  }*/
};
