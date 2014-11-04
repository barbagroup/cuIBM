/***************************************************************************//**
 * \file boundaryCondition.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c boundaryCondition.
 */


#pragma once

#include <string>
#include <sstream>
#include "types.h"
#include "parameterDB.h"


/**
 * \class boundaryCondition
 * \brief Stores the boundary conditions for a given system.
 */
class boundaryCondition
{
public:
	bcType type; ///< type of boundary condition
	real  value; ///< numerical value associated with the boundary condition
	
	/**
	 * \brief Constructor of the class \c boundaryCondition.
	 *
	 * Boundary condition initialized with a Dirichlet-type with
	 * with a value sets to zero.
	 *
	 */
	boundaryCondition() : type(DIRICHLET), value(0) {};

	/**
	 * \brief Other constructor of the class \c boundaryCondition.
	 *
	 * Boundary condition initialized with a given type and a given value.
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
