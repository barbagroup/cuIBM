/**
 * \file boundaryCondition.h
 * \brief Definition of the class \c boundaryCondition.
 */


#pragma once

#include <string>
#include <sstream>
#include "types.h"
#include "parameterDB.h"


/**
 * \class boundaryCondition
 * \brief Stores the type of boundary condition and its value.
 */
class boundaryCondition
{
public:
	bcType type; ///< type of boundary condition
	real value;  ///< numerical value associated with the boundary condition
	
	/**
	 * \brief Constructor. Sets Dirichlet boundary condition with value.
	 */
	boundaryCondition() : type(DIRICHLET), value(0) {};

	/**
	 * \brief Constructor overloading. Sets boundary condition to a given type
	 *        with a given value.
	 *
	 * \param _type the type of boundary condition
	 * \param _value the value at the boundary
	 */
	boundaryCondition(bcType _type, real _value) : type(_type), value(_value) {}; 

};
