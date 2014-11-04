/***************************************************************************//**
 * \file parameterDB.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c property.
 */


#pragma once

#include <cstdio>
#include <typeinfo>
#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include "boundaryCondition.h"


/**
 * \class property
 * \brief Contains information about a property stored in a generic way.
 */
class property
{
	template <typename T> T getInternal();
	
	// hack to get around std::type_info not having a default constructor
	const std::type_info *type;
	
	/// 256 bytes of storage for the name of the folder
	char value[256];

public:
	/**
	 * \brief Constructor of the class \c property.
	 *
	 * It initializes the memory to zero.
	 *
	 */
	property()
	{
		memset(value, 0, 256);
	}
	
	/**
	 * \brief Gets the value of the property as a given type.
	 */
	template <typename T> T get();
	
	/**
	 * \brief Sets the value of the property given a type.
	 */
	template <typename T> void set(T v);

	/**
	 * \brief Returns a string describing the value of property as appropriate.
	 */
	const char *print();
};

/**
 * \typedef componentParameter
 * \brief Map from a string to a \c property.
 */
typedef std::map<std::string, property> componentParameter;

/**
 * \typedef parameterDB
 * \brief Map from a string to a \c componentParameter.
 */
typedef std::map<std::string, componentParameter> parameterDB;

// output a DB
/*
void printDB(parameterDB &DB)
{
  for (parameterDB::iterator it=DB.begin(); it!=DB.end(); ++it)
  {
    printf("%s:\n",it->first.c_str());
    for (componentParameter::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2)
      printf("\t%s: %s\n",it2->first.c_str(),it2->second.print());
  }
}
*/
