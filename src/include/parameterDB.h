/***************************************************************************//**
* \file parameterDB.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Database that contains all the simulation parameters
*/

#pragma once

#include <cstdio>
#include <typeinfo>
#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include <boundaryCondition.h>

/**
* \class property
* \brief Generic property storage
*/
class property
{
	template <typename T> T getInternal();
	
	// hack to get around std::type_info not having a default constructor
	const std::type_info *type;
	
	char value[64]; ///< 64 bytes of storage for the name of the value

public:
	/**
	* \brief Constructor of the class \c property.
	*
	* Initialize the memory to zero.
	*
	*/
	property()
	{
		memset(value, 0, 64);
	}
	
	/**
	* \brief Get the value of the property as a given type
	*/
	template <typename T> T get();
	
	/**
	* \brief Set the value of the property given a type
	*/
	template <typename T> void set(T v);

	/**
	* \brief Return a string describing the value of property as appropriate
	*/
	const char *print();
};

/**
* \typedef componentParameter
* \brief Map from a string to a \c property
*/
typedef std::map<std::string, property> componentParameter;

/**
* \typedef parameterDB
* \brief Map from a string to a \c componentParameter
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
