/**
* @file  parameterDB.h
* @brief Database of all the simluation parameters
*/

#pragma once

#include <cstdio>
#include <typeinfo>
#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include "boundaryCondition.h"

// generic property storage
/**
* @class property
*/
class property
{
	template <typename T> T getInternal();
	
	// hack to get around std::type_info not having a default constructor
	const std::type_info *type;
	
	/// 256 bytes of storage for the name of the folder
	char value[256];

public:
	/// Constructor. Intialises the memory to zero.
	property()
	{
		memset(value, 0, 256);
	}
	
	/// Get a value given a type
	template <typename T> T get();
	
	/**
	* @brief Set a value given a type
	*
	* @param v The value
	*/
	template <typename T> void set(T v);

	// return string describing value of property as appropriate
	const char *print();
};

/**
* @typedef componentParameter
* @typedef Maps from strings to a property 
*/
typedef std::map<std::string, property> componentParameter;

/**
* @typedef parameterDB
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

