#pragma once

#include <cstdio>
#include <typeinfo>
#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include <boundaryCondition.h>

// property.h
// generic property storage
class property
{
public:
	property()
	{
		// initialise the memory to 0
		memset(value,0,64);
	}
	// get a value given a type
	template <typename T> T get();
	// set a value given a type
	template <typename T> void set(T v);

	// return string describing value of property as appropriate
	const char *print();
private:
	template <typename T> T getInternal();
	// hack to get around std::type_info not having a default constructor
	const std::type_info *type;
	// 64 bytes of storage
	char value[64];
};

typedef std::map<std::string, property> componentParameter;
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

