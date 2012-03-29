#pragma once

#include <cstdio>
#include <typeinfo>
#include <string>
#include <map>
#include <sstream>
#include <iostream>

template <typename T>
std::string toString(T num)
{
  std::stringstream ss;
  ss << num;
  return ss.str();
}

template <>
std::string toString(bcType b)
{
  if (b == DIRICHLET)
    return "Dirichlet";
  else if (b == NEUMANN)
    return "Neumann";
  else if (b == CONVECTIVE)
    return "Convective";
  else
    return "Error";
}

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
  const char *print()
  {
    if (*type == typeid(int))
      return toString(this->get<int>()).c_str();
    else if (*type == typeid(double))
      return toString(this->get<double>()).c_str();
    else if (*type == typeid(std::string))
      return (const char *)value;
    else if (*type == typeid(boundaryCondition))
      return reinterpret_cast<boundaryCondition*>(value)->print();
    else
      return "not found";
  }
private:
  // hack to get around std::type_info not having a default constructor
  const std::type_info *type;
  // 64 bytes of storage
  char value[64];
};

template <typename T>
T property::get()
{
  T r = *reinterpret_cast<T*>(&value[0]);
  return r;
}

template <>
std::string property::get()
{
  return std::string(value);
}

template <typename T>
void property::set(T v)
{
  // assume we have enough space (64 bytes)
  type = &typeid(T);
  *reinterpret_cast<T*>(&value[0]) = v;
}

template <>
void property::set(std::string s)
{
  strncpy(value,s.c_str(),64);
}

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

