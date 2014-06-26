/***************************************************************************//**
* \file parameterDB.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the class \c property \c
*/

#include <parameterDB.h>
#include <body.h>
#include <boundaryCondition.h>

// converts a number to a string
template <typename T>
std::string toString(T num)
{
  std::stringstream ss; 
  ss << num;
  return ss.str();
}

// converts a type of boundary condition to a string
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

// returns the value of the property as a string
template <>
std::string property::get()
{
  return std::string(value);
}

// returns the value of the property as a given type
template <typename T>
T property::get()
{
  T r = *reinterpret_cast<T*>(&value[0]);
  return r;
}

// explicit instantiations of the method property::get()
// to be able to define the template function outside the header file
template double property::get<double>();
template float property::get<float>();
template int property::get<int>();
template bool property::get<bool>();
template preconditionerType property::get<preconditionerType>();
template timeScheme property::get<timeScheme>();
template ibmScheme property::get<ibmScheme>();
template std::vector<body> *property::get<std::vector<body>*>();
template boundaryCondition **property::get<boundaryCondition **>();

// returns the value of the property as a string
const char *property::print()
{ 
    if (*type == typeid(int))
      return toString(this->get<int>()).c_str();
    else if (*type == typeid(double))
      return toString(this->get<double>()).c_str();
	else if (*type == typeid(float))
      return toString(this->get<float>()).c_str();
    else if (*type == typeid(std::string))
      return (const char *)value;
    else
      return "not found";
}

// sets the value of the property given a string
template <>
void property::set(std::string s)
{
  strncpy(value,s.c_str(),64);
}

// sets the value of the property given the type
template <typename T>
void property::set(T v)
{
  // assume we have enough space (64 bytes)
  type = &typeid(T);
  *reinterpret_cast<T*>(&value[0]) = v;
}

// explicit instantiations of the method property::set()
// to be able to define the template function outside the header file
template void property::set(int v);
template void property::set(float v);
template void property::set(double v);
template void property::set(bool v);
template void property::set(timeScheme v);
template void property::set(ibmScheme v);
template void property::set(preconditionerType v);
template void property::set(boundaryCondition **v);
template void property::set(std::vector<body> *v);
