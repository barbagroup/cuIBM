/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

#include <parameterDB.h>
#include <body.h>
#include <boundaryCondition.h>

/// convert number to string
template <typename T>
std::string toString(T num)
{
  std::stringstream ss; 
  ss << num;
  return ss.str();
}

/// obtain a string for the boundary conditions
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

template <>
std::string property::get()
{
  return std::string(value);
}

template <typename T>
T property::get()
{
  T r = *reinterpret_cast<T*>(&value[0]);
  return r;
}

template double property::get<double>();
template float property::get<float>();
template int property::get<int>();
template bool property::get<bool>();
template preconditionerType property::get<preconditionerType>();
template timeScheme property::get<timeScheme>();
template ibmScheme property::get<ibmScheme>();
template std::vector<body> *property::get<std::vector<body>*>();
template boundaryCondition **property::get<boundaryCondition **>();

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


template <>
void property::set(std::string s)
{
  strncpy(value,s.c_str(),64);
}

template <typename T>
void property::set(T v)
{
  // assume we have enough space (64 bytes)
  type = &typeid(T);
  *reinterpret_cast<T*>(&value[0]) = v;
}

template void property::set(int v);
template void property::set(float v);
template void property::set(double v);
template void property::set(bool v);
template void property::set(timeScheme v);
template void property::set(ibmScheme v);
template void property::set(preconditionerType v);
template void property::set(boundaryCondition **v);
template void property::set(std::vector<body> *v);

