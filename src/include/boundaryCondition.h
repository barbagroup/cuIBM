#pragma once

#include <types.h>
#include <string>
#include <sstream>
#include <parameterDB.h>

class boundaryCondition
{
public:
  bcType type;
  double value;

  boundaryCondition() : type(DIRICHLET), value(0) {}; 
  boundaryCondition(bcType _type, double _value) : type(_type), value(_value) {}; 


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
