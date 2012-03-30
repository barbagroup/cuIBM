#include <io/io.h>
#include <fstream>
#include <yaml-cpp/yaml.h>

namespace io
{

using std::string;

boundary boundaryFromString(string &s)
{
  if (s == "xMinus")
    return XMINUS;
  else if (s == "xPlus")
    return XPLUS;
  else if (s == "yMinus")
    return YMINUS;
  else if (s == "yPlus")
    return YPLUS;
  else
  {
    printf("[E]: invalid boundary condition location defined\n");
    exit(0);
  }
}

bcType bcTypeFromString(string &s)
{
  if (s == "DIRICHLET")
    return DIRICHLET;
  else if (s == "NEUMANN")
    return NEUMANN;
  else if (s == "CONVECTIVE")
    return CONVECTIVE;
  else if (s == "PERIODIC")
    return PERIODIC;
  else
  {
    printf("[E]: invalid boundary condition type specified\n");
    exit(0);
  }
}

void parseFlow(const YAML::Node &node, parameterDB &DB)
{
  string dbKey = "flow";
  double nu = 0.01, initialU = 1, initialV = 0;
  node["nu"] >> nu;
  node["initialVelocity"][0] >> initialU;
  node["initialVelocity"][1] >> initialV;

  DB[dbKey]["nu"].set<double>(nu);
  DB[dbKey]["uInitial"].set<double>(initialU);
  DB[dbKey]["vInitial"].set<double>(initialV);

  boundaryCondition **bc = 0;
  DB[dbKey]["boundaryConditions"].get<boundaryCondition **>();
  const YAML::Node &BCs = node["boundaryConditions"];
  for (unsigned int i=0; i<BCs.size(); i++)
  {
    string location, uType, vType;
    double uVal, vVal;
    BCs[i]["location"] >> location;
    BCs[i]["u"][0] >> uType;
    BCs[i]["v"][0] >> vType;
    BCs[i]["u"][1] >> uVal;
    BCs[i]["v"][1] >> vVal;

    boundary loc = boundaryFromString(location);
    bc[loc][0] = boundaryCondition(bcTypeFromString(uType),uVal);
    bc[loc][1] = boundaryCondition(bcTypeFromString(vType),vVal);

    // printf("%s: %s, %lg : %s, %lg\n",location.c_str(),uType.c_str(),uVal,vType.c_str(),vVal);
  }
}

void parseFlowFile(std::string &flowFile, parameterDB &DB)
{
  std::ifstream fin(flowFile.c_str());
  YAML::Parser parser(fin);
  YAML::Node doc;
  parser.GetNextDocument(doc);

  for (unsigned int i=0; i<doc.size(); i++)
    parseFlow(doc[i],DB);
}

} // end namespace io
