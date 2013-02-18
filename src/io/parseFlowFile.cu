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
	real nu = 0.01, initialU = 1, initialV = 0;
	real uPerturb = 0, vPerturb = 0;
	node["nu"] >> nu;
	node["initialVelocity"][0] >> initialU;
	node["initialVelocity"][1] >> initialV;
	node["initialPerturbation"][0] >> uPerturb;
	node["initialPerturbation"][1] >> vPerturb;

	DB[dbKey]["nu"].set<real>(nu);
	DB[dbKey]["uInitial"].set<real>(initialU);
	DB[dbKey]["vInitial"].set<real>(initialV);
	DB[dbKey]["uPerturb"].set<real>(uPerturb);
	DB[dbKey]["vPerturb"].set<real>(vPerturb);

	boundaryCondition **bc = 0;
	bc = DB[dbKey]["boundaryConditions"].get<boundaryCondition **>();
	
	if (!bc)
	{
		printf("[E]: BoundaryConditions pointer not initialised\n");
		exit(-1);
	}
	const YAML::Node &BCs = node["boundaryConditions"];
	string location, uType, vType;
	real   uVal, vVal;
	for (unsigned int i=0; i<BCs.size(); i++)
	{
		BCs[i]["location"] >> location;
		BCs[i]["u"][0] >> uType;
		BCs[i]["v"][0] >> vType;
		BCs[i]["u"][1] >> uVal;
		BCs[i]["v"][1] >> vVal;

		boundary loc = boundaryFromString(location);
		bc[loc][0] = boundaryCondition(bcTypeFromString(uType), uVal);
		bc[loc][1] = boundaryCondition(bcTypeFromString(vType), vVal);

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
		parseFlow(doc[i], DB);
}

} // end namespace io
