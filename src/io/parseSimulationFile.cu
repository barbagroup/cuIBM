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
#include <parameterDB.h>
#include <yaml-cpp/yaml.h>
#include <fstream>

namespace io
{

using std::string;

timeScheme timeSchemeFromString(string &s)
{
  if (s == "EULER_EXPLICIT")
    return EULER_EXPLICIT;
  else if (s == "EULER_IMPLICIT")
    return EULER_IMPLICIT;
  else if (s == "ADAMS_BASHFORTH_2")
    return ADAMS_BASHFORTH_2;
  else if (s == "RUNGE_KUTTA_3")
    return RUNGE_KUTTA_3;
  else if (s == "CRANK_NICOLSON")
    return CRANK_NICOLSON;
  else
    return EULER_EXPLICIT;
}

preconditionerType preconditionerTypeFromString(string &s)
{
  if (s == "NONE")
    return NONE;
  else if (s == "DIAGONAL")
    return DIAGONAL;
  else if (s == "SMOOTHED_AGGREGATION")
    return SMOOTHED_AGGREGATION;
  else
    return NONE;
}

ibmScheme ibmSchemeFromString(string &s)
{
  if (s == "NAVIER_STOKES")
    return NAVIER_STOKES;
  else if (s == "SAIKI_BIRINGEN")
    return SAIKI_BIRINGEN;
  else if (s == "FADLUN_ET_AL")
    return FADLUN_ET_AL;
  else if (s == "TAIRA_COLONIUS")
    return TAIRA_COLONIUS;
  else
    return NAVIER_STOKES;
}

void parseSimulation(const YAML::Node &node, parameterDB &DB)
{
	real dt = 0.02;
	int nt = 100, nsave = 100, startStep = 0;
	string convSch = "EULER_EXPLICIT", diffSch = "EULER_IMPLICIT", ibmSch = "NAVIER_STOKES";
	bool restart = false;

	node["dt"] >> dt;
	node["nsave"] >> nsave;
	node["nt"] >> nt;
	node["restart"] >> restart;
	node["startStep"] >> startStep;
	node["ibmScheme"] >> ibmSch;
	node["timeScheme"][0] >> convSch;
	node["timeScheme"][1] >> diffSch;
	node["ibmScheme"] >> ibmSch;

	// write to DB
	string dbKey = "simulation";
	DB[dbKey]["dt"].set<real>(dt);
	DB[dbKey]["nsave"].set<int>(nsave);
	DB[dbKey]["nt"].set<int>(nt);
	DB[dbKey]["restart"].set<bool>(restart);

	DB[dbKey]["convTimeScheme"].set<timeScheme>(timeSchemeFromString(convSch));
	DB[dbKey]["diffTimeScheme"].set<timeScheme>(timeSchemeFromString(diffSch));

	DB[dbKey]["ibmScheme"].set<ibmScheme>(ibmSchemeFromString(ibmSch));

	const YAML::Node &solvers = node["linearSolvers"];
	for (unsigned int i=0; i<solvers.size(); i++)
	{
		string system = "velocity", solver = "CG", PC = "DIAGONAL";
		real tol = 1e-5;
		int maxIter = 10000;
		solvers[i]["system"] >> system;
		solvers[i]["solver"] >> solver;
		solvers[i]["preconditioner"] >> PC;
		solvers[i]["tolerance"] >> tol;
		solvers[i]["maxIterations"] >> maxIter;

		string dbKey;
		if (system == "velocity")
			dbKey = "velocitySolve";
		else if (system == "Poisson")
			dbKey = "PoissonSolve";
		else
		{
			printf("[E]: Unknown solver system found in io::parseSimulation\n");
			exit(0);
		}

		// write to DB
		DB[dbKey]["solver"].set<string>(solver);
		DB[dbKey]["preconditioner"].set<preconditionerType>(preconditionerTypeFromString(PC));
		DB[dbKey]["tolerance"].set<real>(tol);
		DB[dbKey]["maxIterations"].set<int>(maxIter);
	}
}

void parseSimulationFile(std::string &simFile, parameterDB &DB)
{
	std::ifstream fin(simFile.c_str());
	YAML::Parser parser(fin);
	YAML::Node doc;
	parser.GetNextDocument(doc);

	for(unsigned int i=0; i<doc.size(); i++)
		parseSimulation(doc[i], DB);
}

} // end namespace io
