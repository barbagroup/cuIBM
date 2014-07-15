/***************************************************************************//**
* \file parseSimulationFile.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Parse the \a simulation file to fill the databse
*/

#include <io/io.h>
#include <parameterDB.h>
#include <yaml-cpp/yaml.h>
#include <fstream>

/**
* \namespace io
* \brief Contain functions related to I/O tasks
*/
namespace io
{

using std::string;

/**
* \brief Convert a string to a time-scheme type
*
* \param s the string that describes the time-scheme
*
* \return a time-scheme type
*/
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

/**
* \brief Convert a string to a prconditioner type
*
* \param s the string that describes the preconditioner
*
* \return a preconditioner type
*/
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

/**
* \brief Convert a string to an IBM scheme
*
* \param s the string that describes the IBM scheme
*
* \return an IBM-scheme type
*/
ibmScheme ibmSchemeFromString(string &s)
{
  if (s == "NAVIER_STOKES")
    return NAVIER_STOKES;
  else if (s == "SAIKI_BIRINGEN")
    return SAIKI_BIRINGEN;
  else if (s == "DIRECT_FORCING")
    return DIRECT_FORCING;
  else if (s == "TAIRA_COLONIUS")
    return TAIRA_COLONIUS;
  else if (s == "FADLUN_ET_AL")
    return FADLUN_ET_AL;
  else
    return NAVIER_STOKES;
}

/**
* \brief Fill the database with the simulation parameters
*
* \param node the parsed file
* \param DB database that contains the simulation parameters
*/
void parseSimulation(const YAML::Node &node, parameterDB &DB)
{
	real dt = 0.02, scaleCV = 5.0;
	int nt = 100, nsave = 100, startStep = 0;
	string convSch = "EULER_EXPLICIT", diffSch = "EULER_IMPLICIT", ibmSch = "NAVIER_STOKES";
	bool restart = false;

	node["dt"] >> dt;
	node["scaleCV"] >> scaleCV;
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
	DB[dbKey]["scaleCV"].set<real>(scaleCV);
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

/**
* \brief Parse the \a simulation file using YAML
*
* \param simFile the file that contains the simulation parameters
* \param DB the database that will be filled with the simulation parameters
*/
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
