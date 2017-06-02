/**
 * \file parseSimulationFile.cu
 * \brief Parses the file \a simParams.yaml and stores the numerical 
 *        parameters used in the simulation.
 */


#include <fstream>

#include "io.h"
#include "utilities/parameterDB.h"
#include "yaml-cpp/yaml.h"


/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

/**
 * \brief Converts a string to a time-integration scheme type.
 *
 * \param s the string that describes the time-integration scheme
 *
 * \return a time-integration scheme type
 */
timeScheme timeSchemeFromString(std::string s)
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
	{
		printf("Error: Unknown timeScheme '%s'!\n", s.c_str());
		exit(-1);
	}
} // timeSchemeFromString


/**
 * \brief Converts a string to a preconditioner type.
 *
 * \param s the string that describes the preconditioner
 *
 * \return a preconditioner type
 */
preconditionerType preconditionerTypeFromString(std::string s)
{
	if (s == "NONE")
		return NONE;
	else if (s == "DIAGONAL")
		return DIAGONAL;
	else if (s == "SMOOTHED_AGGREGATION")
		return SMOOTHED_AGGREGATION;
	else if (s == "AINV")
		return AINV;
	else
	{
		printf("Error: Unknown preconditioner '%s'!\n", s.c_str());
		exit(-1);
	};
} // preconditionerTypeFromString


/**
 * \brief Converts a string to an IBM scheme.
 *
 * \param s the string that describes the IBM scheme
 *
 * \return an IBM-scheme type
 */
ibmScheme ibmSchemeFromString(std::string s)
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
	else if (s == "DIFFUSION")
		return DIFFUSION;
	else if (s == "DF_MODIFIED")
		return DF_MODIFIED;
	else if (s == "FEA_MODIFIED")
		return FEA_MODIFIED;
	else if (s == "DF_IMPROVED")
		return DF_IMPROVED;
	else
	{
		printf("Error: Unknown ibmScheme '%s'!\n", s.c_str());
		exit(-1);
	}
} // ibmSchemeFromString


/**
 * \brief Converts a string to a interpolation type.
 *
 * \param s the string that describes the type of interpolation
 *
 * \return an interpolation type
 */
interpolationType interpolationTypeFromString(std::string s)
{
	if (s == "CONSTANT")
		return CONSTANT;
	else if (s == "LINEAR")
		return LINEAR;
	else if (s == "QUADRATIC")
		return QUADRATIC;
	else
	{
		printf("Error: Unknown interpolationType '%s'!\n", s.c_str());
		exit(-1);
	};
} // interpolationTypeFromString


/**
 * \brief Fills the database with the simulation parameters.
 *
 * \param node the parsed file
 * \param DB database that contains the simulation parameters
 */
void parseSimulation(const YAML::Node &node, parameterDB &DB)
{
	// write to DB
	std::string dbKey = "simulation";
	DB[dbKey]["dt"].set<real>(node["dt"].as<real>(0.02));
	DB[dbKey]["scaleCV"].set<real>(node["scaleCV"].as<real>(2.0));
	DB[dbKey]["startStep"].set<int>(node["startStep"].as<int>(0));
	DB[dbKey]["nt"].set<int>(node["nt"].as<int>(100));
	DB[dbKey]["nsave"].set<int>(node["nsave"].as<int>(100));
	DB[dbKey]["ibmScheme"].set<ibmScheme>(
	  ibmSchemeFromString(node["ibmScheme"].as<std::string>("NAVIER_STOKES")));
	DB[dbKey]["convTimeScheme"].set<timeScheme>(
	  timeSchemeFromString(node["timeScheme"][0].as<std::string>("EULER_EXPLICIT")));
	DB[dbKey]["diffTimeScheme"].set<timeScheme>(
	  timeSchemeFromString(node["timeScheme"][1].as<std::string>("EULER_IMPLICIT")));
	DB[dbKey]["interpolationType"].set<interpolationType>(
	  interpolationTypeFromString(node["interpolationType"].as<std::string>("LINEAR")));

	const YAML::Node &solvers = node["linearSolvers"];
	for (unsigned int i=0; i<solvers.size(); i++)
	{
		std::string system = solvers[i]["system"].as<std::string>();
		std::string dbKey = system + "Solve";
		DB[dbKey]["solver"].set<std::string>(solvers[i]["solver"].as<std::string>("CG"));
		if (DB[dbKey]["solver"].get<std::string>() == "GMRES")
			DB[dbKey]["restart"].set<int>(solvers[i]["restart"].as<int>(50));
		DB[dbKey]["preconditioner"].set<preconditionerType>(
		  preconditionerTypeFromString(solvers[i]["preconditioner"].as<std::string>("DIAGONAL")));
		DB[dbKey]["rTol"].set<real>(solvers[i]["relTolerance"].as<real>(1.0E-05));
		DB[dbKey]["aTol"].set<real>(solvers[i]["absTolerance"].as<real>(1.0E-50));
		DB[dbKey]["maxIterations"].set<int>(solvers[i]["maxIterations"].as<int>(10000));
		DB[dbKey]["monitor"].set<int>(solvers[i]["monitor"].as<bool>(false));
	}
} // parseSimulation


/**
 * \brief Parses \a simParams.yaml and stores the simulation parameters.
 *
 * \param simFile the file that contains the simulation parameters
 * \param DB the database that will be filled
 */
void parseSimulationFile(std::string &simFile, parameterDB &DB)
{
	printf("Parsing YAML file with simulation parameters ...\n");
	YAML::Node nodes = YAML::LoadFile(simFile);
	for(unsigned int i=0; i<nodes.size(); i++)
		parseSimulation(nodes[i], DB);
} // parseSimulationFile

} // End of namespace io
