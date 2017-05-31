/**
 * \file parseFlowFile.cu
 * \brief Parse the input file \a flow.yaml to get boundary and initial 
 *        conditions of the flow.
 */


#include <fstream>

#include "io.h"
#include "yaml-cpp/yaml.h"


/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

/**
 * \brief Converts a string to a boundary location type.
 *
 * \param s the string that describes the location of the boundary
 *
 * \return a boundary location type
 */
boundary boundaryFromString(std::string s)
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
		printf("Error: Unknown location '%s'!\n", s.c_str());
		exit(-1);
	}
} // boundaryFromString


/**
 * \brief Converts string to a boundary condition type.
 *
 * \param s the string that describes the condition to impose on the boundary
 *
 * \return a boundary condition type
 */
bcType bcTypeFromString(std::string s)
{
	if (s == "DIRICHLET")
		return DIRICHLET;
	else if (s == "NEUMANN")
		return NEUMANN;
	else if (s == "CONVECTIVE")
		return CONVECTIVE;
	else if (s == "PERIODIC")
		return PERIODIC;
	else if (s == "SPECIAL")
		return SPECIAL;
	else
	{
		printf("Error: Unknown boundary type '%s'!\n", s.c_str());
		exit(-1);
	}
} // bcTypeFromString


/**
 * \brief Fills the database with the flow parameters.
 *
 * \param node the parsed file
 * \param DB the database to be filled
 */
void parseFlow(const YAML::Node &node, parameterDB &DB)
{
	std::string dbKey = "flow";
	DB[dbKey]["nu"].set<real>(node["nu"].as<real>(0.01));
	DB[dbKey]["uInitial"].set<real>(node["initialVelocity"][0].as<real>(1.0));
	DB[dbKey]["vInitial"].set<real>(node["initialVelocity"][1].as<real>(0.0));
	DB[dbKey]["uPerturb"].set<real>(node["initialPerturbation"][0].as<real>(0.0));
	DB[dbKey]["vPerturb"].set<real>(node["initialPerturbation"][1].as<real>(0.0));

	boundaryCondition **bc = 0;
	bc = DB[dbKey]["boundaryConditions"].get<boundaryCondition **>();
	if (!bc)
	{
		printf("Error: BoundaryConditions pointer not initialized.\n");
		exit(-1);
	}
	const YAML::Node &BCs = node["boundaryConditions"];
	for (unsigned int i=0; i<BCs.size(); i++)
	{
		boundary loc = boundaryFromString(BCs[i]["location"].as<std::string>());
		bc[loc][0] = boundaryCondition(
			bcTypeFromString(BCs[i]["u"][0].as<std::string>()), BCs[i]["u"][1].as<real>());
		bc[loc][1] = boundaryCondition(
			bcTypeFromString(BCs[i]["v"][0].as<std::string>()), BCs[i]["v"][1].as<real>());
	}
} // parseFlow


/**
 * \brief Parses the \a flow file and stores the parameters in the database.
 *
 * \param flowFile the file that contains information about the flow
 * \param DB the database that will be filled
 */
void parseFlowFile(std::string &flowFile, parameterDB &DB)
{
	YAML::Node nodes = YAML::LoadFile(flowFile);
	for (unsigned int i=0; i<nodes.size(); i++)
		parseFlow(nodes[i], DB);
} // parseFlowFile

} // End of namespace io
