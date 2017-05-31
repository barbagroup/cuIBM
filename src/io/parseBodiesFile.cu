/**
 * \file parseBodiesFile.cu
 * \brief Parses the file \a bodies.yaml to get information about immersed bodies.
 */


#include <fstream>
#include <vector>

#include "io.h"
#include "utilities/body.h"
#include "yaml-cpp/yaml.h" 


/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

/**
 * \brief Overloads the operator >>. Stores information about an immersed body.
 *
 * \param node the parsed file
 * \param Body instance of the class \c body to be filled
 * \param DB database with parameters of the simulation
 */
void parseBodiesNode(const YAML::Node &node, body &Body, parameterDB &DB)
{
	for (int i=0; i<2; i++)
		Body.X0[i] = node["centerRotation"][i].as<real>(0.0);
	// initial configuration
	for (int i=0; i<2; i++)
	{
		Body.Xc0[i] = node["initialOffset"][i].as<real>(0.0);
		Body.Xc[i] = Body.Xc0[i];
	}
	// initial angle of attack
	Body.Theta0 = node["angleOfAttack"].as<real>(0.0) * M_PI / 180.0;
	Body.Theta = Body.Theta0;
	// moving flags
	for (int i=0; i<2; i++)
		Body.moving[i] = node["moving"][i].as<bool>(false);
	// velocity
	for (int i=0; i<2; i++)
		Body.velocity[i] = node["velocity"][i].as<real>(0.0);
	// omega
	Body.omega = node["omega"].as<real>(0.0);
	// oscillation in X
	for (int i=0; i<3; i++)
		Body.xOscillation[i] = node["xOscillation"][i].as<real>(0.0);
	Body.xOscillation[1] *= 2 * M_PI;
	// oscillation in Y
	for (int i=0; i<3; i++)
		Body.yOscillation[i] = node["yOscillation"][i].as<real>(0.0);
	Body.yOscillation[1] *= 2 * M_PI;
	// pitch oscillation
	for (int i=0; i<3; i++)
		Body.pitchOscillation[i] = node["pitchOscillation"][i].as<real>(0.0);
	Body.pitchOscillation[0] *= M_PI / 180.0;
	Body.pitchOscillation[1] *= 2 * M_PI;

	// get the type of body and read in appropriate details
	std::string type = node["type"].as<std::string>();
	if (type == "points")
	{
		std::string fileName = node["pointsFile"].as<std::string>();
		std::string folderPath = DB["inputs"]["caseFolder"].get<std::string>();
		std::ifstream inFile((folderPath + "/" + fileName).c_str());
		inFile >> Body.numPoints;
		Body.X.resize(Body.numPoints);
		Body.Y.resize(Body.numPoints);
		for (int i=0; i<Body.numPoints; i++)
			inFile >> Body.X[i] >> Body.Y[i];
		inFile.close();
	}
	else if (type == "lineSegment")
	{
		real startX = node["segmentOptions"][0].as<real>();
		real endX = node["segmentOptions"][1].as<real>();
		real startY = node["segmentOptions"][2].as<real>();
		real endY = node["segmentOptions"][3].as<real>();
		int numPoints = node["segmentOptions"][4].as<int>();
		Body.numPoints = numPoints;
		// initialize line segment
	}
	else if (type == "circle")
	{
		real cx = node["circleOptions"][0].as<real>();
		real cy = node["circleOptions"][1].as<real>();
		real R = node["circleOptions"][2].as<real>();
		int numPoints = node["circleOptions"][3].as<int>();
		Body.numPoints = numPoints;
		// initialize circle
		Body.X.resize(numPoints);
		Body.Y.resize(numPoints);
		for (int i=0; i<Body.numPoints; i++)
		{
			Body.X[i] = cx + R*cos(i*2*M_PI/Body.numPoints);
			Body.Y[i] = cy + R*sin(i*2*M_PI/Body.numPoints);
		}
	}
	else
	{
		printf("Error: Unknown body type '%s'!\n", type.c_str());
		exit(-1);
	}
} // parseBodiesNode


/**
 * \brief Parses the \a bodies.yaml file and stores information about the immersed bodies.
 *
 * \param DB the database that contains the simulation parameters 
 */
void parseBodiesFile(std::string &bodiesFile, parameterDB &DB)
{
	if (!std::ifstream(bodiesFile.c_str()))
		return;
	printf("Parsing YAML file with bodies info ...\n");
	YAML::Node nodes = YAML::LoadFile(bodiesFile);
	body Body;
	std::vector<body> *B = DB["flow"]["bodies"].get<std::vector<body> *>();
	for (int i=0; i<nodes.size(); i++)
	{
		Body.reset();
		parseBodiesNode(nodes[i], Body, DB);
		B->push_back(Body);
	}
} // parseBodiesFile

} // End of namespace io
