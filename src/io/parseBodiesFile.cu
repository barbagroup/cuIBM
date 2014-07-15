/***************************************************************************//**
* \file parseBodiesFile.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Parse the \a bodies file to get information about immersed bodies
*/

#include <io/io.h>
#include <yaml-cpp/yaml.h>
#include <body.h>
#include <vector>
#include <fstream>

/**
* \namespace io
* \brief Contain functions related to I/O tasks
*/
namespace io
{

using std::string;

/**
* \brief Set up onformation about the immersed bodies
*
* \param node the parsed file
* \param instance of the class \c body to be filled
*/
void operator >> (const YAML::Node &node, body &Body)
{
	// read in center of rotation
	for (int i=0; i<2; i++)
		node["centerRotation"][i] >> Body.X0[i];
	
	// initial configuration
	for(int i=0; i<2; i++)
	{
		node["initialOffset"][i] >> Body.Xc0[i];
		Body.Xc[i] = Body.Xc0[i];
	}
	
	// initial angle of attack
	node["angleOfAttack"] >> Body.Theta0;
	Body.Theta0 *= M_PI/180.0;
	Body.Theta = Body.Theta0;
	
	// moving flags
	for (int i=0; i<2; i++)
		node["moving"][i] >> Body.moving[i];
	
	// velocity
	for (int i=0; i<2; i++)
		node["velocity"][i] >> Body.velocity[i];

	// omega
	node["omega"] >> Body.omega;
	
	// oscillation in X
	for (int i=0; i<3; i++)
		node["xOscillation"][i] >> Body.xOscillation[i];
	Body.xOscillation[1] *= 2*M_PI;
	
	// oscillation in Y
	for (int i=0; i<3; i++)
		node["yOscillation"][i] >> Body.yOscillation[i];
	Body.yOscillation[1] *= 2*M_PI;
	
	// pitch oscillation
	for (int i=0; i<3; i++)
		node["pitchOscillation"][i] >> Body.pitchOscillation[i];
	Body.pitchOscillation[0] *= M_PI/180.0;
	Body.pitchOscillation[1] *= 2*M_PI;

	// get the type of Body and read in appropriate details
	string type;
	node["type"] >> type;
	if (type == "points")
	{
		string fname;
		node["pointsFile"] >> fname;
		fname = "bodyFiles/" + fname;
		std::cout << fname << std::endl;
		// initialise points
		std::ifstream file(fname.c_str());
		file >> Body.numPoints;
		Body.X.resize(Body.numPoints);
		Body.Y.resize(Body.numPoints);
		for(int i=0; i<Body.numPoints; i++)
		{
			file >> Body.X[i] >> Body.Y[i];
		}
		file.close();
	}
	else if (type == "lineSegment")
	{
		real startX, startY, endX, endY;
		int numPoints;
		node["segmentOptions"][0] >> startX;
		node["segmentOptions"][1] >> endX;
		node["segmentOptions"][2] >> startY;
		node["segmentOptions"][3] >> endY;
		node["segmentOptions"][4] >> numPoints;
		Body.numPoints = numPoints;
		// initialise line segment
	}
	else if (type == "circle")
	{
		real cx, cy, R;
		int numPoints;
		node["circleOptions"][0] >> cx;
		node["circleOptions"][1] >> cy;
		node["circleOptions"][2] >> R;
		node["circleOptions"][3] >> numPoints;
		Body.numPoints = numPoints;
		// initialise circle
		Body.X.resize(numPoints);
		Body.Y.resize(numPoints);
		for(int i=0; i<Body.numPoints; i++)
		{
			Body.X[i] = cx + R*cos(i*2*M_PI/Body.numPoints);
			Body.Y[i] = cy + R*sin(i*2*M_PI/Body.numPoints);
		}
	}
	else
		printf("[E]: unknown Body type\n");
	printf("\nnumPoints: %d",Body.numPoints);
}

/**
* \brief Parse the \a bodies file using YAML
*
* \param bodiesFile the file that contains information about the immersed bodies
* \param DB the database that will be filled with information about the immersed bodies
*/
void parseBodiesFile(std::string &bodiesFile, parameterDB &DB)
{
	std::ifstream fin(bodiesFile.c_str());
	YAML::Parser  parser(fin);
	YAML::Node    doc;
	parser.GetNextDocument(doc);
	body Body;
	std::vector<body> *B = DB["flow"]["bodies"].get<std::vector<body> *>();
	
	for (int i=0; i<doc.size(); i++)
	{  
		doc[i] >> Body;
		B->push_back(Body);
	}
}

} // end namespace io

