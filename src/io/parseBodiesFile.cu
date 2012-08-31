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
#include <yaml-cpp/yaml.h>
#include <body.h>
#include <vector>
#include <fstream>

namespace io
{

using std::string;

void operator >> (const YAML::Node &node, body &Body)
{
	// read in center
	for (int i=0; i<2; i++)
		node["centerRotation"][i] >> Body.X0[i];
	
	// moving flags
	for (int i=0; i<2; i++)
		node["moving"][i] >> Body.moving[i];
	
	// velocity
	for (int i=0; i<2; i++)
		node["velocity"][i] >> Body.velocity[i];

	// alpha
	node["alpha"] >> Body.Theta0;
	Body.Theta0 = Body.Theta0 * M_PI/180.0;
	Body.Theta  = Body.Theta0;

	// omega
	node["omega"] >> Body.omega;
	
	// oscillation in X, Y, pitch
	for (int i=0; i<3; i++)
		node["xOscillation"][i] >> Body.xOscillation[i];
	for (int i=0; i<3; i++)
		node["yOscillation"][i] >> Body.yOscillation[i];
	for (int i=0; i<3; i++)
		node["pitchOscillation"][i] >> Body.pitchOscillation[i];

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
	printf("numPoints: %d\n",Body.numPoints);
}

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

