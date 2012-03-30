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
  for (int i=0; i<2; i++) node["centerRotation"][i] >> Body.X0[i];
  // moving flags
  for (int i=0; i<2; i++) node["moving"][i] >> Body.moving[i];
  // velocity
  for (int i=0; i<2; i++) node["velocity"][i] >> Body.velocity[i];
  // omega
  node["omega"] >> Body.Theta0;
  // oscillation in X, Y, pitch
  for (int i=0; i<3; i++) node["xOscillation"][i] >> Body.xOscillation[i];
  for (int i=0; i<3; i++) node["yOscillation"][i] >> Body.yOscillation[i];
  for (int i=0; i<3; i++) node["pitchOscillation"][i] >> Body.pitchOscillation[i];

  // get the type of Body and read in appropriate details
  string type;
  node["type"] >> type;
  if (type == "points")
  {
    string fname;
    node["pointsFile"] >> fname;
    // initialisePoints(fname.c_str(),Body);
  }
  else if (type == "line_segment")
  {
    double startX, startY, endX, endY;
    int nPoints;
    node["segmentOptions"][0] >> startX;
    node["segmentOptions"][1] >> endX;
    node["segmentOptions"][2] >> startY;
    node["segmentOptions"][3] >> endY;
    node["segmentOptions"][4] >> nPoints;
    Body.numPoints = nPoints;
    // initialise line segment
  }
  else if (type == "circle")
  {
    double cx, cy, R;
    int nPoints;
    node["circleOptions"][0] >> cx;
    node["circleOptions"][1] >> cy;
    node["circleOptions"][2] >> R;
    node["circleOptions"][3] >> nPoints;
    Body.numPoints = nPoints;
    printf("numPoints: %d\n",nPoints);
    // initialise circle
  }
  else
    printf("[E]: unknown Body type\n");
}

void parseBodiesFile(std::string &bodiesFile, parameterDB &DB)
{
  std::ifstream fin(bodiesFile.c_str());
  YAML::Parser parser(fin);
  YAML::Node doc;
  parser.GetNextDocument(doc);

  std::vector<body> *B = DB["flow"]["bodies"].get<std::vector<body> *>();
  for (int i=0; i<doc.size(); i++)
  {
    body Body;
    doc >> Body;
    B->push_back(Body);
  }
}

} // end namespace io

