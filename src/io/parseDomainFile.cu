#include <io/io.h>
#include <yaml-cpp/yaml.h>
#include <fstream>

namespace io
{

using std::string;

void operator >> (const YAML::Node &node, domain &d)
{
  string dir;
  node["direction"] >> dir;
  double start;
  node["start"] >> start;
  //if (dir == "x") d.startX = start;
  //else            d.startY = start;

  const YAML::Node &subDomains = node["subDomains"];
  // first pass
  for (unsigned int i=0; i<subDomains.size(); i++)
  {
    int nPoints;
    subDomains[i]["points"] >> nPoints;
    if (dir == "x") d.nx += nPoints;
    else            d.ny += nPoints;
  }
  // allocate memory <ANUSH>

  // second pass
  for (unsigned int i=0; i<subDomains.size(); i++)
  {
    double end, stretchRatio;
    int nPoints;
    subDomains[i]["end"] >> end;
    subDomains[i]["points"] >> nPoints;
    subDomains[i]["stretchRatio"] >> stretchRatio;

    /* <ANUSH>
    if (dir == "x") do_something_with_X
    else            do_something_with_Y
    */
  }
}

void parseDomainFile(std::string &domFile, domain &D)
{
  std::ifstream fin(domFile.c_str());
  YAML::Parser parser(fin);
  YAML::Node doc;
  parser.GetNextDocument(doc);

  for (unsigned int i=0; i<doc.size(); i++)
    doc[i] >> D;
}

} // end namespace io
