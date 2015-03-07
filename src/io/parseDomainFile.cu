/***************************************************************************//**
 * \file parseDomainFile.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Parse the input file \a domain.yaml to obtain information about the 
 *        computational grid.
 */


#include <fstream>
#include <yaml-cpp/yaml.h>
#include "io.h"


/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

using std::string;

/**
 * \brief Overloads the operator >>. Gets information from the parsed domain file.
 *
 * \param node the parsed file
 * \param D instance of the class \c domain to be filled
 */
void operator >> (const YAML::Node &node, domain &D)
{
	string dir;
	real start;
	int  numCells;
	
	node["direction"] >> dir;
	node["start"] >> start;
	
	if (dir=="x")
		D.nx = 0;
	else if(dir=="y")
		D.ny = 0;
	
	const YAML::Node &subDomains = node["subDomains"];
	// first pass
	for (unsigned int i=0; i<subDomains.size(); i++)
	{
		subDomains[i]["cells"] >> numCells;
		if (dir=="x")
			D.nx += numCells;
		else if(dir=="y")
			D.ny += numCells;
	}
	
	// allocate memory
	int  beg = 0;
	if(dir=="x")
	{
		D.x.resize(D.nx+1);
		D.dx.resize(D.nx);
		D.xD.resize(D.nx+1);
		D.dxD.resize(D.nx);
		D.x[beg] = start;
	}
	if(dir=="y")
	{
		D.y.resize(D.ny+1);
		D.dy.resize(D.ny);	
		D.yD.resize(D.ny+1);
		D.dyD.resize(D.ny);
		D.y[beg] = start;
	}
	
	// second pass
	real end, stretchRatio, h;
	for (unsigned int i=0; i<subDomains.size(); i++)
	{
		subDomains[i]["end"] >> end;
		subDomains[i]["cells"] >> numCells;
		subDomains[i]["stretchRatio"] >> stretchRatio;
		
		if(fabs(stretchRatio-1.0) < 1.0e-6)
		{
			h = (end - start)/numCells;
			for(int j=beg; j<beg+numCells; j++)
			{
				if(dir=="x")
				{
					D.dx[j]  = h;
					D.x[j+1] = D.x[j] + D.dx[j];
				}
				else if(dir=="y")
				{
					D.dy[j]  = h;
					D.y[j+1] = D.y[j] + D.dy[j];
				}	
			}
		}
		else
		{
			h = (end - start)*(stretchRatio-1)/(pow(stretchRatio, numCells)-1);
			for(int j=beg; j<beg+numCells; j++)
			{
				if(dir=="x")
				{
					D.dx[j]  = h*pow(stretchRatio, j-beg);
					D.x[j+1] = D.x[j] + D.dx[j];
				}
				else if(dir=="y")
				{
					D.dy[j]  = h*pow(stretchRatio, j-beg);
					D.y[j+1] = D.y[j] + D.dy[j];
				}
			}
		}
		beg += numCells;
		start = end;
	}
	
	if(dir=="x")
	{
		D.xD  = D.x;
		D.dxD = D.dx;
	}
	else if(dir=="y")
	{
		D.yD  = D.y;
		D.dyD = D.dy;
	}
}

/**
 * \brief Parses the \a domain file and generates the computational grid.
 *
 * \param domFile the file that contains information about the computational grid
 * \param D instance of the class \c domain that will be filled with information about the computational grid
 */
void parseDomainFile(std::string &domFile, domain &D)
{
	std::ifstream fin(domFile.c_str());
	YAML::Parser  parser(fin);
	YAML::Node    doc;
	parser.GetNextDocument(doc);

	for (unsigned int i=0; i<doc.size(); i++)
		doc[i] >> D;
		
	D.xu.resize(D.nx-1);
	D.yu.resize(D.ny);
	D.xv.resize(D.nx);
	D.yv.resize(D.ny-1);
	
	int i, j;
	for(i=0; i<D.nx-1; i++)
	{
		D.xu[i] = D.x[i+1];
		D.xv[i] = (D.x[i]+D.x[i+1])/2.0;
	}
	D.xv[i] = (D.x[i]+D.x[i+1])/2.0;
	
	for(j=0; j<D.ny-1; j++)
	{
		D.yu[j] = (D.y[j]+D.y[j+1])/2.0;
		D.yv[j] = D.y[j+1];
	}
	D.yu[j] = (D.y[j]+D.y[j+1])/2.0;
}

} // end namespace io
