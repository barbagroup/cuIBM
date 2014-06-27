/***************************************************************************//**
* \file io.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of functions related to input and output
*/

#pragma once

#include <types.h>
#include <domain.h>
#include <parameterDB.h>
#include <bodies.h>
#include <vector>

/***************************************************************************//**
* \namespace io
* \brief Contains all the functions related to I/O
*/
namespace io
{
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	std::vector<std::string> split(const std::string &s, char delim);
	void makeDirectory(const std::string s);
	
	// input
	void readInputs(int argc, char **argv, parameterDB &DB, domain &D); 

	// read in all different config files
	void parseDomainFile(std::string &domFile, domain &D);
	void parseFlowFile(std::string &flowFile, parameterDB &DB);
	void parseSimulationFile(std::string &simFile, parameterDB &DB);
	void parseBodiesFile(std::string &bodiesFile, parameterDB &DB);

	void initialiseDefaultDB(parameterDB &DB);
	void commandLineParse1(int argc, char **argv, parameterDB &DB);
	void commandLineParse2(int argc, char **argv, parameterDB &DB);

	// output
	void printSimulationInfo(parameterDB &DB, domain &D);
	void printTimingInfo(Logger &logger);
	
	void writeInfoFile(parameterDB &DB, domain &D);
	void writeGrid(std::string &caseFolder, domain &D);
	template <typename Vector>
	void writeData(std::string &caseFolder, int n, Vector &q, Vector &lambda, domain &D);//, bodies &B);
	void printDeviceMemoryUsage(char *label);
}
