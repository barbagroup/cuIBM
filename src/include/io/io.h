#pragma once

#include <types.h>
#include <domain.h>
#include <parameterDB.h>
#include <bodies.h>
#include <vector>

namespace io
{
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
	
	void writeGrid(std::string &folderName, domain &D);
	template <typename Vector>
	void writeData(std::string &folderName, int n, Vector &q, Vector &lambda, domain &D);//, bodies &B);
	void writeForce(std::string &folderName, real t, real forceX, real forceY);
	void writeIterations(std::string &folderName, size_t n, size_t iter1, size_t iter2);
//  void writeDataDevice(std::string &s, int n, vecD &q, vecD &lambda, domain &D);//, bodies &B);
}
