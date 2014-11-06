/***************************************************************************//**
 * \file io.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the functions of the namespace \c io.
 */


#pragma once

#include <vector>
#include <types.h>
#include <domain.h>
#include <parameterDB.h>
#include <bodies.h>


/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{
	// split a string given a delimiter
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	
	// split a string given a delimiter
	std::vector<std::string> split(const std::string &s, char delim);

	// create a directory
	void makeDirectory(const std::string s);
	
	// read data inputs from the command-line and the simulation files
	void readInputs(int argc, char **argv, parameterDB &DB, domain &D); 

    // parse the \a domain file and generate the computational grid
	void parseDomainFile(std::string &domFile, domain &D);
	
	// parse the \a flow file and store the parameters in the database
	void parseFlowFile(std::string &flowFile, parameterDB &DB);
	
	// parse the \a simulation file and store the parameters in the database
	void parseSimulationFile(std::string &simFile, parameterDB &DB);
	
	// parse the \a bodies file and store information about the immersed bodies
	void parseBodiesFile(std::string &bodiesFile, parameterDB &DB);

	// initialize the database with default values
	void initialiseDefaultDB(parameterDB &DB);

	// parse command-line to get simulation directory and device number
	void commandLineParse1(int argc, char **argv, parameterDB &DB);

	// overwrite parameters with additional arguments of the command-line
	void commandLineParse2(int argc, char **argv, parameterDB &DB);

	// print the parameters of the simulation
	void printSimulationInfo(parameterDB &DB, domain &D);

	// print the time spent to execute tasks
	void printTimingInfo(Logger &logger);
	
	// write information about the run into the file run.info
	void writeInfoFile(parameterDB &DB, domain &D);

	// write grid-points coordinates into the file grid
	void writeGrid(std::string &caseFolder, domain &D);

	// write numerical data at a given time-step
	template <typename Vector>
	void writeData(std::string &caseFolder, int n, Vector &q, Vector &lambda, domain &D);//, bodies &B);
	
	// print device memory usage
	void printDeviceMemoryUsage(char *label);
}
