/***************************************************************************//**
* \file io.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the functions of the namespace \c io
*/

#pragma once

#include <vector>
#include <types.h>
#include <domain.h>
#include <parameterDB.h>
#include <bodies.h>

/**
* \namespace io
* \brief Contain functions related to I/O tasks
*/
namespace io
{
	/**
	* \brief Split a string given a delimiter
	*/
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	
	/**
	* \brief Split a string given a delimiter
	*/
	std::vector<std::string> split(const std::string &s, char delim);

	/**
	* \brief Make a directory
	*/
	void makeDirectory(const std::string s);
	
	/**
	* \brief Read inputs by parsing the command-line and simulation files
	*/
	void readInputs(int argc, char **argv, parameterDB &DB, domain &D); 

	/**
	* \brief Parse the \a domain file using YAML
	*/
	void parseDomainFile(std::string &domFile, domain &D);
	
	/**
	* \brief Parse the \a flow file using YAML
	*/
	void parseFlowFile(std::string &flowFile, parameterDB &DB);
	
	/**
	* \brief Parse the \a simulation file using YAML
	*/
	void parseSimulationFile(std::string &simFile, parameterDB &DB);
	
	/**
	* \brief Parse the \a bodies file using YAML
	*/
	void parseBodiesFile(std::string &bodiesFile, parameterDB &DB);

	/**
	* \brief Initialize the database with default values
	*/
	void initialiseDefaultDB(parameterDB &DB);

	/**
	* \brief Parse the command-line to get the case folder name and the device number
	*/
	void commandLineParse1(int argc, char **argv, parameterDB &DB);

	/**
	* \brief Overwrite values in the database with additional arguments of the command-line
	*/
	void commandLineParse2(int argc, char **argv, parameterDB &DB);

	/**
	* \brief Print the parameters of the simulation
	*/
	void printSimulationInfo(parameterDB &DB, domain &D);

	/**
	*\brief Print the time spent on certain tasks
	*/
	void printTimingInfo(Logger &logger);
	
	/**
	* \brief Write information about the run into a file
	*/
	void writeInfoFile(parameterDB &DB, domain &D);

	/**
	* \brief Write the mesh points into a file named \a grid
	*/
	void writeGrid(std::string &caseFolder, domain &D);

	/**
	* \brief Write numerical data at a given iteration number
	*/
	template <typename Vector>
	void writeData(std::string &caseFolder, int n, Vector &q, Vector &lambda, domain &D);//, bodies &B);
	
	/**
	* \brief Print the memory usage on the device
	*/
	void printDeviceMemoryUsage(char *label);
}
