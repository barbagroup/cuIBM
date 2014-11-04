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
	/**
	 * \brief Splits a string given a delimiter.
	 */
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	
	/**
	 * \brief Splits a string given a delimiter.
	 */
	std::vector<std::string> split(const std::string &s, char delim);

	/**
	 * \brief Creates a directory.
	 */
	void makeDirectory(const std::string s);
	
	/**
	 * \brief Reads data inputs by parsing the command-line and simulation files.
	 */
	void readInputs(int argc, char **argv, parameterDB &DB, domain &D); 

	/**
	 * \brief Parses the \a domain file using YAML.
.	 */
	void parseDomainFile(std::string &domFile, domain &D);
	
	/**
	 * \brief Parses the \a flow file using YAML.
	 */
	void parseFlowFile(std::string &flowFile, parameterDB &DB);
	
	/**
	 * \brief Parses the \a simulation file using YAML.
	 */
	void parseSimulationFile(std::string &simFile, parameterDB &DB);
	
	/**
	 * \brief Parses the \a bodies file using YAML.
	 */
	void parseBodiesFile(std::string &bodiesFile, parameterDB &DB);

	/**
	 * \brief Initializes the database with default values.
	 */
	void initialiseDefaultDB(parameterDB &DB);

	/**
	 * \brief Parses the command-line to get the case folder name 
	 *        and the device number.
	 */
	void commandLineParse1(int argc, char **argv, parameterDB &DB);

	/**
	 * \brief Overwrites values in the database 
	 *        with additional arguments of the command-line.
	 */
	void commandLineParse2(int argc, char **argv, parameterDB &DB);

	/**
	 * \brief Prints the parameters of the simulation.
	 */
	void printSimulationInfo(parameterDB &DB, domain &D);

	/**
	 * \brief Prints the time spent on certain tasks.
	 */
	void printTimingInfo(Logger &logger);
	
	/**
	 * \brief Writes information about the run into a file.
	 */
	void writeInfoFile(parameterDB &DB, domain &D);

	/**
	 * \brief Writes the mesh points into a file named \a grid.
	 */
	void writeGrid(std::string &caseFolder, domain &D);

	/**
	 * \brief Writes numerical data at a given iteration number.
	 */
	template <typename Vector>
	void writeData(std::string &caseFolder, int n, Vector &q, Vector &lambda, domain &D);//, bodies &B);
	
	/**
	 * \brief Prints the memory usage on the device.
	 */
	void printDeviceMemoryUsage(char *label);
}
