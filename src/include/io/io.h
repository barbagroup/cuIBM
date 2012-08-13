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

/**
* @file io.h
* @brief Functions for input and output.
*/

#pragma once

#include <types.h>
#include <domain.h>
#include <parameterDB.h>
#include <bodies.h>
#include <vector>

/**
* @namespace io
* @brief     Contains all the functions related to I/O
*/
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
	
	void writeInfoFile(parameterDB &DB, domain &D);
	void writeGrid(std::string &folderName, domain &D);
	template <typename Vector>
	void writeData(std::string &folderName, int n, Vector &q, Vector &lambda, domain &D);//, bodies &B);
}
