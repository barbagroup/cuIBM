/**
* \file
* \brief
*/
#pragma once

#include <sys/stat.h>

#include <types.h>
#include <flowDescription.h>
#include <simulationParameters.h>
#include <domain.h>
#include <integrationScheme.h>
#include <options.h>

#include <cusp/print.h>

namespace io
{
	void readInputs(int argc, char **argv, options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom);
	void readDomainInfo(std::string domFile, domain &D);
	void printSimulationInfo(options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom);
	void writeGrid(std::string folderName, domain &D);
	template <typename Vector>
	void writeData(std::string s, int n, Vector q, Vector lambda, domain &D);//, bodies &B);
}