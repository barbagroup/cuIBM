/**
* \file
* \brief
*/
#pragma once

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
	void printSimulationInfo(options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom);
}