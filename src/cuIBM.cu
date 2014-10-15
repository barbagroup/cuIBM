#include "types.h"
#include "helpers.h"
#include "domain.h"
#include "io/io.h"
#include "solvers/NavierStokes/createSolver.h"

int main(int argc, char **argv)
{
	domain dom_info;
	parameterDB paramDB;

	io::readInputs(argc, argv, paramDB, dom_info);
	io::printSimulationInfo(paramDB, dom_info);

	/// choose the appropriate flow solver
	NavierStokesSolver<device_memory> *solver = createSolver<device_memory>(paramDB, dom_info);
	solver->initialise();
	io::printDeviceMemoryUsage("Initialisation complete");
	
	io::writeInfoFile(paramDB, dom_info);
	
	while (!solver->finished())
	{
		solver->stepTime();
		solver->writeData();
	}
	solver->shutDown();
}
