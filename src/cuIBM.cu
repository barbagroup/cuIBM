#include <types.h>
#include <options.h>
#include <flowDescription.h>
#include <simulationParameters.h>
#include <domain.h>
#include <io/io.h>
#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>

int main(int argc, char **argv)
{
	domain dom_info;
  parameterDB paramDB;
	
	cudaSetDevice(4);
	
	io::readInputs(argc, argv, paramDB, dom_info);
	// io::printSimulationInfo(opts, flow_desc, paramDB); // sim_par, dom_info);
	
	NavierStokesSolver<device_memory> *solver = NavierStokesSolver<device_memory>::createSolver(paramDB, dom_info);
	solver->initialise();

	while (!solver->finished())
	{	 
		solver->stepTime();
		solver->writeData();
	}
}
