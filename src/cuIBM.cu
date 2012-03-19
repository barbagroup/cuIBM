#include <types.h>
#include <options.h>
#include <flowDescription.h>
#include <simulationParameters.h>
#include <domain.h>
#include <io.h>
#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>

int main(int argc, char **argv)
{
	options opts;
	flowDescription flow_desc;
	simulationParameters sim_par;
	domain dom_info;
	
	io::readInputs(argc, argv, opts, flow_desc, sim_par, dom_info);
	io::printSimulationInfo(opts, flow_desc, sim_par, dom_info);
	
	NavierStokesSolver<matrix,vector> *solver = NavierStokesSolver<matrix,vector>::createSolver(flow_desc, sim_par, dom_info);
	solver->initialise();
	
	while (!solver->finished())
	{	
		solver->stepTime();
		solver->writeData();
	}
}