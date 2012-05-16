#include <types.h>
#include <domain.h>
#include <io/io.h>
#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>

int main(int argc, char **argv)
{
	std::cout << std::endl;
	
	domain dom_info;
	parameterDB paramDB;

	cudaSetDevice(0);

	io::readInputs(argc, argv, paramDB, dom_info);
	io::printSimulationInfo(paramDB, dom_info);

	NavierStokesSolver<device_memory> *solver = NavierStokesSolver<device_memory>::createSolver(paramDB, dom_info);
	solver->initialise();

	while (!solver->finished())
	{
		solver->stepTime();
		solver->writeData();
	}
}
