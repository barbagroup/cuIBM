/***************************************************************************//**
* \mainpage cuIBM
*
*		A GPU-based Immersed Boundary Method
*
* \author Krishnan, A. (anush@bu.edu)
*/

/***************************************************************************//**
* \file cuIBM.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Main source-file of <B>cuIBM</B>
*/

#include <types.h>
#include <helpers.h>
#include <domain.h>
#include <io/io.h>
#include <solvers/NavierStokes/createSolver.h>

int main(int argc, char **argv)
{
	// initializes the computational domain
	domain dom_info;

	// initializes the parameters of the simulation
	parameterDB paramDB;

	// reads input files
	io::readInputs(argc, argv, paramDB, dom_info);
	
	io::printSimulationInfo(paramDB, dom_info);

	// chooses the appropriate flow solver
	NavierStokesSolver<device_memory> *solver = createSolver<device_memory>(paramDB, dom_info);
	solver->initialise();
	io::printDeviceMemoryUsage("Initialisation complete");
	
	io::writeInfoFile(paramDB, dom_info);
	
	// time-step loop
	while (!solver->finished())
	{
		solver->stepTime();
		solver->writeData();
	}

	solver->shutDown();
}
