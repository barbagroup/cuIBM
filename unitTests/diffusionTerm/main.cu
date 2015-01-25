#include <types.h>
#include <helpers.h>
#include <domain.h>
#include <io/io.h>
#include <solvers/NavierStokes/NavierStokesSolver.h>
#include "diffusionTerm.h"

int main(int argc, char **argv)
{
	domain dom_info;
	parameterDB paramDB;

	io::readInputs(argc, argv, paramDB, dom_info);

	/// create an object to test convection kernel
	diffusionTerm<device_memory> *solver = new diffusionTerm<device_memory>(paramDB, dom_info);
	
	solver->initialise();
	
	io::writeInfoFile(paramDB, dom_info);
	
	solver->writeData();

	solver->shutDown();
}
