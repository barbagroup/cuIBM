/**
 * \file main.cu
 * \brief Unit-test to check the spatial order of convergence of the convective terms.
 */


#include "utilities/types.h"
#include "utilities/helpers.h"
#include "utilities/domain.h"
#include "io/io.h"
#include "solvers//NavierStokesSolver.h"
#include "convectionTerm.h"


int main(int argc, char **argv)
{
	domain dom_info;
	parameterDB paramDB;

	io::readInputs(argc, argv, paramDB, dom_info);

	/// create an object to test convection kernel
	convectionTerm<device_memory> *solver = new convectionTerm<device_memory>(paramDB, dom_info);
	
	solver->initialise();
	
	solver->writeData();

	solver->shutDown();
}
