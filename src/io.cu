#include <io.h>

namespace io
{
	void readInputs(int argc, char **argv, options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
	{
		flow_desc = flowDescription();
		sim_par   = simulationParameters();
		dom_info  = domain();
	}
	
	void printSimulationInfo(options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
	{
		std::cout << std::endl;
		std::cout << "Simulation parameters" << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << "dt = " << sim_par.dt         << std::endl;
		std::cout << "startStep = " << sim_par.startStep << std::endl;
		std::cout << "nt = "    << sim_par.nt    << std::endl;
		std::cout << "nsave = " << sim_par.nsave << std::endl;
		std::cout << "T (start) = " << sim_par.startStep*sim_par.dt << std::endl;
		std::cout << "T (end)   = " << (sim_par.startStep + sim_par.nt)*sim_par.dt << std::endl << std::endl;
	
		std::cout << "Flow parameters" << std::endl;
		std::cout << "---------------" << std::endl;
		std::cout << "nu = " << flow_desc.nu << std::endl << std::endl;
		
		std::cout << "Domain" << std::endl;
		std::cout << "------" << std::endl;
		std::cout << dom_info.nx << " x " << dom_info.ny << std::endl << std::endl;
	}
}
