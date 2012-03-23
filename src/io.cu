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
	
	/*static int createDirectory(const char *path, mode_t mode)
	{
		stat   st;
		int    status = 0;

		if (stat(path, &st) != 0)
		{
			// Directory does not exist
			if (mkdir(path, mode) != 0)
				status = -1;
		}
		else if (!S_ISDIR(st.st_mode))
		{
			errno = ENOTDIR;
			status = -1;
		}

	    return(status);
	}*/
	
	/*template <typename Vector>
	void writeArray()
	{
		
	}*/
	
	void writeGrid(std::string folderName, domain &D)
	{
		std::stringstream out;
		out << folderName << "/grid";
		std::ofstream f(out.str().c_str());
		
		f << D.nx << std::endl; 
		for(int i=0; i<D.nx+1; i++)
			f << D.x[i] << std::endl;
		f << std::endl;
		f << D.ny << std::endl;
		for(int j=0; j<D.ny+1; j++)
			f << D.y[j] << std::endl;
			
		f.close();
	}

	template <>
	void writeData<vecH>(std::string folderName, int n, vecH q, vecH lambda, domain &D//, bodies &B)
	{
		std::string path;
		std::stringstream out;
		int numUV = D.nx*(D.ny-1) + D.ny*(D.nx-1),
		    numP  = D.nx*D.ny,
		    numB  = 0;

		out << folderName << '/' << std::setfill('0') << std::setw(7) << n;
		path = out.str();

		//createDirectory(path.c_str(), S_IRWXO);
		mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
		out.str("");
		out << path << "/q";
		//out << "q";
		std::ofstream file(out.str().c_str());
		file << numUV << std::endl;
		for(int i=0; i<numUV; i++)
			file << q[i] << std::endl;
		file.close();

		/*
		out.str("");
		out << path << "/bodies";
		file.open(out.str().c_str());
		file << '#' << std::setw(19) << "x-coordinate" << std::setw(20) << "y-coordinate" << std::endl;
		for (int i=0; i < B.no_of_bodies; i++)
		{
			for (int j=0; j < B.Body[i].no_of_points; j++)
			{
				file << std::setw(20) << B.Body[i].bp[j].x << std::setw(20) << B.Body[i].bp[j].y << std::endl;
				nb++;
			}
			file << std::endl;
		}
		file.close();*/

		out.str("");
		out << path << "/lambda";
		//out << "lambda";
		file.open(out.str().c_str());
		file << numP+2*numB << std::endl;
		for(int i=0; i<numP+2*numB; i++)
			file << lambda[i] << std::endl;
		file.close();
	
		std::cout << "Data saved to folder " << path << std::endl;
	}
	
	template <>
	void writeData<vecD>(std::string folderName, int n, vecD q, vecD lambda, domain &D)//, bodies &B)
	{
		vecH qH = q,
		     lambdaH = lambda;
		writeData(folderName, n, qH, lambdaH, D);
	}
}