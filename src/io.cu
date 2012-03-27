#include <io.h>

namespace io
{

void readInputs(int argc, char **argv, options &opts, flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
{
	flow_desc = flowDescription();
	sim_par   = simulationParameters();
	readDomainInfo(opts.domFile, dom_info);
}

void readDomainInfo(std::string domFile, domain &D)
{
	int  beg = 0, end, numCells;
	int  nsubx, nsuby;
	real h, r, a, b;
	std::ifstream f(domFile.c_str());
	
	// x-coordinates of the grid
	f	>> D.nx;					// read the total number of cells along the edge
	f	>> nsubx;				// read the number of sections which the edge is divided into
	D.x.resize(D.nx+1);
	D.dx.resize(D.nx);
	f >> D.x[beg];				// read the coordinate of the first point
	for(int i=0; i<nsubx; i++)
	{
		a	= D.x[beg];	// coordinate of the first node of the section
		f	>> b;		// coordinate of the last node of the section (which is the first node of the next section)
		f	>> numCells;		// number of cells in the current section
		f	>> r;		// stretching ratio
		end = beg + numCells;			

		if(fabs(r-1.0) < 1e-6)					// if there is no stretching (i.e. stretching factor = 1)
		{	
			h = (b - a)/numCells;						// calculate the width of each cell
			for(int j=beg; j<end; j++)
			{
				D.dx[j] = h;						// assign the width to the array
				D.x[j+1] = D.x[j] + D.dx[j];		// calculate the coordinate of the node
			}
		}
		else									// if the stretching factor is not equal to 1
		{
			h = (b-a)*(r-1)/(pow(r,numCells)-1);		// calculate the width of the first cell
			for(int j=beg; j<end; j++)
			{
				D.dx[j] = h*pow(r, j-beg);		// calculate the width of each cell in the section
				D.x[j+1] = D.x[j] + D.dx[j];			// calculate the coordinate of each node
			}
		}
		beg = end;
		D.x[end] = b;
	}
	
	// y-coordinates of the grid
	beg = 0;
	f >> D.ny >> nsuby;
	D.y.resize(D.ny+1);
	D.dy.resize(D.ny);
	f >> D.y[beg];
	for(int i=0; i<nsuby; i++)
	{
		a	= D.y[beg];
		f >> b >> numCells >> r;
		end = beg + numCells;
		if(fabs(r-1.0) < 1e-6)
		{
			h = (b - a)/numCells;
			for(int j=beg; j<beg+numCells; j++)
			{
				D.dy[j] = h;
				D.y[j+1] = D.y[j] + D.dy[j];
			}
		}
		else
		{
			h = (b-a)*(r-1)/(pow(r,numCells)-1);
			for(int j=beg; j<beg+numCells; j++)
			{
				D.dy[j] = h*pow(r, j-beg);
				D.y[j+1] = D.y[j] + D.dy[j];
			}
		}
		beg = end;
		D.y[end] = b;
	}
	D.xD.resize(D.nx+1);
	D.dxD.resize(D.nx);
	D.yD.resize(D.ny+1);
	D.dyD.resize(D.ny);
	D.xD = D.x;
	D.dxD = D.dx;
	D.yD = D.y;
	D.dyD = D.dy;
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
void writeData<vecH>(std::string folderName, int n, vecH q, vecH lambda, domain &D)//, bodies &B)
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

} // end of namespace io