/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

#include <io/io.h>
#include <types.h>
#include <sys/stat.h>
#include <boundaryCondition.h>

using std::string;

/// convert string to real number or integer
template <typename T>
T toNumber(string str)
{
     T num;
     std::stringstream ss(str); //turn the string into a stream
     ss >> num; //convert
     return num;
}

namespace io
{

//##############################################################################
//                                 INPUT
//##############################################################################
	
void readInputs(int argc, char **argv, parameterDB &DB, domain &D)
{
	// get a default database
	initialiseDefaultDB(DB);
	
	// first pass of command line arguments
	commandLineParse1(argc, argv, DB);
	
	// read the simulation file
	string fname = DB["inputs"]["simulationFile"].get<string>();
	parseSimulationFile(fname, DB);
	
	// read the flow file
	fname = DB["inputs"]["flowFile"].get<string>();
	parseFlowFile(fname, DB);

	// read the domain file
	fname = DB["inputs"]["domainFile"].get<string>();
	parseDomainFile(fname, D);
	
	// read the body file
	fname = DB["inputs"]["bodyFile"].get<string>();
	parseBodiesFile(fname, DB);
	
	// second pass of command line -- overwrite values in DB
	commandLineParse2(argc, argv, DB);
}

void initialiseDefaultDB(parameterDB &DB)
{
	DB["inputs"] = componentParameter();
	DB["flow"] = componentParameter();
	DB["simulation"] = componentParameter();
	DB["velocitySolve"] = componentParameter();
	DB["PoissonSolve"] = componentParameter();

	// default input files
	string inputs = "inputs";
	DB[inputs]["flowFile"].set<string>("flows/openFlow.yaml");
	DB[inputs]["simulationFile"].set<string>("simParams/openFlow.yaml");
	DB[inputs]["bodyFile"].set<string>("bodies/cylinder.yaml");
	DB[inputs]["domainFile"].set<string>("domains/openFlow.yaml");
	DB[inputs]["folderName"].set<string>("new");
	DB[inputs]["deviceNumber"].set<int>(0);

	// flow parameters
	string flow = "flow";
	DB[flow]["nu"].set<real>(0.01);
	DB[flow]["uInitial"].set<real>(1.0);
	DB[flow]["vInitial"].set<real>(0.0);
	DB[flow]["numBodies"].set<int>(0);
	std::vector<body> *bodyVec = new std::vector<body>;
	DB[flow]["bodies"].set<std::vector<body> *>(bodyVec);

	boundaryCondition **bc = new boundaryCondition*[4];
	for (int i=0; i<4; i++)
		bc[i] = new boundaryCondition[2];
	DB[flow]["boundaryConditions"].set<boundaryCondition **>(bc);

	// simulation parameters
	string sim = "simulation";
	DB[sim]["dt"].set<real>(0.02);
	DB[sim]["nt"].set<int>(100);
	DB[sim]["nsave"].set<int>(100);
	DB[sim]["restart"].set<bool>(false);
	DB[sim]["startStep"].set<bool>(0);
	DB[sim]["convTimeScheme"].set<timeScheme>(EULER_EXPLICIT);
	DB[sim]["diffTimeScheme"].set<timeScheme>(EULER_IMPLICIT);
	DB[sim]["ibmScheme"].set<ibmScheme>(TAIRA_COLONIUS);

	// velocity solver
	string solver = "velocitySolve";
	DB[solver]["solver"].set<string>("CG");
	DB[solver]["preconditioner"].set<preconditionerType>(DIAGONAL);
	DB[solver]["tolerance"].set<real>(1e-5);
	DB[solver]["maxIterations"].set<int>(10000);

	// Poisson solver
	solver = "PoissonSolve";
	DB[solver]["solver"].set<string>("CG");
	DB[solver]["preconditioner"].set<preconditionerType>(DIAGONAL);
	DB[solver]["tolerance"].set<real>(1e-5);
	DB[solver]["maxIterations"].set<int>(20000);
}

// first pass - only get the files to continue
void commandLineParse1(int argc, char **argv, parameterDB &DB)
{
	for (int i=1; i<argc; i++)
	{
		if (strcmp(argv[i],"-flowFile")==0)
		{
			i++;
			DB["inputs"]["flowFile"].set<string>(string(argv[i]));
		}
		else if (strcmp(argv[i],"-simulationFile")==0)
		{
			i++;
			DB["inputs"]["simulationFile"].set<string>(string(argv[i]));
		}
		else if (strcmp(argv[i],"-bodyFile")==0)
		{
			i++;
			DB["inputs"]["bodyFile"].set<string>(string(argv[i]));
		}
		else if (strcmp(argv[i],"-domainFile")==0)
		{
			i++;
			DB["inputs"]["domainFile"].set<string>(string(argv[i]));
		}
		else if (strcmp(argv[i],"-folderName")==0)
		{
			i++;
			DB["inputs"]["folderName"].set<string>(string(argv[i]));
		}
		else if (strcmp(argv[i],"-deviceNumber")==0)
		{
			i++;
			int devNum = toNumber<int>(string(argv[i]));
			DB["inputs"]["deviceNumber"].set<int>(devNum);
			cudaSetDevice(devNum);
		}
	}
}

// overwrite values in the DB from command line
void commandLineParse2(int argc, char **argv, parameterDB &DB)
{
	for (int i=1; i<argc; i++)
	{
		// ignore these - already parsed in pass 1
		/*
		if (strcmp(argv[i],"-flowFile")==0 ||
		strcmp(argv[i],"-simulationFile")==0 ||
		strcmp(argv[i],"-bodyFile")==0 ||
		strcmp(argv[i],"-domainFile")==0 ||
		strcmp(argv[i],"-folderName")==0 ||
		strcmp(argv[i],"-deviceNumber")==0)
			continue;
		*/
			
		// Add code here
		// kinematic viscosity
		if ( strcmp(argv[i],"-nu")==0 )
		{
			i++;
			DB["flow"]["nu"].set<real>(toNumber<real>(string(argv[i])));
		}
		// perturbation in the x-velocity
		if ( strcmp(argv[i],"-uPerturb")==0 )
		{
			i++;
			DB["flow"]["uPerturb"].set<real>(toNumber<real>(string(argv[i])));
		}
		// scale the CV with respect to the body
		if ( strcmp(argv[i],"-scaleCV")==0 )
		{
			i++;
			DB["simulation"]["scaleCV"].set<real>(toNumber<real>(string(argv[i])));
		}
		// frequency of saving the data
		if ( strcmp(argv[i],"-nsave")==0 )
		{
			i++;
			DB["simulation"]["nsave"].set<int>(toNumber<int>(string(argv[i])));
		}
		// total number of time steps
		if ( strcmp(argv[i],"-nt")==0 )
		{
			i++;
			DB["simulation"]["nt"].set<int>(toNumber<int>(string(argv[i])));
		}
		// IBM Scheme
		if ( strcmp(argv[i],"-ibmScheme")==0 )
		{
			i++;
			if ( strcmp(argv[i],"NavierStokes")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(NAVIER_STOKES);
			else
			if ( strcmp(argv[i],"TairaColonius")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(TAIRA_COLONIUS);
			else 
			if ( strcmp(argv[i],"FadlunEtAl")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(FADLUN_ET_AL);
		}
	}
}

//##############################################################################
//                                OUTPUT
//##############################################################################

string stringFromPreconditionerType(preconditionerType s)
{
  if (s == NONE)
    return "None";
  else if (s == DIAGONAL)
    return "Diagonal";
  else if (s == SMOOTHED_AGGREGATION)
    return "Smoothed Aggregation";
  else
    return "Unrecognised preconditioner";
}

// output
void printSimulationInfo(parameterDB &DB, domain &D)
{
	real dt = DB["simulation"]["dt"].get<real>(),
	     scaleCV = DB["simulation"]["scaleCV"].get<real>();
	int  nt = DB["simulation"]["nt"].get<int>(),
	     nsave = DB["simulation"]["nsave"].get<int>(),
	     startStep = DB["simulation"]["startStep"].get<int>();

    std::cout << '\n';
	
	std::cout << "\nFlow parameters" << '\n';
	std::cout << "---------------" << '\n';
	std::cout << "nu = " << DB["flow"]["nu"].get<real>() << '\n';

	std::cout << "\nDomain" << '\n';
	std::cout << "------" << '\n';
	std::cout << D.nx << " x " << D.ny << '\n';
	
	std::cout << "\nSimulation parameters" << '\n';
	std::cout << "---------------------" << '\n';
	std::cout << "dt = " << dt << '\n';
	std::cout << "scaleCV = " << scaleCV << '\n';
	std::cout << "startStep = " << startStep << '\n';
	std::cout << "nt = "    << nt << '\n';
	std::cout << "nsave = " << nsave << '\n';
	
	std::cout << "\nVelocity Solve" << '\n';
	std::cout << "--------------" << '\n';
	std::cout << "Solver = " << DB["velocitySolve"]["solver"].get<string>() << '\n';
	std::cout << "Preconditioner = " << stringFromPreconditionerType(DB["velocitySolve"]["preconditioner"].get<preconditionerType>()) << '\n';
	std::cout << "Tolerance = " << DB["velocitySolve"]["tolerance"].get<real>() << '\n';
	
	std::cout << "\nPoisson Solve" << '\n';
	std::cout << "-------------" << '\n';
	std::cout << "Solver = " << DB["PoissonSolve"]["solver"].get<string>() << '\n';
	std::cout << "Preconditioner = " << stringFromPreconditionerType(DB["PoissonSolve"]["preconditioner"].get<preconditionerType>()) << '\n';
	std::cout << "Tolerance = " << DB["PoissonSolve"]["tolerance"].get<real>() << '\n';
	
	std::cout << "\nOutput parameters" << '\n';
	std::cout << "-----------------" << '\n';
	std::cout << "Output folder = " << DB["inputs"]["folderName"].get<string>() << '\n';
	std::cout << "nsave = " << DB["simulation"]["nsave"].get<int>() << '\n';
	
	cudaDeviceProp deviceProp;
	int gpu = DB["inputs"]["deviceNumber"].get<int>();
	cudaGetDeviceProperties(&deviceProp, gpu);
	std::cout << "\nDevice Properties" << '\n';
	std::cout << "-----------------" << '\n';
	std::cout << "Name = " << deviceProp.name << '\n';
	std::cout << "Number = " << gpu << '\n';
	std::string ecc = deviceProp.ECCEnabled ? "yes" : "no";
	std::cout << "Compute capability = " << deviceProp.major << "." << deviceProp.minor << '\n';
	std::cout << "ECC Enabled = " << ecc << std::endl;
}

void printTimingInfo(Logger &logger)
{
	//logger.writeLegend();
	logger.printAllTime();
	std::cout << std::endl;
}

void writeInfoFile(parameterDB &DB, domain &D)
{
	std::string   foldername = DB["inputs"]["folderName"].get<string>();
	std::ofstream infofile((foldername+"/run.info").c_str());
	infofile << std::setw(20) << "--nx"  << "\t" << D.nx << '\n';
	infofile << std::setw(20) << "--ny"  << "\t" << D.ny << '\n';
	infofile << std::setw(20) << "--start_step" << "\t" << DB["simulation"]["startStep"].get<int>() << '\n';
	infofile << std::setw(20) << "--nt"     << "\t" << DB["simulation"]["nt"].get<int>() << '\n';
	infofile << std::setw(20) << "--nsave"  << "\t" << DB["simulation"]["nsave"].get<int>() << '\n';
	infofile << std::setw(20) << "--dt"     << "\t" << DB["simulation"]["dt"].get<real>() << '\n';
	infofile << std::setw(20) << "--vortlim"<< "\t" << 15 << '\n';
	infofile << std::setw(20) << "--folder" << "\t" << foldername << '\n';
	infofile << std::setw(20) << "--nu"     << "\t" << DB["flow"]["nu"].get<real>() << '\n';
	infofile << std::setw(20) << "--flow_file" << "\t" << DB["inputs"]["flowFile"].get<string>() << '\n';
	infofile << std::setw(20) << "--sim_file" << "\t"  << DB["inputs"]["simulationFile"].get<string>() << '\n';
	infofile << std::setw(20) << "--dom_file" << "\t"  << DB["inputs"]["domainFile"].get<string>() << '\n';
	infofile << std::setw(20) << "--body_file" << "\t" << DB["inputs"]["bodyFile"].get<string>() << '\n';
	infofile.close();
}

void writeGrid(std::string &folderName, domain &D)
{
	std::stringstream out;
	out << folderName << "/grid";
	std::ofstream f(out.str().c_str());

	f << D.nx << std::endl;
	for(int i=0; i<D.nx+1; i++)
		f << D.x[i] << '\n';
	f << '\n';
	
	f << D.ny << '\n';
	for(int j=0; j<D.ny+1; j++)
		f << D.y[j] << '\n';

	f.close();
}

template <>
void writeData<vecH>(std::string &folderName, int n, vecH &q, vecH &lambda, domain &D)//, bodies &B)
{
	std::string path;
	std::stringstream out;

	out << folderName << '/' << std::setfill('0') << std::setw(7) << n;
	path = out.str();

	//createDirectory(path.c_str(), S_IRWXO);
	mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	out.str("");
	out << path << "/q";
	std::ofstream file(out.str().c_str());
	file << q.size() << std::endl;
	for(int i=0; i<q.size(); i++)
		file << q[i] << '\n';
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
	file.open(out.str().c_str());
	file << lambda.size() << '\n';
	for(int i=0; i<lambda.size(); i++)
		file << lambda[i] << '\n';
	file.close();
	
	std::cout << "Data saved to folder " << path << std::endl;
}

template <>
void writeData<vecD>(std::string &folderName, int n, vecD &q, vecD &lambda, domain &D)//, bodies &B)
{
	vecH qH = q,
	     lambdaH = lambda;
	     
	writeData(folderName, n, qH, lambdaH, D);
}

void printDeviceMemoryUsage(char *label)
{
	size_t _free, _total;
	cudaMemGetInfo(&_free, &_total);
	std::cout << '\n' << label << ": Memory Usage " << std::setprecision(3) << (_total-_free)/(1024.0*1024*1024) \
	          << " / " << std::setprecision(3) << _total/(1024.0*1024*1024) << " GB" << std::setprecision(6) << '\n' << std::endl;
}

} // end namespace io
