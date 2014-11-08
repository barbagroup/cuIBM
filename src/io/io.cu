#include "io.h"
#include <types.h>
#include <sys/stat.h>
#include <boundaryCondition.h>

using std::string;
using std::ios;

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

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

//##############################################################################
//                                 INPUT
//##############################################################################
	
void readInputs(int argc, char **argv, parameterDB &DB, domain &D)
{
	// get a default database
	initialiseDefaultDB(DB);
	
	// first pass of command line arguments
	commandLineParse1(argc, argv, DB);
	
	// case folder
	string folder = DB["inputs"]["caseFolder"].get<string>();
	
	// read the simulation file
	string fname = folder + "/simParams.yaml";
	parseSimulationFile(fname, DB);
	
	// read the flow file
	fname = folder + "/flow.yaml";
	parseFlowFile(fname, DB);

	// read the domain file
	fname = folder + "/domain.yaml";
	parseDomainFile(fname, D);
	
	// read the body file
	fname = folder + "/bodies.yaml";
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
	DB[inputs]["caseFolder"].set<string>("cases/cylinder/Re40");
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
	DB[sim]["interpolationType"].set<interpolationType>(LINEAR);

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
		if (strcmp(argv[i],"-caseFolder")==0)
		{
			i++;
			DB["inputs"]["caseFolder"].set<string>(string(argv[i]));
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
		// kinematic viscosity
		if ( strcmp(argv[i],"-nu")==0 )
		{
			i++;
			DB["flow"]["nu"].set<real>(toNumber<real>(string(argv[i])));
		}
		//// angle of attack
		//if ( strcmp(argv[i],"-alpha")==0 )
		//{
		//	i++;
		//	DB["flow"]["nu"].set<real>(toNumber<real>(string(argv[i])));
		//}
		// perturbation in the x-velocity
		if ( strcmp(argv[i],"-uPerturb")==0 )
		{
			i++;
			DB["flow"]["uPerturb"].set<real>(toNumber<real>(string(argv[i])));
		}
		// perturbation in the y-velocity
		if ( strcmp(argv[i],"-vPerturb")==0 )
		{
			i++;
			DB["flow"]["vPerturb"].set<real>(toNumber<real>(string(argv[i])));
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
		// size of time increment
		if ( strcmp(argv[i],"-dt")==0 )
		{
			i++;
			DB["simulation"]["dt"].set<real>(toNumber<real>(string(argv[i])));
		}
		// tolerance for the velocity solve
		if ( strcmp(argv[i],"-velocityTol")==0 )
		{
			i++;
			DB["velocitySolve"]["tolerance"].set<real>(toNumber<real>(string(argv[i])));
		}
		// tolerance for the Poisson solve
		if ( strcmp(argv[i],"-poissonTol")==0 )
		{
			i++;
			DB["PoissonSolve"]["tolerance"].set<real>(toNumber<real>(string(argv[i])));
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
			if ( strcmp(argv[i],"DirectForcing")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(DIRECT_FORCING);
			else 
			if ( strcmp(argv[i],"FadlunEtAl")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(FADLUN_ET_AL);
			else 
			if ( strcmp(argv[i],"SLL0")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(SLL0);
			else 
			if ( strcmp(argv[i],"SLL1")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(SLL1);
			else 
			if ( strcmp(argv[i],"SLL2")==0 )
				DB["simulation"]["ibmScheme"].set<ibmScheme>(SLL2);
		}
		// interpolation type for Eulerian direct forcing methods
		if ( strcmp(argv[i],"-interpolationType")==0 )
		{
			i++;
			if ( strcmp(argv[i],"constant")==0 )
				DB["simulation"]["interpolationType"].set<interpolationType>(CONSTANT);
			else
			if ( strcmp(argv[i],"linear")==0 )
				DB["simulation"]["interpolationType"].set<interpolationType>(LINEAR);
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

string stringFromTimeScheme(timeScheme s)
{
	if (s == EULER_EXPLICIT)
		return "Explicit Euler Method";
	else if (s == EULER_IMPLICIT)
		return "Implicit Euler Method";
	else if (s == ADAMS_BASHFORTH_2)
		return "2nd Order Adams-Bashforth";
	else if (s == CRANK_NICOLSON)
		return "Crank-Nicolson";
	else
		return "Unknown";
}

// output
void printSimulationInfo(parameterDB &DB, domain &D)
{
	real dt = DB["simulation"]["dt"].get<real>(),
	     scaleCV = DB["simulation"]["scaleCV"].get<real>();
	int  nt = DB["simulation"]["nt"].get<int>(),
	     nsave = DB["simulation"]["nsave"].get<int>(),
	     startStep = DB["simulation"]["startStep"].get<int>();
	interpolationType interpType = DB["simulation"]["interpolationType"].get<interpolationType>();
	ibmScheme ibmSchm = DB["simulation"]["ibmScheme"].get<ibmScheme>();


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
	std::cout << "Convection time scheme = " << stringFromTimeScheme(DB["simulation"]["convTimeScheme"].get<timeScheme>()) << '\n';
	std::cout << "Diffusion time scheme  = " << stringFromTimeScheme(DB["simulation"]["diffTimeScheme"].get<timeScheme>()) << '\n';
	if(ibmSchm==FADLUN_ET_AL || ibmSchm==DIRECT_FORCING)
	{
		std::cout << "Interpolation type: ";
		switch(interpType)
		{
			case CONSTANT: std::cout << "Constant\n"; break;
			case LINEAR  : std::cout << "Linear\n"; break;
			default : std::cout << "Unknown\n"; break;
		}
	}
	
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
	std::cout << "Output folder = " << DB["inputs"]["caseFolder"].get<string>() << '\n';
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
	logger.printAllTime();
	std::cout << std::endl;
}

void writeInfoFile(parameterDB &DB, domain &D)
{
	std::string   folder = DB["inputs"]["caseFolder"].get<string>();
	std::ofstream infofile((folder+"/run.info").c_str());
	infofile << std::setw(20) << "--nx"  << "\t" << D.nx << '\n';
	infofile << std::setw(20) << "--ny"  << "\t" << D.ny << '\n';
	infofile << std::setw(20) << "--startStep" << "\t" << DB["simulation"]["startStep"].get<int>() << '\n';
	infofile << std::setw(20) << "--nt"     << "\t" << DB["simulation"]["nt"].get<int>() << '\n';
	infofile << std::setw(20) << "--nsave"  << "\t" << DB["simulation"]["nsave"].get<int>() << '\n';
	infofile << std::setw(20) << "--dt"     << "\t" << DB["simulation"]["dt"].get<real>() << '\n';
	infofile << std::setw(20) << "--vortlim"<< "\t" << 15 << '\n';
	infofile << std::setw(20) << "--folder" << "\t" << folder << '\n';
	infofile << std::setw(20) << "--nu"     << "\t" << DB["flow"]["nu"].get<real>() << '\n';
	infofile.close();
}

void writeGrid(std::string &caseFolder, domain &D)
{
	std::stringstream out;
	out << caseFolder << "/grid";
	std::ofstream file(out.str().c_str(), ios::binary);
	file.write((char*)(&D.nx), sizeof(int));
	file.write((char*)(&D.x[0]), (D.nx+1)*sizeof(real));
	file.write((char*)(&D.ny), sizeof(int));
	file.write((char*)(&D.y[0]), (D.ny+1)*sizeof(real));
	file.close();
}

template <>
void writeData<vecH>(std::string &caseFolder, int n, vecH &q, vecH &lambda, domain &D)//, bodies &B)
{
	std::string path;
	std::stringstream out;
	int N;

	out << caseFolder << '/' << std::setfill('0') << std::setw(7) << n;
	path = out.str();

	mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	out.str("");
	out << path << "/q";
	std::ofstream file(out.str().c_str(), ios::binary);
	N = q.size();
	file.write((char*)(&N), sizeof(int));
	file.write((char*)(&q[0]), N*sizeof(real));
	file.close();

	out.str("");
	out << path << "/lambda";
	file.open(out.str().c_str(), ios::binary);
	N = lambda.size();
	file.write((char*)(&N), sizeof(int));
	file.write((char*)(&lambda[0]), N*sizeof(real));
	file.close();
	
	std::cout << "Data saved to folder " << path << std::endl;
}

template <>
void writeData<vecD>(std::string &caseFolder, int n, vecD &q, vecD &lambda, domain &D)//, bodies &B)
{
	vecH qH = q,
	     lambdaH = lambda;
	     
	writeData(caseFolder, n, qH, lambdaH, D);
}

void printDeviceMemoryUsage(char *label)
{
	size_t _free, _total;
	cudaMemGetInfo(&_free, &_total);
	std::cout << '\n' << label << ": Memory Usage " << std::setprecision(3) << (_total-_free)/(1024.0*1024*1024) \
	          << " / " << std::setprecision(3) << _total/(1024.0*1024*1024) << " GB" << std::setprecision(6) << '\n' << std::endl;
}

} // end namespace io
