#include <io/io.h>
#include <types.h>
#include <sys/stat.h>
#include <boundaryCondition.h>

using std::string;

namespace io
{

void initialiseDefaultDB(parameterDB &DB)
{
  DB["inputs"] = componentParameter();
  DB["flow"] = componentParameter();
  DB["simulation"] = componentParameter();
  DB["velocitySolve"] = componentParameter();
  DB["PoissonSolve"] = componentParameter();

  // default input files
  string inputs = "inputs";
  DB[inputs]["flowFile"].set<string>("flows/open_flow.yaml");
  DB[inputs]["simulationFile"].set<string>("test.yaml");
  DB[inputs]["bodyFile"].set<string>("inputs/cylinder.yaml");
  DB[inputs]["domainFile"].set<string>("domains/open_flow.yaml");
  DB[inputs]["folderName"].set<string>("new");

  // flow parameters
  string flow = "flow";
  DB[flow]["nu"].set<double>(0.01);
  DB[flow]["uInitial"].set<double>(1);
  DB[flow]["vInitial"].set<double>(0);
  DB[flow]["numBodies"].set<int>(0);
  std::vector<body> *bodyVec = new std::vector<body>;
  DB[flow]["bodies"].set<std::vector<body> *>(bodyVec);

  boundaryCondition **bc = new boundaryCondition*[4];
  for (int i=0; i<4; i++) bc[i] = new boundaryCondition[2];
  DB[flow]["boundaryConditions"].set<boundaryCondition **>(bc);

  // simulation parameters
  string sim = "simulation";
  DB[sim]["dt"].set<double>(0.02);
  DB[sim]["nt"].set<int>(100);
  DB[sim]["nsave"].set<int>(100);
  DB[sim]["restart"].set<bool>(false);
  DB[sim]["startStep"].set<bool>(0);
  DB[sim]["convTimeScheme"].set<timeScheme>(EULER_EXPLICIT);
  DB[sim]["diffTimeScheme"].set<timeScheme>(EULER_IMPLICIT);

  // velocity solver
  string solver = "velocitySolve";
  DB[solver]["solver"].set<string>("CG");
  DB[solver]["preconditioner"].set<preconditionerType>(DIAGONAL);
  DB[solver]["tolerance"].set<double>(1e-5);
  DB[solver]["maxIterations"].set<int>(10000);

  // Poisson solver
  solver = "PoissonSolve";
  DB[solver]["solver"].set<string>("CG");
  DB[solver]["preconditioner"].set<preconditionerType>(DIAGONAL);
  DB[solver]["tolerance"].set<double>(1e-5);
  DB[solver]["maxIterations"].set<int>(20000);

}

// first pass -- only get the files to continue
void commandLineParse1(int argc, char **argv, parameterDB &DB)
{
  for (int i=1; i<argc; i++)
  {
    if (strcmp(argv[i],"-flowFile")==0)
    {
      i++;
      DB["inputs"]["flowFile"].set<string>(string(argv[i]));
    }
    else if (strcmp(argv[i],"simulationFile")==0)
    {
      i++;
      DB["inputs"]["simulationFile"].set<string>(string(argv[i]));
    }
    else if (strcmp(argv[i],"bodyFile")==0)
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
  }
}

// overwrite values in the DB from command line
void commandLineParse2(int argc, char **argv, parameterDB &DB)
{
  for (int i=1; i<argc; i++)
  {
    // ignore these -- already parsed in pass 1
    if (strcmp(argv[i],"-flowFile")==0 ||
        strcmp(argv[i],"-simulationFile")==0 ||
        strcmp(argv[i],"-bodyFile")==0 ||
        strcmp(argv[i],"-domainFile")==0 ||
        strcmp(argv[i],"-folderName")==0) continue;
  }
}

void readInputs(int argc, char **argv, parameterDB &DB, domain &D)
{
  std::vector<body> B;
  B.reserve(5);
  // get a default database
  initialiseDefaultDB(DB);
  // first pass of command line arguments
  commandLineParse1(argc,argv,DB);
  // read the simulation file
  string fname = DB["inputs"]["domainFile"].get<string>();
  parseSimulationFile(fname,DB);
  // read the flow file
  fname = DB["inputs"]["flowFile"].get<string>();
  parseFlowFile(fname,DB);
  // read the domain file
  fname = DB["inputs"]["domainFile"].get<string>();
  parseDomainFile(fname,D);
  // read the body file
  fname = DB["inputs"]["bodyFile"].get<string>();
  parseBodiesFile(fname,DB);
  // second pass of command line -- overwrite values in DB
  commandLineParse2(argc,argv,DB);
}

// output
void printSimulationInfo(parameterDB &DB, domain &D)
{
  /*
  std::cout << std::endl;
  std::cout << "Simulation parameters" << std::endl;
  std::cout << "---------------------" << std::endl;
  double dt = DB["simulation"]["dt"].get<double>();
  std::cout << "dt = " << dt         << std::endl;
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
  */
}

void writeGrid(std::string &folderName, domain &D)
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
void writeData<vecH>(std::string &folderName, int n, vecH &q, vecH &lambda, domain &D)//, bodies &B)
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
void writeData<vecD>(std::string &folderName, int n, vecD &q, vecD &lambda, domain &D)//, bodies &B)
{
  vecH qH = q,
       lambdaH = lambda;
  writeData(folderName, n, qH, lambdaH, D);
}

} // end namespace io
