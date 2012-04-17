#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <sys/stat.h>
/*
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateA()
{
}
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateA()
{
}*/
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateBodies()
{
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	initialiseBodies();
	std::cout << "Initialised bodies!" << std::endl;
	
	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	
	NavierStokesSolver<memoryType>::timeStep = db["simulation"]["startStep"].get<int>();
	
	// create directory
	std::string folderName = db["inputs"]["folderName"].get<std::string>();
	mkdir(folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	// write grid data to file
	io::writeGrid(folderName, *NavierStokesSolver<memoryType>::domInfo);

	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	
	tags.resize(numUV);
	tagsD.resize(numUV);
	coeffs.resize(numUV);
	coeffsD.resize(numUV);
	tagPoints();
	
	NavierStokesSolver<memoryType>::assembleMatrices();
	std::cout << "Assembled matrices!" << std::endl;
}
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateSolverState()
{
	NavierStokesSolver<memoryType>::updateBoundaryConditions();
	if (B.bodiesMove)
	{
		updateBodies();
		//updateA();
	}
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateRN()
{
	NavierStokesSolver<memoryType>::generateRNFull(1.0, 0.0, 0.0);
	updateRN();
	
	/**
	* Added for testing. REMOVE SOON!!!
	*/
	//cusp::blas::fill(NavierStokesSolver<memoryType>::rn, 0.0);
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateBC1()
{
	NavierStokesSolver<memoryType>::generateBC1Full(1.0);
	updateBC1();
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateBC1()
{
}

#include "FadlunEtAl/initialiseBodies.inl"
#include "FadlunEtAl/tagPoints.inl"
#include "FadlunEtAl/generateL.inl"
#include "FadlunEtAl/updateRN.inl"

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;