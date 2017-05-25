/**
 * \file TairaColoniusSolver.cu
 * \brief Implementation of the methods of the class \c TairaColoniusSolver.
 */


#include <sys/stat.h>

#include "TairaColoniusSolver.h"


/**
 * \brief Initializes the solvers.
 *
 * Initializes bodies, arrays and matrices of the intermediate flux solver
 * and the Poisson solver.
 *
 */
template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
      ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1),
	    numP = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	
	NSWithBody<memoryType>::initialiseBodies();
	int totalPoints  = NSWithBody<memoryType>::B.totalPoints; 
	
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP+2*totalPoints);
	if(totalPoints > 0)
	{
		// for each Lagragian point, 24 velocity nodes around it
		// are influenced by the delta function at the point
		// A 4x3 grid of u-nodes, and a 3x4 grid of v-nodes
		// which gives 12*2 velocity nodes influence by a boundary point in 2D
		E.resize(2*totalPoints, numUV, 24*totalPoints);
	}
	NavierStokesSolver<memoryType>::assembleMatrices();
}


/**
 * \brief Calculates and writes forces acting on each immersed body at current time.
 */
template <typename memoryType>
void TairaColoniusSolver<memoryType>::writeData()
{	
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	NSWithBody<memoryType>::writeCommon();

	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	real dt = db["simulation"]["dt"].get<real>();
	int numBodies = NSWithBody<memoryType>::B.numBodies;

	// calculate forces using the T&C method
	calculateForce();
	
	// print to file
	NSWithBody<memoryType>::forceFile << NavierStokesSolver<memoryType>::timeStep*dt << '\t';
	for(int l=0; l<numBodies; l++)
	{
		NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::B.forceX[l] << '\t' << NSWithBody<memoryType>::B.forceY[l] << '\t';
	}
	NSWithBody<memoryType>::forceFile << std::endl;
	
	// print forces calculated using the CV approach
	//NSWithBody<memoryType>::calculateForce();
	//NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << std::endl;
	
	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}


/**
 * \brief Updates the location of the bodies and re-generates appropriate matrices.
 *
 * Updates only done if the bodies are moving.
 * It re-generates the interpolation matrix, therefore the Poisson matrix too.
 *
 */
template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateSolverState()
{
	if (NSWithBody<memoryType>::B.bodiesMove)
	{	
		NSWithBody<memoryType>::updateBodies();
		updateQT();
		NavierStokesSolver<memoryType>::generateC();
		
		NavierStokesSolver<memoryType>::logger.startTimer("preconditioner2");
		NavierStokesSolver<memoryType>::PC2->update(NavierStokesSolver<memoryType>::C);
		NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner2");
	}
}


/**
 * \brief Constructor. Copies the database and information about the computational grid.
 */
template <typename memoryType>
TairaColoniusSolver<memoryType>::TairaColoniusSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}


// include inline files located in "./TairaColonius/"
#include "TairaColonius/generateQT.inl"
#include "TairaColonius/generateBC2.inl"
#include "TairaColonius/calculateForce.inl"


// specialization of the class TairaColoniusSolver
template class TairaColoniusSolver<host_memory>;
template class TairaColoniusSolver<device_memory>;
