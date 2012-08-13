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

#include <types.h>
#include <domain.h>
#include <io/io.h>
#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>

int main(int argc, char **argv)
{
	std::cout << std::endl;
	
	domain dom_info;
	parameterDB paramDB;

	cudaSetDevice(0);

	io::readInputs(argc, argv, paramDB, dom_info);
	io::printSimulationInfo(paramDB, dom_info);

	/// choose the appropriate flow solver
	NavierStokesSolver<device_memory> *solver = NavierStokesSolver<device_memory>::createSolver(paramDB, dom_info);
	solver->initialise();
	
	io::writeInfoFile(paramDB, dom_info);
	
	while (!solver->finished())
	{
		solver->stepTime();
		solver->writeData();
	}
	solver->shutDown();
}
