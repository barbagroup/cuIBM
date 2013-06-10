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

#include <solvers/NavierStokes/SLL0Solver.h>
#include <sys/stat.h>
#include <cusp/io/matrix_market.h>

template <typename memoryType>
void SLL0Solver<memoryType>::assembleRHS1()
{
	NavierStokesSolver<memoryType>::logger.startTimer("assembleRHS1");
	
	cusp::blas::axpby(NavierStokesSolver<memoryType>::rn, NavierStokesSolver<memoryType>::bc1, NavierStokesSolver<memoryType>::rhs1, 1.0, 1.0);
	cusp::multiply(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::lambda, NavierStokesSolver<memoryType>::temp1);
	cusp::blas::axpby(NavierStokesSolver<memoryType>::rhs1, NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rhs1, 1.0, -1.0);
	
	NavierStokesSolver<memoryType>::logger.stopTimer("assembleRHS1");
}

template class SLL0Solver<host_memory>;
template class SLL0Solver<device_memory>;
