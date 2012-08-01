/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
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

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>

/**
* Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class TairaColoniusSolver : public NSWithBody<memoryType>
{
private:
	coo_matrix<int, real, memoryType> E, ET;
	
	array1d<real, memoryType> FxX, FxY, FxU;
	
	void generateQT(int *QTRows, int *QTCols, real *QTVals){}
	void generateQT();
	void updateQT();
	void generateBC2();
	
	void updateSolverState();
	
	void calculateForce1();
	void calculateForceTC();
	void calculateForce();
	void generateE();
	
public:
	void initialise();
	std::string name()
	{
		return "Taira & Colonius";
	}
};
