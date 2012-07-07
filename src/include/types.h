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

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>

#include <cusp/csr_matrix.h>
#include <cusp/array1d.h>
#include <cusp/print.h>

///
#include <cusp/transpose.h>
#include <cusp/blas.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/gmres.h>

#include <cusp/wrapped/add.h>
#include <cusp/wrapped/subtract.h>
#include <cusp/wrapped/multiply.h>

#include <thrust/device_ptr.h>

#include <logger.h>

enum bcType {DIRICHLET, NEUMANN, CONVECTIVE, PERIODIC};
enum boundary {XMINUS, XPLUS, YMINUS, YPLUS};
enum timeScheme {EULER_EXPLICIT, EULER_IMPLICIT, ADAMS_BASHFORTH_2, RUNGE_KUTTA_3, CRANK_NICOLSON};
enum ibmScheme {NAVIER_STOKES, SAIKI_BIRINGEN, FADLUN_ET_AL, TAIRA_COLONIUS};
enum preconditionerType {NONE, DIAGONAL, SMOOTHED_AGGREGATION};

typedef double real;

using cusp::device_memory;
using cusp::host_memory;
using cusp::array1d;
using cusp::coo_matrix;
using cusp::csr_matrix;

typedef coo_matrix<int, real, host_memory> cooH;
typedef coo_matrix<int, real, device_memory> cooD;
typedef csr_matrix<int, real, host_memory> csrH;
typedef csr_matrix<int, real, device_memory> csrD;

typedef array1d<real, host_memory> vecH;
typedef array1d<real, device_memory> vecD;
