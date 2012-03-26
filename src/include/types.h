#pragma once

#include <fstream>
#include <cmath>
#include <cstdio>

#include <cusp/csr_matrix.h>
#include <cusp/array1d.h>
#include <cusp/transpose.h>
#include <cusp/print.h>

#include <cusp/blas.h>
#include <cusp/krylov/cg.h>

#include <cusp/wrapped/add.h>
#include <cusp/wrapped/subtract.h>
#include <cusp/wrapped/multiply.h>

#include <thrust/device_ptr.h>

enum bcType { DIRICHLET, NEUMANN, CONVECTIVE, PERIODIC };
enum boundary { XMINUS, XPLUS, YMINUS, YPLUS };
enum timeScheme { EULER_EXPLICIT, EULER_IMPLICIT, ADAMS_BASHFORTH_2, RUNGE_KUTTA_3, CRANK_NICOLSON };
enum ibmScheme { NAVIER_STOKES, SAIKI_BIRINGEN, FADLUN_ET_AL, TAIRA_COLONIUS };

using cusp::device_memory;
using cusp::host_memory;
using cusp::array1d;
using cusp::csr_matrix;
using cusp::coo_matrix;

typedef float real;
//typedef csr_matrix<int, real, host_memory> matrix;
typedef coo_matrix<int, real, host_memory> cooH;
typedef csr_matrix<int, real, host_memory> csrH;
typedef array1d<real, host_memory> vecH;
typedef coo_matrix<int, real, device_memory> cooD;
typedef csr_matrix<int, real, device_memory> csrD;
typedef array1d<real, device_memory> vecD;
//typedef memory_space cusp::device_memory
typedef cooD matrix;
typedef vecD vector;
typedef array1d<int, host_memory> intVecH;
typedef array1d<int, device_memory> intVecD;
