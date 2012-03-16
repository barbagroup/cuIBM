#pragma once

#include <cusp/csr_matrix.h>
#include <cusp/array1d.h>

enum bcType { DIRICHLET, NEUMANN, CONVECTIVE, PERIODIC };
enum boundary { XMINUS, XPLUS, YMINUS, YPLUS };
enum timeScheme { EULER_EXPLICIT, EULER_IMPLICIT, ADAMS_BASHFORTH_2, RUNGE_KUTTA_3, CRANK_NICOLSON };
enum ibmScheme { NAVIER_STOKES, SAIKI_BIRINGEN, FADLUN_ET_AL, TAIRA_COLONIUS };

using cusp::device_memory;
using cusp::host_memory;
using cusp::array1d;
using cusp::csr_matrix;

typedef float real;
typedef csr_matrix<int, real, host_memory> matrix;
typedef array1d<real, host_memory> vector;
//typedef memory_space cusp::device_memory
