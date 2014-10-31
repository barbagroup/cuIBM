/**
* @file  types.h
* @brief Custom types required by the code.
*/

#pragma once

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>

#include <cusp/csr_matrix.h>
#include <cusp/array1d.h>
#include <cusp/print.h>

#include <cusp/transpose.h>
#include <cusp/blas.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/gmres.h>

#include <cusp/multiply.h>
#include <cusp/elementwise.h>

#include <thrust/device_ptr.h>

#include "logger.h"

/**
* @enum  bcType
* @brief Specifies the type of boundary condition
*/
enum bcType
{
	DIRICHLET,  ///< Dirichlet boundary condition
	NEUMANN,    ///< Neumann boundary condition
	CONVECTIVE, ///< Convective boundary condition
	PERIODIC,   ///< Periodic boundary condition
	SPECIAL
};

/**
* @enum  boundary
* @brief Specifies the boundary of concern
*/
enum boundary
{
	XMINUS, ///< Left boundary
	XPLUS,  ///< Right boundary
	YMINUS, ///< Bottom boundary
	YPLUS   ///< Top boundary
};

/**
* @enum  timeScheme
* @brief Numerical scheme used to discretise the time derivative
*/
enum timeScheme
{
	EULER_EXPLICIT,    ///< Explicit Euler method (first order)
	EULER_IMPLICIT,    ///< Implicit Euler method (first order)
	ADAMS_BASHFORTH_2, ///< Second-order Adams-Bashforth scheme
	RUNGE_KUTTA_3,     ///< Third-order low storage Runge-Kutta method
	CRANK_NICOLSON     ///< Crank-Nicolson scheme (second order)
};

/**
* @enum  ibmScheme
* @brief The immersed boundary method used to solve for the flow.
*/
enum ibmScheme
{
	NAVIER_STOKES,  ///< No immersed bodies. Perot (1993)
	SAIKI_BIRINGEN, ///< Saiki & Biringen (1996)
	DIRECT_FORCING, ///< Fully discrete direct forcing method
	FADLUN_ET_AL,   ///< Fadlun et al (2000)
	TAIRA_COLONIUS, ///< Taira & Colonius (2007)
	SLL0,           ///< SLL0
	SLL1,           ///< SLL1
	SLL2            ///< SLL2
};

enum interpolationType
{
	CONSTANT,
	LINEAR,
	QUADRATIC
};

/**
* @enum  preconditionerType
* @brief Specify the type of preconditioner
*/
enum preconditionerType
{
	NONE,                ///< No preconditioner
	DIAGONAL,            ///< Diagonal preconditioner
	SMOOTHED_AGGREGATION ///< Smoothed aggregation preconditioner
};

// Shorthand for various types

/**
* @typedef real
* @brief Can be \c float or \c double
*/
typedef double real;

using cusp::device_memory;
using cusp::host_memory;

/**
* @typedef cooH
* @brief   COO matrix stored on the host.
*/
typedef cusp::coo_matrix<int, real, host_memory>   cooH;

/**
* @typedef cooD
* @brief   COO matrix stored on the device.
*/
typedef cusp::coo_matrix<int, real, device_memory> cooD;

typedef cusp::csr_matrix<int, real, host_memory>   csrH;
typedef cusp::csr_matrix<int, real, device_memory> csrD;

/**
* @typedef vecH
* @brief   Cusp array of real numbers stored in the host memory.
*/
typedef cusp::array1d<real, host_memory>   vecH;
/**
* @typedef vecD
* @brief   Cusp array of real numbers stored in the device memory.
*/
typedef cusp::array1d<real, device_memory> vecD;
