/***************************************************************************//**
 * \file types.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of custom types required by the code.
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
#include <cusp/blas/blas.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/gmres.h>
#include <cusp/multiply.h>
#include <cusp/elementwise.h>

#include <thrust/device_ptr.h>

#include "logger.h"


/**
 * \enum  bcType
 * \brief Specifies the type of boundary condition.
 */
enum bcType
{
	DIRICHLET,  ///< Dirichlet boundary condition
	NEUMANN,    ///< Neumann boundary condition
	CONVECTIVE, ///< convective boundary condition
	PERIODIC,   ///< periodic boundary condition
	SPECIAL
};

/**
 * \enum  boundary
 * \brief Specifies the location of the boundary.
 */
enum boundary
{
	XMINUS, ///< left boundary
	XPLUS,  ///< right boundary
	YMINUS, ///< bottom boundary
	YPLUS   ///< top boundary
};

/**
 * \enum  timeScheme
 * \brief Specifies the numerical scheme used for time-integration.
 */
enum timeScheme
{
	EULER_EXPLICIT,    ///< explicit Euler method (first order)
	EULER_IMPLICIT,    ///< implicit Euler method (first order)
	ADAMS_BASHFORTH_2, ///< second-order Adams-Bashforth scheme
	RUNGE_KUTTA_3,     ///< third-order low storage Runge-Kutta method
	CRANK_NICOLSON     ///< Crank-Nicolson scheme (second order)
};

/**
 * \enum  ibmScheme
 * \brief Specifies the immersed boundary method used to solve the flow.
 */
enum ibmScheme
{
	NAVIER_STOKES,  ///< no immersed bodies - Perot (1993)
	SAIKI_BIRINGEN, ///< Saiki & Biringen (1996)
	DIRECT_FORCING, ///< fully discrete direct forcing method
	FADLUN_ET_AL,   ///< Fadlun et al (2000)
	TAIRA_COLONIUS, ///< Taira & Colonius (2007)
	DIFFUSION,      ///< Diffusion
	DF_MODIFIED,    ///< Direct Forcing modified
	FEA_MODIFIED,   ///< Fadlun et al. modified
	DF_IMPROVED     ///< Direct Forcing Improved
};

/**
 * \enum interpolationType
 * \brief Specifies the type of interpolation.
 */
enum interpolationType
{
	CONSTANT,
	LINEAR,
	QUADRATIC
};

/**
 * \enum  preconditionerType
 * \brief Specifies the type of preconditioner.
 */
enum preconditionerType
{
	NONE,                 ///< no preconditioner
	DIAGONAL,             ///< diagonal preconditioner
	SMOOTHED_AGGREGATION, ///< smoothed aggregation preconditioner
	AINV                  ///< approximate inverse preconditioner
};

/**
 * \typedef real
 * \brief Is a \c float or a \c double depending on the machine precision.
 */
typedef double real;

using cusp::device_memory;
using cusp::host_memory;

/**
 * \typedef cooH
 * \brief COO matrix stored on the host.
 */
typedef cusp::coo_matrix<int, real, host_memory>   cooH;

/**
 * \typedef cooD
 * \brief COO matrix stored on the device.
 */
typedef cusp::coo_matrix<int, real, device_memory> cooD;

/**
 * \typedef csrH
 * \brief CSR matrix stored on the host.
 */
typedef cusp::csr_matrix<int, real, host_memory>   csrH;

/**
 * \typedef csrD
 * \brief CSR matrix stored on the device.
 */
typedef cusp::csr_matrix<int, real, device_memory> csrD;

/**
 * \typedef vecH
 * \brief Cusp 1D array stored in the host.
 */
typedef cusp::array1d<real, host_memory>   vecH;

/**
 * \typedef vecD
 * \brief Cusp 1D array stored in the device.
 */
typedef cusp::array1d<real, device_memory> vecD;
