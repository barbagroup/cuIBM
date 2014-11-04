/***************************************************************************//**
 * \file helpers.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the discrete Delta function.
 */


#pragma once

#include "types.h"


/**
 * \brief Discrete Delta function from Roma et al. (1999).
 */
real dhRoma(real x, real h);

/**
 * \brief Two-dimensional discrete Delta function.
 */
real delta(real x, real y, real h);
