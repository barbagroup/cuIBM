/***************************************************************************//**
 * \file helpers.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the discrete delta function.
 */


#pragma once

#include "types.h"


// discrete delta function used by Roma et al. (1999)
real dhRoma(real x, real h);

// two-dimensional discrete delta function
real delta(real x, real y, real h);
