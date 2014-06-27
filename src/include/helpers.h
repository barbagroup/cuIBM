/***************************************************************************//**
* \file helpers.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of discrete Delta function
*/

#pragma once

#include <types.h>

/********************//**
\brief Discrete Delta function from Roma et al. (1999)
*/
real dhRoma(real x, real h);

/********************//**
\brief Discrete Delta function in 2D
*/
real delta(real x, real y, real h);
