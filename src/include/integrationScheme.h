/**
* \file
* \brief
*/
#pragma once

#include <types.h>

/**
* Class description
*/
class integrationScheme
{
public:
	int substeps;
	timeScheme conv_t, diff_t;
};