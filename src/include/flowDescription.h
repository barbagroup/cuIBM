#pragma once

#include <types.h>
#include <body.h>
#include <utility>
/**
* Contains the description of the flow problem.
*/
class flowDescription
{
public:
	real nu;
	real initialU,
	     initialV;
	
	int numBodies;
	std::vector<body> B;
	std::pair <bcType, real> bcInfo[4];
};