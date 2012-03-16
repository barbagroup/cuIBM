#pragma once

#include <types.h>
#include <integrationScheme.h>
/**
* Class description
*/
class simulationParameters
{
public:
	real  dt;
	int   nt, nsave;
	int   restart;
	int   startStep;
	integrationScheme intSch;
	ibmScheme ibmSch;
};
