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
	bool  restart;
	int   startStep;
	integrationScheme intSch;
	ibmScheme ibmSch;
	simulationParameters(real DT=0.02, int NT=100, int NSAVE=100, bool RESTART=false, int START_STEP=0)
	{
		dt = DT;
		nt = NT;
		nsave = NSAVE;
		restart = RESTART;
		startStep = START_STEP;
	}
};