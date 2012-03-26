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
	simulationParameters(real DT=0.01, int NT=1000, int NSAVE=1000, bool RESTART=false, int START_STEP=0, ibmScheme IBM_SCH=TAIRA_COLONIUS)
	{
		dt = DT;
		nt = NT;
		nsave = NSAVE;
		restart = RESTART;
		startStep = START_STEP;
		ibmSch = IBM_SCH;
	}
};