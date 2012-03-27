#pragma once
/**
* Defines the domain and grid on which the problem is solved.
*/

class domain
{
public:
	int   nx, ny;
	vecH  x, y,
	      dx, dy;
	vecD  xD, yD,
	      dxD, dyD;
};