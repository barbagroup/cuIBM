#pragma once
/**
* Defines the domain and grid on which the problem is solved.
*/

class domain
{
public:
	int    nx, ny;
	vecH  x, y;
	vecH  dx, dy;
	vecD  xD, yD;
	vecD  dxD, dyD;
	domain(real X0=-5.0, real X1=5.0, real Y0=-5.0, real Y1=5.0, int NX=100, int NY=100)
	{
		nx = NX;
		ny = NY;
		x.resize(nx+1);
		y.resize(ny+1);
		dx.resize(nx);
		dy.resize(ny);
		real DX = (X1-X0)/NX;
		real DY = (Y1-Y0)/NY;
		x[0] = X0;
		for(int i=0; i<nx; i++)
		{
			dx[i]  = DX;
			x[i+1] = x[i] + dx[i];
		}
		y[0] = Y0;
		for(int j=0; j<ny; j++)
		{
			dy[j]  = DY;
			y[j+1] = y[j] + dy[j];
		}
		xD = x;
		yD = y;
		dxD = dx;
		dyD = dy;
	}
};