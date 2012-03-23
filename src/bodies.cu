#include <bodies.h>

void bodies::initialise(flowDescription &F, domain &D)
{
	std::cout << "Entered B.initialise" << std::endl;
	numBodies = F.numBodies;
	
	X0_x = new real[numBodies];
	X0_y = new real[numBodies];
	Theta0 = new real[numBodies];
	Xc_x = new real[numBodies];
	Xc_y = new real[numBodies];
	numPoints = new int[numBodies];
	offsets   = new int[numBodies];
	
	totalPoints = 0;
	for(int k=0; k<numBodies; k++)
	{
		offsets[k] = totalPoints;
		numPoints[k] = F.B[k].numPoints;
		totalPoints += numPoints[k];
	}
	X  = new real[totalPoints];
	Y  = new real[totalPoints];
	ds = new real[totalPoints];
	for(int k=0; k<numBodies; k++)
	{
		X0_x[k] = F.B[k].X0_x;
		X0_y[k] = F.B[k].X0_y;
		Theta0[k] = 0.0;
		for(int i=0; i<numPoints[k]; i++)
		{
			X[i+offsets[k]] = F.B[k].X[i];
			Y[i+offsets[k]] = F.B[k].Y[i];
		}
	}
	x  = new real[totalPoints];
	y  = new real[totalPoints];
	uB = new real[totalPoints];
	vB = new real[totalPoints];
	I  = new int[totalPoints];
	J  = new int[totalPoints];
	
	for(int k=0; k<numBodies; k++)
	{
		Xc_x[k] = X0_x[k]; ////
		Xc_y[k] = X0_y[k]; ////
		for(int i=offsets[k], j = offsets[k]+numPoints[k]-1; i<offsets[k]+numPoints[k];)
		{
			ds[i] = sqrt( (X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j]) );
			x[i] = X[i];
			y[i] = Y[i];
			uB[i] = 0.0;
			vB[i] = 0.0;
			j = i++;
		}
	}
	bodiesMove = false;
	calculateCellIndices(D);
}

void bodies::calculateCellIndices(domain &D)
{
	std::cout << "Entered calculateCellIndices" << std::endl;
	int	i=0, j=0;

	/// find the cell for the zeroth point
	while(D.x[i+1] < x[0])
		i++;
	while(D.y[j+1] < y[0])
		j++;	
	I[0] = i;
	J[0] = j;

	for(int k=1; k<totalPoints; k++)
	{
		/// if the next boundary point is to the left of the current boundary point
		if(x[k] < x[k-1])
		{
			while(D.x[i] > x[k])
				i--;
		}
		/// if the next boundary point is to the right of the current boundary point
		else
		{
			while(D.x[i+1] < x[k])
				i++;
		}
		/// if the next boundary point is below the current boundary point
		if(y[k] < y[k-1])
		{
			while(D.y[j] > y[k])
				j--;
		}
		/// if the next boundary point is above the current boundary point
		else
		{
			while(D.y[j+1] < y[k])
				j++;
		}
		I[k] = i;
		J[k] = j;
	}
}

void bodies::update()
{
}