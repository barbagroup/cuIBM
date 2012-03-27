#include <bodies.h>

template <typename memoryType>
void bodies<memoryType>::initialise(flowDescription &F, domain &D)
{
	std::cout << "Entered B.initialise" << std::endl;
	numBodies = F.numBodies;
	
	X0_x.resize(numBodies);
	X0_y.resize(numBodies);
	Xc_x.resize(numBodies);
	Xc_y.resize(numBodies);
	Theta0.resize(numBodies);
	numPoints.resize(numBodies);
	offsets.resize(numBodies);
	
	totalPoints = 0;
	for(int k=0; k<numBodies; k++)
	{
		offsets[k] = totalPoints;
		numPoints[k] = F.B[k].numPoints;
		totalPoints += numPoints[k];
	}
	X.resize(totalPoints);
	Y.resize(totalPoints);
	ds.resize(totalPoints);
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
	x.resize(totalPoints);
	y.resize(totalPoints);
	uB.resize(totalPoints);
	vB.resize(totalPoints);
	I.resize(totalPoints);
	J.resize(totalPoints);
	
	for(int k=0; k<numBodies; k++)
	{
		Xc_x[k] = X0_x[k];
		Xc_y[k] = X0_y[k];
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
	if(numBodies)
		calculateCellIndices(D);
	std::cout << "Completed B.initialise" << std::endl;
}

/*
template <typename memoryType>
void bodies<memoryType>::initialise(flowDescription &F, domain &D)
{
}
*/

template <typename memoryType>
void bodies<memoryType>::calculateCellIndices(domain &D)
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
/*
template <>
void bodies<device_memory>::calculateCellIndices(domain &D)
{
}*/

template <typename memoryType>
void bodies<memoryType>::update()
{
}

template class bodies<host_memory>;
template class bodies<device_memory>;