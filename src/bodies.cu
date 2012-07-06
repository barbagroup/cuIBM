#include <bodies.h>

template <typename memoryType>
void bodies<memoryType>::initialise(parameterDB &db, domain &D)
{
	std::cout << "Entered B.initialise" << std::endl;
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();

	numBodies = B->size();

	X0_x.resize(numBodies);
	X0_y.resize(numBodies);
	Xc_x.resize(numBodies);
	Xc_y.resize(numBodies);
	Theta0.resize(numBodies);
	numPoints.resize(numBodies);
	offsets.resize(numBodies);
	
	startI.resize(numBodies);
	startJ.resize(numBodies);
	numCellsX.resize(numBodies);
	numCellsY.resize(numBodies);
	
	xmin.resize(numBodies);
	xmax.resize(numBodies);
	ymin.resize(numBodies);
	ymax.resize(numBodies);

	totalPoints = 0;
	for(int k=0; k<numBodies; k++)
	{
		offsets[k] = totalPoints;
		numPoints[k] = (*B)[k].numPoints;
		totalPoints += numPoints[k];
	}
	
	X.resize(totalPoints);
	Y.resize(totalPoints);
	ds.resize(totalPoints);
	for(int k=0; k<numBodies; k++)
	{
		X0_x[k] = (*B)[k].X0[0];
		X0_y[k] = (*B)[k].X0[1];
		Theta0[k] = (*B)[k].Theta0;
		for(int i=0; i<numPoints[k]; i++)
		{
			X[i+offsets[k]] = (*B)[k].X[i];
			Y[i+offsets[k]] = (*B)[k].Y[i];
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
	{
		calculateCellIndices(D);
		calculateBoundingBoxes(D);
	}
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
/**
* \brief This function determines the bounding box around each body.
*/
template <typename memoryType>
void bodies<memoryType>::calculateBoundingBoxes(domain &D)
{
	real scale=2.0, dx, dy;
	int  i, j;
	for(int k=0; k<numBodies; k++)
	{
		xmin[k] = x[offsets[k]];
		xmax[k] = xmin[k];
		ymin[k] = y[offsets[k]];
		ymax[k] = ymin[k];
		for(int l=offsets[k]+1; l<offsets[k]+numPoints[k]; l++)
		{
			if(x[l] < xmin[k]) xmin[k] = x[l];
			if(x[l] > xmax[k]) xmax[k] = x[l];
			if(y[l] < ymin[k]) ymin[k] = y[l];
			if(y[l] > ymax[k]) ymax[k] = y[l];
		}
		dx = xmax[k]-xmin[k];
		dy = ymax[k]-ymin[k];
		xmax[k] += 0.5*dx*(scale-1.0);
		xmin[k] -= 0.5*dx*(scale-1.0);
		ymax[k] += 0.5*dy*(scale-1.0);
		ymin[k] -= 0.5*dy*(scale-1.0);
		
		i=0; j=0;
		while(D.x[i+1] < xmin[k])
			i++;
		while(D.y[j+1] < ymin[k])
			j++;
		startI[k] = i;
		startJ[k] = j;
		
		while(D.x[i] < xmax[k])
			i++;
		while(D.y[j] < ymax[k])
			j++;
		numCellsX[k] = i - startI[k];
		numCellsY[k] = j - startJ[k];
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
