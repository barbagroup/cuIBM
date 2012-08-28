/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

#include <bodies.h>

template <typename memoryType>
void bodies<memoryType>::initialise(parameterDB &db, domain &D)
{
	std::cout << "Initialising bodies... ";
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();

	numBodies = B->size();

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
	ones.resize(totalPoints);
	cusp::blas::fill(ones, 1.0);
	for(int k=0; k<numBodies; k++)
	{
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

	bodiesMove = false;
	for(int k=0; k<numBodies; k++)
	{
		for(int i=offsets[k], j = offsets[k]+numPoints[k]-1; i<offsets[k]+numPoints[k];)
		{
			// calculate the boundary segments lengths
			ds[i] = sqrt( (X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j]) );
			// set the initial x-coordinates of the boundary points
			x[i] = (*B)[k].Xc[0] + ((*B)[k].X[i]-(*B)[k].X0[0])*cos((*B)[k].Theta0) - ((*B)[k].Y[i] - (*B)[k].X0[1])*sin((*B)[k].Theta0);
			// set the initial y-coordinates of the boundary points
			y[i] = (*B)[k].Xc[1] + ((*B)[k].X[i]-(*B)[k].X0[0])*sin((*B)[k].Theta0) + ((*B)[k].Y[i] - (*B)[k].X0[1])*cos((*B)[k].Theta0);
			uB[i] = 0.0;
			vB[i] = 0.0;
			j = i++;
		}
		// if the body is moving, set bodiesMove to true
		bodiesMove = bodiesMove || (*B)[k].moving[0] || (*B)[k].moving[1];
	}
	if(numBodies)
	{
		calculateCellIndices(D);
		calculateBoundingBoxes(D);
	}
	std::cout << "DONE!" << std::endl;
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
		// if the next boundary point is to the left of the current boundary point
		if(x[k] < x[k-1])
		{
			while(D.x[i] > x[k])
				i--;
		}
		// if the next boundary point is to the right of the current boundary point
		else
		{
			while(D.x[i+1] < x[k])
				i++;
		}
		// if the next boundary point is below the current boundary point
		if(y[k] < y[k-1])
		{
			while(D.y[j] > y[k])
				j--;
		}
		// if the next boundary point is above the current boundary point
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
void bodies<memoryType>::update(parameterDB &db, domain &D, real Time)
{
	
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();
	
	// Update the location and velocity of the body
	(*B)[0].update(Time);

	// Update postitions	
	// x-coordinates
	cusp::blas::axpbypcz( ones, X, ones, x, (*B)[0].Xc[0],  cos((*B)[0].Theta), -(*B)[0].X0[0]*cos((*B)[0].Theta) );
	cusp::blas::axpbypcz( x,    Y, ones, x,           1.0, -sin((*B)[0].Theta),  (*B)[0].X0[1]*sin((*B)[0].Theta) );
	/// y-coordinates
	cusp::blas::axpbypcz( ones, X, ones, y, (*B)[0].Xc[1],  sin((*B)[0].Theta), -(*B)[0].X0[0]*sin((*B)[0].Theta) );
	cusp::blas::axpbypcz( y,    Y, ones, y,           1.0,  cos((*B)[0].Theta), -(*B)[0].X0[1]*cos((*B)[0].Theta) );
	
	// Update velocities
	// x-velocities
	cusp::blas::axpbypcz(ones, y, ones, uB, (*B)[0].vel[0], -(*B)[0].angVel,  (*B)[0].angVel*(*B)[0].Xc[1]);
	// y-velocities
	cusp::blas::axpbypcz(ones, x, ones, vB, (*B)[0].vel[1],  (*B)[0].angVel, -(*B)[0].angVel*(*B)[0].Xc[0]);
	
	calculateCellIndices(D);
}

template class bodies<host_memory>;
template class bodies<device_memory>;
