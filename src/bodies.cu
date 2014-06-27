/***************************************************************************//**
* \file bodies.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the class \c bodies
*/

#include <bodies.h>

/***************************************************************************//**
* \brief Initializes the arrays in the class with information from \c body instances.
*
* Information regarding the coordinates of the body points and motion of the bodies
* is stored on the host as an array of instances of the class \c body.
* This function transfers that information to arrays on the device,
* where they are stored as a structur of arrays.
* This makes computation more efficient.
*
* \param db database that contains all the simulation parameters
* \param D Information about the computational grid
*/
template <typename memoryType>
void bodies<memoryType>::initialise(parameterDB &db, domain &D)
{
	std::cout << "Initialising bodies... ";
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();

	// number of bodies in the flow
	numBodies = B->size();

	// set the sizes of all the arrays
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
	
	forceX.resize(numBodies);
	forceY.resize(numBodies);

	// calculate offsets, number of points in each body and the total number of points
	totalPoints = 0;
	for(int k=0; k<numBodies; k++)
	{
		offsets[k] = totalPoints;
		numPoints[k] = (*B)[k].numPoints;
		totalPoints += numPoints[k];
	}
	
	// fill up coordinates of body points
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
		// ASSUMING CLOSED LOOPS
		for(int i=offsets[k], j = offsets[k]+numPoints[k]-1; i<offsets[k]+numPoints[k];)
		{
			// calculate the lengths of the boundary segments
			ds[i] = sqrt( (X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j]) );
			
			j = i++;
		}
		// if the body is moving, set bodiesMove to true
		bodiesMove = bodiesMove || (*B)[k].moving[0] || (*B)[k].moving[1];
	}
	// set initial position of the body
	update(db, D, 0.0);
		
	if(numBodies)
	{
		calculateCellIndices(D);
		calculateBoundingBoxes(db, D);
	}
	std::cout << "DONE!" << std::endl;
}

/***************************************************************************//**
* \brief Calculates the indices of the cells in which the boundary points are present.
*
* It calculates the index of the x-coordinate and the index of the y-coordinate
* of the bottom-left node of the containing cell.
* This information is useful when transferring data between the boundary points
* and the computational grid.
*
* \param D information of the computational grid
*/
template <typename memoryType>
void bodies<memoryType>::calculateCellIndices(domain &D)
{
	int	i=0, j=0;

	// find the cell for the zeroth point
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

/***************************************************************************//**
* \brief Calculates the bounding boxes for each body
* \param D information about the computational grid
*/
template <typename memoryType>
void bodies<memoryType>::calculateBoundingBoxes(parameterDB &db, domain &D)
{
	real scale = db["simulation"]["scaleCV"].get<real>(),
	     dx, dy;
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

/***************************************************************************//**
* \brief Updates the locations of the body points.
*
* This is done using the formulae:
* \f$x_{i,m} = X^c_m + (X_{i,m} - X^0_m) \cos\theta - (Y_{i,m} - Y^0_m) \sin\theta\f$
* and
* \f$y_{i,m} = Y^c_m + (X_{i,m} - X^0_m) \sin\theta + (Y_{i,m} - Y^0_m) \cos\theta\f$
*
* \param db Database containing all the simulation parameters
* \param D Information related to the computational grid
* \param Time the time
*/
template <typename memoryType>
void bodies<memoryType>::update(parameterDB &db, domain &D, real Time)
{
	typedef typename cusp::array1d<real, memoryType> Array;
	typedef typename Array::iterator                 Iterator;
	typedef cusp::array1d_view<Iterator>             View;

	// Views of the vectors that store the coordinates and velocities of all the body points
	View    XView, YView, xView, yView, onesView, uBView, vBView;
	
	// body data
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();
	
	for(int l=0; l<numBodies; l++)
	{
		// Update the location and velocity of the body
		(*B)[l].update(Time);
		
		// create the views for the current body
		if(l < numBodies-1)
		{
			XView    = View(X.begin()+offsets[l], X.begin()+offsets[l+1]);
			YView    = View(Y.begin()+offsets[l], Y.begin()+offsets[l+1]);
			xView    = View(x.begin()+offsets[l], x.begin()+offsets[l+1]);
			yView    = View(y.begin()+offsets[l], y.begin()+offsets[l+1]);
			onesView = View(ones.begin()+offsets[l], ones.begin()+offsets[l+1]);
			uBView   = View(uB.begin()+offsets[l], uB.begin()+offsets[l+1]);
			vBView   = View(vB.begin()+offsets[l], vB.begin()+offsets[l+1]);
		}
		else
		{
			XView    = View(X.begin()+offsets[l], X.end());
			YView    = View(Y.begin()+offsets[l], Y.end());
			xView    = View(x.begin()+offsets[l], x.end());
			yView    = View(y.begin()+offsets[l], y.end());
			onesView = View(ones.begin()+offsets[l], ones.end());
			uBView   = View(uB.begin()+offsets[l], uB.end());
			vBView   = View(vB.begin()+offsets[l], vB.end());
		}
		
		// Update postitions	
		// x-coordinates
		cusp::blas::axpbypcz( onesView, XView, onesView, xView, (*B)[l].Xc[0],  cos((*B)[l].Theta), -(*B)[l].X0[0]*cos((*B)[l].Theta) );
		cusp::blas::axpbypcz( xView,    YView, onesView, xView,           1.0, -sin((*B)[l].Theta),  (*B)[l].X0[1]*sin((*B)[l].Theta) );
		/// y-coordinates
		cusp::blas::axpbypcz( onesView, XView, onesView, yView, (*B)[l].Xc[1],  sin((*B)[l].Theta), -(*B)[l].X0[0]*sin((*B)[l].Theta) );
		cusp::blas::axpbypcz( yView,    YView, onesView, yView,           1.0,  cos((*B)[l].Theta), -(*B)[l].X0[1]*cos((*B)[l].Theta) );
	
		// Update velocities
		// x-velocities
		cusp::blas::axpbypcz(onesView, yView, onesView, uBView, (*B)[l].vel[0], -(*B)[l].angVel,  (*B)[l].angVel*(*B)[l].Xc[1]);
		// y-velocities
		cusp::blas::axpbypcz(onesView, xView, onesView, vBView, (*B)[l].vel[1],  (*B)[l].angVel, -(*B)[l].angVel*(*B)[l].Xc[0]);
	}
	
	if(numBodies)
		calculateCellIndices(D);
}

/***************************************************************************//**
*/
template <>
void bodies<host_memory>::writeToFile(std::string &caseFolder, int timeStep)
{
	 real *bx = thrust::raw_pointer_cast(&(x[0])),
	      *by = thrust::raw_pointer_cast(&(y[0]));
	 writeToFile(bx, by, caseFolder, timeStep);
}

/***************************************************************************//**
*/
template <>
void bodies<device_memory>::writeToFile(std::string &caseFolder, int timeStep)
{
	vecH xHost = x,
	     yHost = y;
	real *bx = thrust::raw_pointer_cast(&(xHost[0])),
	     *by = thrust::raw_pointer_cast(&(yHost[0]));
	writeToFile(bx, by, caseFolder, timeStep);
}

/***************************************************************************//**
*/
template <typename memoryType>
void bodies<memoryType>::writeToFile(real *bx, real *by, std::string &caseFolder, int timeStep)
{
	std::string       path;
	std::stringstream out;
	out << caseFolder << '/' << std::setfill('0') << std::setw(7) << timeStep << "/bodies";
	std::ofstream file(out.str().c_str());
	file << '#' << std::setw(19) << "x-coordinate" << std::setw(20) << "y-coordinate" << std::endl;
	for (int l=0; l < totalPoints; l++)
	{
		file << bx[l] << '\t' << by[l] << '\n';
	}
	file.close();
}

template class bodies<host_memory>;
template class bodies<device_memory>;
