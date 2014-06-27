/***************************************************************************//**
* \file body.h
* \author Krishnan, A. (anus@bu.edu)
* \brief Declaration and definition of the class \c body
*/

#pragma once

/***************************************************************************//**
* \class body
* \brief Information about an immersed body in the flow
*
* The number of boundary points and the reference x- and y-coordinates of the
* boundary points are provided. The motion of the body is also described.
* For uniform motion, the translational velocity needs to be specified, and for
* uniform rotation, the angular velocity has to be specified.
*/
class body
{
  public:
    int  numPoints; ///< number of boundary points
	
	vecH X,         ///< Reference x-coordinate of boundary points
	     Y;         ///< Reference y-coordinate of boundary points
         
	real X0[2];     ///< Reference centre of rotation
	     
	real Xc0[2],    ///< Initial position of center of rotation
	     Theta0;    ///< Initial angle of attack
    
    real Xc[2],     ///< Actual centre of rotation (x- and y-coordinates)
	     Theta,     ///< Actual angle of attack (counterclockwise is positive)
	     vel[2],    ///< Translational velocity (x- and y- components)
	     angVel;    ///< Angular velocity (counterlockwise is positive)

	bool moving[2]; ///< Flag to indicate if the body is moving (translating or rotating)

	real velocity[2], ///< Uniform translational velocity (x- and y-components)
	     omega;       ///< Uniform rotational velocity
         
	real xOscillation[3],     ///< Amplitude, angular frequency and phase difference of oscillation in the x-direction
         yOscillation[3],     ///< Amplitude, angular frequency and phase difference of oscillation in the y-direction
         pitchOscillation[3]; ///< Amplitude, angular frequency and phase difference of pitch oscillation
    
	/********************//**
	* \brief Update the position of the center of rotation of the body and
	*        the translational and rotational velocities of the boundary points.
	*
	* Requires a detailed description.
	*/
	void update(real Time)
	{
		// If the body is translating
		if(moving[0])
		{
			Xc[0] = Xc0[0] + velocity[0]*Time + xOscillation[0]*sin(xOscillation[1]*Time + xOscillation[2]);
			Xc[1] = Xc0[1] + velocity[1]*Time + yOscillation[0]*sin(yOscillation[1]*Time + yOscillation[2]);
			vel[0]= velocity[0] + xOscillation[0]*xOscillation[1]*cos(xOscillation[1]*Time + xOscillation[2]);
			vel[1]= velocity[1] + yOscillation[0]*yOscillation[1]*cos(yOscillation[1]*Time + yOscillation[2]);
		}
		// If the body is rotating
		if(moving[1])
		{
			Theta = Theta0 + omega*Time + pitchOscillation[0]*sin(pitchOscillation[1]*Time + pitchOscillation[2]);
			angVel = omega + pitchOscillation[0]*pitchOscillation[1]*cos(pitchOscillation[1]*Time + pitchOscillation[2]);
		}
	}
};
