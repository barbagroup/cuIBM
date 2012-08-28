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

/**
* @file body.h
* @brief Stores the locations of the body points and the motion parameters
*/

#pragma once

/**
* @class body
* @brief Information about an immersed body in the flow
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

    real X0[2],     ///< Reference centre of rotation
         Theta0;    ///< Reference angle of attack
    
    real Xc[2],     ///< Actual centre of rotation (x- and y-coordinates)
         Theta,     ///< Actual angle of attack (counterclockwise is positive)
         vel[2],    ///< Uniform translational velocity (x- and y- components)
         angVel;    ///< Uniform angular velocity (counterlockwise is positive)

    vecH X,         ///< Reference x-coordinate of boundary points
         Y;         ///< Reference y-coordinate of boundary points

    bool moving[2]; ///< Flag to indicate if the body is moving (translating or rotating)

    real velocity[2], ///< uniform translational velocity (x- and y-components)
         omega;       ///< uniform rotational velocity
         
    real xOscillation[3],     ///< amplitude, angular frequency and phase difference of oscillation in the x-direction
         yOscillation[3],     ///< amplitude, angular frequency and phase difference of oscillation in the y-direction
         pitchOscillation[3]; ///< amplitude, angular frequency and phase difference of pitch oscillation
    
    /**
    * @brief Update the position of the center of rotation of the body and
    *        the translational and rotational velocities of the boundary points
    *
    * Requires a detailed description.
    */
	void update(real Time)
	{
		// If the body is translating
		if(moving[0])
		{
			Xc[0] = X0[0] + velocity[0]*Time + xOscillation[0]*sin(xOscillation[1]*Time + xOscillation[2]);
			Xc[1] = X0[1] + velocity[1]*Time + yOscillation[0]*sin(yOscillation[1]*Time + yOscillation[2]);
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
