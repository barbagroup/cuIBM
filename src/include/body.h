/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
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

#pragma once
/**
* \brief Stores the locations of the body points and the motion parameters
*/
class body
{
  public:
    int  numPoints; ///< number of boundary points

    real X0[2], ///< Reference centre of rotation
         Theta0; ///< Reference angle
    
    real Xc[2],
         Theta,
         vel[2],
         angVel;

    vecH X, Y; ///< Reference location of boundary points

    bool moving[2];

    real velocity[2], omega;
    real xOscillation[3],
         yOscillation[3],
         pitchOscillation[3];
         
	void update(real Time)
	{
		/// If the body is translating
		//if(moving[0])
		//{
			Xc[0] = X0[0] + velocity[0]*Time + xOscillation[0]*sin(xOscillation[1]*Time + xOscillation[2]);
			Xc[1] = X0[1] + velocity[1]*Time + yOscillation[0]*sin(yOscillation[1]*Time + yOscillation[2]);
			vel[0]= velocity[0] + xOscillation[0]*xOscillation[1]*cos(xOscillation[1]*Time + xOscillation[2]);
			vel[1]= velocity[1] + yOscillation[0]*yOscillation[1]*cos(yOscillation[1]*Time + yOscillation[2]);
		//}
		/// If the body is rotating
		//if(moving[1])
		//{
			Theta = Theta0 + omega*Time + pitchOscillation[0]*sin(pitchOscillation[1]*Time + pitchOscillation[2]);
			angVel = omega + pitchOscillation[0]*pitchOscillation[1]*cos(pitchOscillation[1]*Time + pitchOscillation[2]);
		//}
	}
};
