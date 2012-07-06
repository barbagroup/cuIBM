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

    vecH X, Y; ///< Reference location of boundary points

    bool moving[2];

    real velocity[2], omega;
    real xOscillation[3],
         yOscillation[3],
         pitchOscillation[3];
};
