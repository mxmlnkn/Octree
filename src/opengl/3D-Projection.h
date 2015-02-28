#pragma once

#include <cmath>
#include "math/Matrix.h"
#include "math/TVector.h"

struct CameraData {
    typedef Vec<double,3> VecD;
	VecD Eye;       // Position of the Camera (initialize: {0,0,0})
	VecD Center;    // Point where we look at (initialize: {0,0,-1})
    /* specifiying what up is (relative!!!) to be able to tile the camera *
     * along the viewing-axis (KEY_LEFT, KEY_RIGHT) (initialize: {0,1,0}) */
	VecD Up;
    /* spherial coordinates with theta ranging from -pi/2 (-z-axis) to    *
     * +pi/2 (+z-axis) (initialize: theta=-M_PI, PSI=0)                   */
	double Phi,Theta,Psi;
	/* the initial values let the camera look to the -z-axis                  *
	 * => Center = Eye + [Cos Phi Sin Theta, Sin Phi Sin Theta, Cos Theta]    *
	 *    Psi is the angle from the y-axis and the vector specified by Theta  *
     *    and Phi ... It's to tilt the Camera.Up                              *
	 * => Camera.Up = RotateMatrixAround(Camera.Center-Camera.Eye, Degree)*   *
     *    (0,1,0)                                                             */
     /*           y                                                   *
      *           ^                                                   *
      *           |                                                   *
      *      Up   |         /                                         *
      *      ^    |       /                                           *
      *      |    |     /      o                                      *
      *      |    |   /      Center                                   *
      *      o    | /                                                 *
      *     Eye   +--------------------------> x                      *
      *         / |                                                   *
      *       /   |                                                   *
      *     /     |                                                   *
      *    z                                                          */
     /* Is the scaling of Up relevant ??? */
    CameraData( VecD pEye, VecD pCenter, VecD pUp = VecD(0,1,0) )
    : Eye(pEye), Center(pCenter), Up(pUp), Phi(0),
      Theta( atan( (Center[1]-Eye[1]) / fabs(Center[2]-Eye[2]) )), Psi(0)
    {}
};

Matrix CalcProjection(double fovy, double aspect, double zNear, double zFar);

Matrix CreateRotationMat4f(const float ux, const float uy, const float uz, const float theta, float** m);
Matrix CalcView(const CameraData & Camera);
Matrix CalcModelMatrix( const Vec<double,3> translationvector,
    const double w_transl, const double obj_scale, const double obj_rot );

std::ostream& operator<<( std::ostream& out, const Matrix mat );

#include "3D-Projection.cpp"
