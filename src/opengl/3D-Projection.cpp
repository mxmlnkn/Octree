
Matrix CreateRotationMat4f( const Vec<double,3> pu, const double theta )
{
    const int ROT_DIM = 4;
	assert( ROT_DIM >= 3 );
    Matrix m(ROT_DIM,ROT_DIM); m = 0;
    m(0,0) =      0;   m(0,1) = -pu[2];   m(0,2) =  pu[1];
    m(1,0) =  pu[2];   m(1,1) =      0;   m(1,2) = -pu[0];
    m(2,0) = -pu[1];   m(2,1) =  pu[0];   m(2,2) =      0;

    Matrix I(ROT_DIM,ROT_DIM);
    I = 0;
    I.SetDiagonal(1);
    m = cos(theta)*I + sin(theta)*m
        + (1-cos(theta)) * ( Matrix(pu)*Matrix(pu).Transpose() );
    // if ROT_DIM>3 then the matrix should contain a 3x3 matrix in the upper left corner and be zero in the other elements
    for ( int i=4; i < ROT_DIM; i++ )
        m(i,i) = 1;
	return m;
}

Matrix CalcProjection(double fovy, double aspect, double zNear, double zFar)
{	//zNear and zFar always positive! aspect=w/h
	//Note: OpenGL normally views orthogonal from infinity and only shows -w<[x,y,z]<+w but depth buffer from 0 to 1 oO?
	#ifdef DEBUG_SIM
		if(zNear<0 || zFar<0)
			fprintf(stderr, "zNear and zFar aren't positive for this CalcProjecton-Call! zNear:%f, zFar:%f\n",zNear,zFar);
	#endif
	/*****************************************************
	 * Assumptions: -Camera is at(0,0,0) and looks to -z *
	 *  zRar > zNear > 0, w=1                            *
	 * This Matrix does the following to a 4x1 Vector:   *
	 ****** FOVY *****************************************
	 * Wanted: Scale x,y                                 *
	 * f = cos(fovy/2)/sin(fovy/2) = A/O = span_Z/span_Y *
	 * y = f*y; x=f*x; z=z; w=w;                         *
	 * => y is scaled to match the fovy-angle at z=1     *
	 *    (and x is scaled to keep the aspect ratio)     *
	 ****** ASPECT-RATIO: ********************************
	 * (aspect = w/h = span_X/span_Y)                    *
	 * x = x/aspect; y=y; z=z; w=w;                      *
	 * => if w>h then the x-axis will be "elongated"     *
	 ****** Depth Scale **********************************
	 * Wanted: -a 2D object takes the same percent area  *
	 *   on the zNear-Plane as on the zFar-Plane         *
	 *  -for z=pyramide-tip=0 => w=0(scale=inf)          *
	 * height at zF|zN shall be tan theta*zF|zN          *
	 * => h(z)=tan theta*z*h0 => y=y/w; w(z)=1/tan th*z  *
	 * Increase/z_Span = [tan th*zF-tan th*zN]/(zF-zN) = *
	 *  tan theta = f; => w = -z*f;                      *
	 * <=> view_m[3][2]=-f; (other in this row is 0)     *
	 ****** INTERIM RESULTS ******************************
	 * |y|<tan(fovy/2), |x|<f*h/w => -1<(x,y)<+1         *
	 *****************************************************/
	Matrix view_m(4,4);

	double f = 1.0 / tan(0.5*fovy);
	view_m(0,0) = f/aspect;	// x-scale
	view_m(1,1) = f;		// y-scale
	view_m(3,2) = -1;		// w=-z
	/****** CLIPPING: ************************************
	 * Unfortunately we can't really think in scaling    *
	 *  and translating anymore, because of W which      *
	 *  scales our vector by 1/z and is a fucking pain   *
	 *  to calculate. Therefore we simply look at it in  *
	 *  an abstract way, calculating a linear fit with   *
	 *  two points: z'(-zNear)=-1 and z'(-zFar)=+1       *
	 * We have two variables: a:=view_m[2*4+2],          *
	 *  b:=view_m[3*4+3] which gives us z'=1/w*(a*z+b)   *
	 * => z'(-zN)=-1/zN*(-a*zN+b)=a+b/zN=-1  |           *
	 * => z'(-zF)=               =a+b/zF=+1  v -         *
	 *  => 2=b/zF-b/zN <=> b=2*zF*zN/(zN-zF)             *
	 *  => input b into equation 2 => a=(zN-zF-2*zN)/... *
	 *      a=(zF+zN)/(zN-zF)                            *
	 *****************************************************/
	view_m(2,2) = (zFar+zNear)   /(zNear-zFar); //a: zscale
	view_m(2,3) = (2*zFar*zNear) /(zNear-zFar); //b: z-translate (heavily interdependent on w=-z!)

	return view_m;
}

Matrix CalcView(const CameraData & Camera) {
	/*****************************************************
	 * This function rotates and translates the world to *
	 * fit our camera position                           *
	 *****************************************************
	 * new z-Axis = -(Camera.Center-Camera.Eye)          *
	 *  no minus because for Camera.Center=-10,          *
	 *  Camera.Eye=0 => z-axis = (0,0,-1) but in this    *
	 *  case z should be (0,0,1) because we look to its  *
	 *  negative end                                     *
	 * we can't use: new y-Axis = Camera.Up , because    *
	 *  then the user could specify Camera.Up to get a   *
	 *  non-orthonormal system => use only a projection  *
	 *  of Camera.Up. Start by:                          *
	 * new x-Axis = Camera.Up(kinda y-axis) x z_new      *
	 * new y-Axis = z_new cross x_new                    *
	 * => write these into the rows of a matrix to get   *
	 *  the wanted transfomration                        *
	 *     ( x-axis, 0 )   This is the case, because when*
	 *     ( y-axis, 0 )   multiplying it is like SCP    *
	 *  m= ( z-axis, 0 )   for every row in the vector   *
	 *     ( 0,0,0,  1 )                                 *
	 * Now we have the correct rotation, but we need to  *
	 * translate it using Eye. If Eye=(1,1,1) we will    *
	 * translate everything including the camera by      *
	 * (-1,-1,-1) to get our Eye to 0 like anticipated   *
	 * when calculating the perspective and Rotation!    *
	 * Multiplying our matrix to a translation-matrix    *
	 *  does the trick. Only the last column may change, *
	 *  but not easily e.g.: view[4*0+3]=-Camera->Eye[0] *
	 *   *x[0]-Camera->Eye[1]*x[1]-Camera->Eye[2]*x[2]   *
	 *          ( 1,0,0, -Eye[0] )     ( x-axis, T[0] )  *
	 *          ( 0,1,0, -Eye[1] )     ( y-axis, T[1] )  *
	 *  view = m*( 0,0,1, -Eye[2] )  =  ( z-axis, T[2] ) *
	 *          ( 0,0,0,     1   )     ( 0,0,0,   1   )  *
	 * T[0] = <x-Axis,-Eye>; T[1] = <y-Axis,-Eye>;       *
	 *****************************************************/

    Matrix view(4,4);
    Vec<double,3> T,x,y,z;
    z = (-1)*(Camera.Center - Camera.Eye) / (Camera.Center - Camera.Eye).abs();
    x = ( Camera.Up / Camera.Up.norm() ).cross( z );
    y = z.cross( x );

    T[0] = -x.scp(Camera.Eye);
    T[1] = -y.scp(Camera.Eye);
    T[2] = -z.scp(Camera.Eye);

	view(0,0) = x[0];   view(0,1) = x[1];   view(0,2) = x[2];   view(0,3) = T[0];
	view(1,0) = y[0];   view(1,1) = y[1];   view(1,2) = y[2];   view(1,3) = T[1];
	view(2,0) = z[0];   view(2,1) = z[1];   view(2,2) = z[2];   view(2,3) = T[2];
	view(3,0) = 0;      view(3,1) = 0;      view(3,2) = 0;      view(3,3) = 1;
	return view;
}

Matrix CalcModelMatrix
( const Vec<double,3> translationvector, const double w_transl,
  const double obj_scale, const double obj_rot )
{
	Matrix tmp(4,4), model_matrix(4,4);
    tmp.SetDiagonal(1);
	tmp(0,3) = translationvector[0];    // x-translation
	tmp(1,3) = translationvector[1];    // y-translation
	tmp(2,3) = translationvector[2];    // z-translation
	tmp(3,3) = w_transl;                // 1/Scale

    /* more general, because this is only rotation about z-axis !!! */
    model_matrix = tmp * CreateRotationMat4f( Vec<double,3>(0,0,1), obj_rot);
    tmp = 0;
    tmp.SetDiagonal(obj_scale);
	tmp(3,3) = 1; // except w

    return model_matrix * tmp;
}

std::ostream& operator<<( std::ostream& out, const Matrix mat ) {
	out << mat.GetDim().m << "x" << mat.GetDim().n << "-Matrix: {";
	for ( int i=0; i < mat.GetDim().m; i++ ) {
		out << "{";
		for ( int j=0; j < mat.GetDim().n; j++ ) {
			out << mat(i,j);
			if ( j != mat.GetDim().n - 1 )
				out << ",";
		}
		out << "}";
		if ( i != mat.GetDim().m - 1 )
			out << ",";
	}
	out << "}";
	return out;
}
