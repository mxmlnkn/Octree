#pragma once

double t_spawn_func( double t_SI ) {
	double T_SI  = 40e-9 /* m */ / SPEED_OF_LIGHT_SI;
	double sigmaE = T_SI/2.;
	double t0    = 40;
	return sin( 2.*M_PI*t_SI / T_SI );
	return ( t_SI < T_SI ? 1.0 : 0.0 );
	return exp(-pow(t_SI-t0,2)/(2.*sigmaE*sigmaE));
}

namespace TIME_SPAWN_FUNCTIONS {
	double sinewave( double T, double t, double lambda = 1, double x = 0) {
		return std::sin( 2.*M_PI*( x/lambda + t/T ) );
	}
	double sinewave2d( double T, double t, double kx = 0, double x = 0, double ky = 0, double y = 0) {
		return std::sin(  kx*x + ky*y - 2.*M_PI*t/T );
	}
	double PSQ_STEP( double T, double t ) {
		if ( t < 0 )
			return 0;
		else if (t < T/2)
			return 0.5* pow( 1 + (t-T/2.)/(T/2.), 2 );
		else if (t < T)
			return 1 - 0.5* pow( 1 - (t-T/2.)/(T/2.), 2 );
		else
			return 1;
	}
	double gauss( double x, double mu=0, double sigmaE=1 ) {
		return 1./(sigmaE*sqrt(2.*M_PI))*exp(-pow(x-mu,2)/(2.*sigmaE*sigmaE));
	}
	double GaussNotNormed( double x, double mu=0, double sigmaE=1 ) {
		return exp(-pow(x-mu,2)/(2.*sigmaE*sigmaE));
	}
	double GaussianBeam( double z, double r, double t, double w0, double lambda,
                         double T, double E0 = 1 )
    {
        assert( w0 > 2*lambda/M_PI );
        double zR = M_PI*w0*w0/lambda;
        double wz = w0*sqrt( 1 + z*z/(zR*zR) );
        double Rz = z + zR*zR/z; // diverges for z->0 :S?
		return E0*w0/wz * exp(-r*r/(wz*wz)) * cos( 2*M_PI/lambda*( z + r*r/(2*Rz) ) - 2.*M_PI*t/T );
	}
}
