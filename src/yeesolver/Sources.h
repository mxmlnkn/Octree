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
}
