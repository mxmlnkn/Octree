#pragma once

#include <cassert>
#include "math/TVector.h"
#include "math/TVector.tpp"

#ifndef YEE_CELL_TIMESTEPS_TO_SAVE
	#define YEE_CELL_TIMESTEPS_TO_SAVE 3
#endif

class YeeCell {
public:
	typedef Vec<double, 3> VecD;

	VecD E[YEE_CELL_TIMESTEPS_TO_SAVE];
	VecD H[YEE_CELL_TIMESTEPS_TO_SAVE];
	//int internaltimestep;
	double epsilon;  // relativ electric permittivity in this cell
	double mu;       // magnetic permeability in this cell
	double sigmaE;    // electric resistivity in this cell
	double sigmaM; // magnetic resistivity in this cell

	YeeCell(void) : epsilon(1), mu(1), sigmaE(0), sigmaM(0)	{}
	~YeeCell(void) {}

    /* Methods needed for interpolation (copy and assignment constructor are  *
     * automatically provided and should work. For the array of E and H they  *
     * call the correct assignment operator= of Vec, so no worry there!       */
    YeeCell(double init)
    : epsilon(init), mu(init), sigmaE(init), sigmaM(init) {
        for ( int i=0; i<YEE_CELL_TIMESTEPS_TO_SAVE; i++ ) {
            E[i] = VecD(init);
            H[i] = VecD(init);
        }
    }
    YeeCell& operator*=(const YeeCell & a) {
        for ( int i=0; i<YEE_CELL_TIMESTEPS_TO_SAVE; i++ ) {
            E[i] *= a.E[i];
            H[i] *= a.E[i];
        }
        epsilon *= a.epsilon;
        mu      *= a.mu     ;
        sigmaE  *= a.sigmaE ;
        sigmaM  *= a.sigmaM ;
        return *this;
    }
    YeeCell& operator+=(const YeeCell & a) {
        for ( int i=0; i<YEE_CELL_TIMESTEPS_TO_SAVE; i++ ) {
            E[i] += a.E[i];
            H[i] += a.E[i];
        }
        epsilon += a.epsilon;
        mu      += a.mu     ;
        sigmaE  += a.sigmaE ;
        sigmaM  += a.sigmaM ;
        return *this;
    }
};


/* Enables cout << Cell; This also works with fstream and therefore with tout */
/*ostream& operator<<( ostream& out, const SimulationBox::YeeCell cell ) {
    out << cell.value;
    return out;
}*/
