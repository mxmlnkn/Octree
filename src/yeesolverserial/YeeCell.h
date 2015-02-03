#pragma once

#include <cassert>
#include "math/TVector.h"
#include "math/TVector.tpp"

class YeeCell {
public:
	typedef Vec<double, 3> VecD;

	VecD E;
	VecD H;
	//int internaltimestep;
	double epsilon;  // relativ electric permittivity in this cell
	double mu;       // magnetic permeability in this cell
	double sigmaE;   // electric resistivity in this cell
	double sigmaM;   // magnetic resistivity in this cell

	YeeCell(void) : E(0.0), H(0.0), epsilon(1), mu(1), sigmaE(0), sigmaM(0)	{}
	~YeeCell(void) {}

    /* Methods needed for interpolation (copy and assignment constructor are  *
     * automatically provided and should work. For the array of E and H they  *
     * call the correct assignment operator= of Vec, so no worry there!       */
    YeeCell(double init)
    : E(init), H(init), epsilon(init), mu(init), sigmaE(init), sigmaM(init)
    {}
    YeeCell& operator*=(const YeeCell & a) {
        E       *= a.E;
        H       *= a.E;
        epsilon *= a.epsilon;
        mu      *= a.mu     ;
        sigmaE  *= a.sigmaE ;
        sigmaM  *= a.sigmaM ;
        return *this;
    }
    YeeCell& operator+=(const YeeCell & a) {
        E       += a.E;
        H       += a.E;
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
