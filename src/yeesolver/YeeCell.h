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
	
	YeeCell(void) 
	: epsilon(1), mu(1), sigmaE(0), sigmaM(0)
	{}
	~YeeCell(void) {}
};


/* Enables cout << Cell; This also works with fstream and therefore with tout */
/*ostream& operator<<( ostream& out, const SimulationBox::YeeCell cell ) {
    out << cell.value;
    return out;
}*/
