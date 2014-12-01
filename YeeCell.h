#pragma once

#include <cassert>
#include "Vector.h"
#include "Vector.tpp"

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
	double sigma;    // electric resistivity in this cell
	double rhoprime; // magnetic resistivity in this cell
	
	YeeCell(void) 
	: epsilon(1), mu(1), sigma(0), rhoprime(0)
	{}
	~YeeCell(void) {}
	/*VecD & E( int time ) {
		assert( time <= 1 );
		assert( time > 1-YEE_CELL_TIMESTEPS_TO_SAVE );
		return Edata[ (time + internaltimestep) % YEE_CELL_TIMESTEPS_TO_SAVE ];
	}
	VecD & H( int time ) {
		assert( time <= 1 );
		assert( time > 1-YEE_CELL_TIMESTEPS_TO_SAVE );
		return Hdata[ (time + internaltimestep) % YEE_CELL_TIMESTEPS_TO_SAVE ];
	}*/
	/* mod not really necessary, but if we have more timesteps than int   *
	 * can specify, than an overflow could occur and the timesteps could  *
	 * get scrambled, therefore better use modulo after incrementing      */
	/*void incTimestep( void ) {
		internaltimestep = ++internaltimestep % YEE_CELL_TIMESTEPS_TO_SAVE;
	}*/
	/* Problem with this is that it would be stored for every cell and therefore wouldn't be too much redundant -.- */
};


/* Enables cout << Cell; This also works with fstream and therefore with tout */
/*ostream& operator<<( ostream& out, const SimulationBox::YeeCell cell ) {
    out << cell.value;
    return out;
}*/
