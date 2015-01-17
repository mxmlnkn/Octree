#pragma once

#include "SimulationBoxDefines.h"
#include "SimulationBox.h"
#include "math/TBaseMatrix.h"
#include "math/TVector.h"

namespace SimulationBox {

template<int T_DIM, typename T_CELLDATA>
class Iterator {
public:
    /* Definitions and Aliases */
    const int dim = T_DIM;
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    typedef BaseMatrix<T_CELLDATA,T_DIM> CellMatrix;

private:
    int  area;   // next 3 are derived from this one
    bool core;
    bool border;
    bool guard;
    CellMatrix * srcmat;

public:
    /* Only public for convenience, shouldn't be tempered with !!! */
    VecI icell;  // stores current cell index. This is where we are
    VecI ncells; // dimension of cellmatrix this iterator works on with Guard
    int guardsize;

    Iterator( const int area, const VecI ncells, const int guardsize, CellMatrix & srcmat );
    ~Iterator( void );
    Iterator( const Iterator & src );
    int getCellIndex( void ) const;
    Iterator begin( void ) const;
    Iterator end( void ) const;
    Iterator& operator<<( const int index[T_DIM] );
    Iterator& operator++( void );
    bool operator==( const Iterator & it );
    bool operator!=( const Iterator & it );
    void operator=( const Iterator & src );
    T_CELLDATA& operator*(void) const;
    T_CELLDATA* operator->(void) const;

};

}
