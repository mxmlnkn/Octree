#pragma once

#include "SimulationBoxDefines.h"
#include "SimulationBox.h"

namespace SimulationBox {

template<int T_DIM, int T_GUARDSIZE>
class Iterator {
private:
    const void * simbox;
public:
    /* Only public for convenience, shouldn't be tempered with !!! */
    bool core;
    bool border;
    bool guard;

    int area;

    VecI icell;  // stores current cell index. This is where we are
    VecI ncells; // with Guard !!!

    Iterator( const int area, const VecI ncells, const void* simbox );
    ~Iterator( void );
    Iterator( const Iterator & src );
    int getCellIndex( void ) const;
    Iterator begin( void ) const;
    Iterator end( void ) const;
    Iterator& operator<<( const int index[T_DIM] );
    Iterator & operator++( void );
    bool operator==( const Iterator & it );
    bool operator!=( const Iterator & it );
    void operator=( const Iterator & src );

};

}