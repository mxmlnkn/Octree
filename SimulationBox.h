#pragma once

#include <cstdlib>  // malloc
#include <iostream>  // malloc
#include "Vector.h"
#include "BaseMatrix.h"
#include "TeeStream.h"
#include "Cell.h"
#include "SimulationBoxDefines.h"
#include "SimulationBoxIterator.h"


#define DEBUG_SIMBOX 1


namespace SimulationBox {


bool InArea( const VecI & pos, const int & area, const VecI & localcells, int guardsize );


template<int T_DIM, int T_GUARDSIZE>
class SimulationBox {
    SimulationBox( void );
public:
    typedef BaseMatrix<CellData,T_DIM> CellMatrix;
    typedef Iterator<T_DIM,T_GUARDSIZE> IteratorType;

    BaseMatrix<CellData,T_DIM> cells;
    static const int dim       = T_DIM;
    static const int guardsize = T_GUARDSIZE;
    
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;

    //static SimulationBox& getInstance(void);
    static SimulationBox& getInstance(void) {
        static SimulationBox instance;
        return instance;
    }
    ~SimulationBox(void);
    bool inArea( const VecI & pos, const int & area ) const;
    IteratorType getIterator( const int area ) const;
    CellData & operator[]( const IteratorType & it );
    void Init( VecD globsize,
               VecI globalcells,
               VecI localcells,
               int mpicoords[T_DIM] );
#if DEBUG_SIMBOX >= 1
    void PrintValues( void ) const;
#endif

    /* global physics simulationBox properties */
    VecD globsize;
    VecD cellsize;

    /* global cell-algorithm properties */
    VecI globalcells;
    VecI localcells;
    VecD abspos;  // of lower (or better upper ???) left point  (without Guard)

    VecD getGlobalPosition ( const IteratorType it ) const;

};

}