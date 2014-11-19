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
#define TIMESTEPS_NEEDED_FOR_CALCULATION 2


namespace SimulationBox {


bool InArea( const VecI & pos, const int & area, const VecI & localcells, int guardsize );


template<int T_DIM, int T_GUARDSIZE>
class SimulationBox {
    SimulationBox( void );
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    typedef BaseMatrix<CellData,T_DIM> CellMatrix;
    typedef Iterator<T_DIM,T_GUARDSIZE> IteratorType;
    static const int dim       = T_DIM;
    static const int guardsize = T_GUARDSIZE;
    
    /* Current Time on which Calculations are to be done or where done are    *
     * held in t[0]. t[1] is the prior time step and so on                    */
    const int ntimesteps = TIMESTEPS_NEEDED_FOR_CALCULATION+1;
    typedef struct TimeDataStruct {
        BaseMatrix<CellData,T_DIM> cells;
    } TimeData;
    TimeData * t[TIMESTEPS_NEEDED_FOR_CALCULATION+1]; // holds array to different worldmatrices at different times
                  // use array of pointers instead of simple array so I can
                  // cyclically swap the time pointers, instead of the data !

    /* global physics simulationBox properties */
    VecD globsize;
    VecD cellsize;

    /* global cell-algorithm properties */
    VecI globalcells;
    VecI localcells;
    VecD abspos;  // of lower (or better upper ???) left point  (without Guard)
    
    /******************************** Methods *********************************/
    //static SimulationBox& getInstance(void);
    static SimulationBox& getInstance(void) {
        static SimulationBox instance;
        return instance;
    }
    ~SimulationBox(void);
    bool inArea( const VecI & pos, const int & area ) const;
    IteratorType getIterator( const int area ) const;
    void Init( VecD globsize, VecI globalcells, VecI localcells, int mpicoords[T_DIM] );
#if DEBUG_SIMBOX >= 1
    void PrintValues( int timestep = 0 ) const;
#endif

    VecD getGlobalPosition ( const IteratorType it ) const;
    void CopyCurrentToPriorTimestep( void );
    
};

}