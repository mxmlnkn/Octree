#pragma once

#include <cstdlib>  // malloc
#include <iostream>  // malloc
#include "Vector.h"
#include "BaseMatrix.h"

using namespace std;

#define DEBUG_SIMBOX 1


typedef struct CellDataStruct {
    int value;
} CellData;

template<int T_DIMENSION, int T_GUARDSIZE>
class SimulationBox {
    SimulationBox( void ) {};


public:
    typedef BaseMatrix<CellData,T_DIMENSION> CellMatrix;

    BaseMatrix<CellData,T_DIMENSION> cells;
    static const int dim       = T_DIMENSION;
    static const int guardsize = T_GUARDSIZE;
    typedef Vec<double,T_DIMENSION> VecD;
    typedef Vec<int   ,T_DIMENSION> VecI;

    static SimulationBox& getInstance(void) {
        static SimulationBox instance; // Guaranteed to be destroyed.
        return instance;
    }
    ~SimulationBox(void) {
    }


    /**************************************************************************
     * The 'Halo' in Musubi is here being called Guard for consistency with   *
     * PicOnGPU. Musbu doesn't differ between Core and Border. In PicOnGPU    *
     * these two are differentiated in order to do calculations on the Core,  *
     * while sending the data from the Border to the respective Guard of      *
     * processes.                                                             *
     * The Guard must be as wide as we need neighbor cells for one calcula-   *
     * tion in one cell. The Border is as wide as the Guard by definition.    *
     * Therefore cells in the Core can be calculated by only including the    *
     * annexing Core- or Bordercells!                                         *
     *                                             Rank 1                     *
     *                                        -----------------               *
     *                                       | G   G G G G   G |              *
     *                   Rank 0              |    ---------    |              *
     *              -----------------        | G | B B B B | G |              *
     *             | G   G G G G   G |       | G | B C C B | G |              *
     *             |    ---------    |       | G | B C C B | G |              *
     *             | G | B B B B | G |       | G | B B B B | G |              *
     *             | G | B C C B | G |       |    ---------    |              *
     *             | G | B C C B | G |       | G   G G G G   G |              *
     *             | G | B B B B | G |        -----------------               *
     *             |    ---------    |                                        *
     *             | G   G G G G   G |                                        *
     *              -----------------                                         *
     * In Order to do calculations on these areas easily, iterators and       *
     * relative cell-directions need to be included                           *
     * The internal data matrix will begin with the upper left corner like    *
     * any normal matrix                                                      *
     **************************************************************************/

    /* This are predefined constant for addressing. E.g. if we want to        *
     * iterate only over the CORE or we want a value on the X_AXIS. These are *
     * are not sizes or widths!                                               */
    static const int GUARD  = 1;
    static const int BORDER = 2;
    static const int CORE   = 4;
    
    static const int X_AXIS = 0;
    static const int Y_AXIS = 1;
    static const int Z_AXIS = 2;

    #include "SimulationBoxIterator.h"
    
    Iterator getIterator( const int area ) const {
        return Iterator( area, localcells + VecI(2*T_GUARDSIZE) );
    }

    CellData & operator[]( const Iterator & it ) {
        return cells[ it.icell ];
    }

    // Later the coordinates shouldn't come from mpi, but from octotree
    // should maybe in constructor / getInstance
    void Init( VecD globsize,
               VecI globalcells,
               VecI localcells,
               int mpicoords[T_DIMENSION] )
    {
        this->globsize    = globsize;
        this->globalcells = globalcells;
        this->localcells  = localcells;
        this->cellsize    = globsize / VecD(globalcells);
        VecD coords( mpicoords );
        // Assuming localcells is the same for every process !!!
        this->abspos = coords * localcells / VecD( globalcells ) * globsize;
        this->cells = BaseMatrix<CellData,T_DIMENSION>( localcells + VecI(2*T_GUARDSIZE) );
    }
    
#if DEBUG_SIMBOX >= 1
    void PrintValues( void ) {
        cout << "My Raw Matrix with Guard(size=" << guardsize << ") is" << endl;
        VecI size = this->cells.size;
        VecI ind(0);
        for ( ind[Y_AXIS]=0; ind[Y_AXIS]<size[Y_AXIS]; ind[Y_AXIS]++) {
            for ( ind[X_AXIS]=0; ind[X_AXIS]<size[X_AXIS]; ind[X_AXIS]++)
                cout << cells[ ind ].value << " ";
            cout << endl;
        }
        cout << endl;
    }
#endif

    /* global physics simulationBox properties */
    VecD globsize;
    VecD cellsize;

    /* global cell-algorithm properties */
    VecI globalcells;
    VecI localcells;
    VecD abspos;  // of lower (or better upper ???) left point  (without Guard)


    VecD getGlobalPosition ( Iterator it ) {
        return abspos + ( VecD(it.icell) * cellsize );
    }
    //delinCoreCoord

};
