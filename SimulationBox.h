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

bool InArea( const VecI & pos, const int & area, const VecI & localcells, int guardsize ) {
    assert( area != 0 ); // would be trivial iterator
    assert( area <= CORE+BORDER+GUARD );

    bool inCore   = (     ( pos >= VecI(2*guardsize) )
                      and ( pos <  VecI(guardsize) + localcells - VecI(guardsize) )
                    );
    bool inBorder = (     ( pos >= VecI(guardsize) )
                      and ( pos <  VecI(guardsize) + localcells )
                      and ( not inCore ) 
                    );
    bool inGuard  = (     ( pos >= VecI(0)) 
                      and ( pos <  VecI(guardsize) + localcells + VecI(guardsize) )
                      and ( not inBorder ) 
                      and ( not inCore   )
                    );

    bool allor = false;
    if ( area & GUARD )
        allor = allor or inGuard;
    if ( area & BORDER )
        allor = allor or inBorder;
    if ( area & CORE )
        allor = allor or inCore;

    return allor;
}


template<int T_DIM, int T_GUARDSIZE>
class SimulationBox {
    SimulationBox( void ) {};
public:
    typedef BaseMatrix<CellData,T_DIM> CellMatrix;
    typedef Iterator<T_DIM,T_GUARDSIZE> IteratorType;

    BaseMatrix<CellData,T_DIM> cells;
    static const int dim       = T_DIM;
    static const int guardsize = T_GUARDSIZE;
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;

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
    
    bool inArea( const VecI & pos, const int & area ) const {
        return InArea( pos, area, this->localcells, T_GUARDSIZE );
    }
/*
    bool inArea( const VecI & pos, const int & area ) const {
        assert( area != 0 ); // would be trivial iterator
        assert( area <= CORE+BORDER+GUARD );

        bool inCore   = (     ( pos >= VecI(2*T_GUARDSIZE) )
                          and ( pos <  VecI(T_GUARDSIZE) + localcells - VecI(T_GUARDSIZE) )
                        );
        bool inBorder = (     ( pos >= VecI(T_GUARDSIZE) )
                          and ( pos <  VecI(T_GUARDSIZE) + localcells )
                          and ( not inCore ) 
                        );
        bool inGuard  = (     ( pos >= VecI(0)) 
                          and ( pos <  VecI(T_GUARDSIZE) + localcells + VecI(T_GUARDSIZE) )
                          and ( not inBorder ) 
                          and ( not inCore   )
                        );

        bool allor = false;
        if ( area & GUARD )
            allor = allor or inGuard;
        if ( area & BORDER )
            allor = allor or inBorder;
        if ( area & CORE )
            allor = allor or inCore;

        #if DEBUG_SIMBOX >= 2
            terr << std::endl << "[InArea: (" ;
            for ( int i=0; i<T_DIM-1; i++ )
                terr << pos[i] << ",";
            terr << pos[T_DIM-1] << ") ";
            terr << "is in ";
            if (inCore  ) terr << "CORE ";
            if (inBorder) terr << "BORDER ";
            if (inGuard ) terr << "GUARD ";
            if (!inCore and !inBorder and !inGuard) terr << "None ";
            terr << ", localcells=(" ;
            for ( int i=0; i<T_DIM-1; i++ )
                terr << localcells[i] << ",";
            terr << localcells[T_DIM-1] << ") ";
            terr << " GuardSize=" << guardsize;
            terr << "]" << std::endl;
        #endif
        
        return allor;
    }
*/

    IteratorType getIterator( const int area ) const {
        return IteratorType( area, localcells + VecI(2*T_GUARDSIZE), this );
    }

    CellData & operator[]( const IteratorType & it ) {
        return cells[ it.icell ];
    }

    // Later the coordinates shouldn't come from mpi, but from octotree
    // should maybe in constructor / getInstance
    void Init( VecD globsize,
               VecI globalcells,
               VecI localcells,
               int mpicoords[T_DIM] )
    {
        this->globsize    = globsize;
        this->globalcells = globalcells;
        this->localcells  = localcells;
        this->cellsize    = globsize / VecD(globalcells);
        VecD coords( mpicoords );
        // Assuming localcells is the same for every process !!!
        this->abspos = coords * localcells / VecD( globalcells ) * globsize;
        this->cells = BaseMatrix<CellData,T_DIM>( localcells + VecI(2*T_GUARDSIZE) );
    }

#if DEBUG_SIMBOX >= 1
    void PrintValues( void ) {
        tout << "My Raw Matrix with Guard(size=" << guardsize << ") is" << std::endl;
        VecI size = this->cells.getSize();
        VecI ind(0);
        for ( ind[Y_AXIS]=0; ind[Y_AXIS]<size[Y_AXIS]; ind[Y_AXIS]++) {
            for ( ind[X_AXIS]=0; ind[X_AXIS]<size[X_AXIS]; ind[X_AXIS]++)
                tout << cells[ ind ].value << " ";
            tout << std::endl;
        }
        tout << std::endl << std::flush;
        
        /* Wait a bit till everything is flushed out, to not get scrambled output */
        double t0 = MPI_Wtime();
        while( MPI_Wtime() - t0 < 0.1*( double(rand())/double(RAND_MAX) ) ) {}
    }
#endif

    /* global physics simulationBox properties */
    VecD globsize;
    VecD cellsize;

    /* global cell-algorithm properties */
    VecI globalcells;
    VecI localcells;
    VecD abspos;  // of lower (or better upper ???) left point  (without Guard)


    VecD getGlobalPosition ( Iterator<T_DIM,T_GUARDSIZE> it ) {
        return abspos + ( VecD(it.icell) * cellsize );
    }

};

}