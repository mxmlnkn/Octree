#pragma once

#include "SimulationBox.h"


namespace SimulationBox {

/* static members can't be defined in a class declaration, that's why we need *
 * these seemingly useless lines to actually allocate the memory space for    *
 * those variables                                                            */
template<int T_DIMENSION, int T_GUARDSIZE>
const int SimulationBox<T_DIMENSION, T_GUARDSIZE>::dim ;
template<int T_DIMENSION, int T_GUARDSIZE>
const int SimulationBox<T_DIMENSION, T_GUARDSIZE>::guardsize;

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
SimulationBox<T_DIM,T_GUARDSIZE>::SimulationBox( void ) {
    for (int i=0; i<this->ntimesteps; i++)
        this->t[i] = (TimeData*) malloc( sizeof(TimeData) );
};

template<int T_DIM, int T_GUARDSIZE>
SimulationBox<T_DIM,T_GUARDSIZE>::~SimulationBox(void) {
    for (int i=0; i<this->ntimesteps; i++)
        free( this->t[i] );
}

template<int T_DIM, int T_GUARDSIZE>
bool SimulationBox<T_DIM,T_GUARDSIZE>::inArea( const VecI & pos, const int & area ) const {
    return InArea( pos, area, this->localcells, T_GUARDSIZE );
}

template<int T_DIM, int T_GUARDSIZE>
typename SimulationBox<T_DIM,T_GUARDSIZE>::IteratorType SimulationBox<T_DIM,T_GUARDSIZE>::getIterator( const int area ) const {
    return IteratorType( area, localcells + VecI(2*T_GUARDSIZE) );
}

// Later the coordinates shouldn't come from mpi, but from octotree
// should maybe in constructor / getInstance
template<int T_DIM, int T_GUARDSIZE>
void SimulationBox<T_DIM,T_GUARDSIZE>::Init( VecD globsize,
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
    for ( int timestep=0; timestep < this->ntimesteps; timestep++ )
        this->t[timestep]->cells = BaseMatrix<CellData,T_DIM>( localcells + VecI(2*T_GUARDSIZE) );
}

#if DEBUG_SIMBOX >= 1
template<int T_DIM, int T_GUARDSIZE>
void SimulationBox<T_DIM,T_GUARDSIZE>::PrintValues( int timestep ) const {
    tout << "My Raw Matrix with Guard(size=" << guardsize << ") is" << std::endl;
    VecI size = this->t[timestep]->cells.getSize();
    VecI ind(0);
    for ( ind[Y_AXIS]=0; ind[Y_AXIS]<size[Y_AXIS]; ind[Y_AXIS]++) {
        for ( ind[X_AXIS]=0; ind[X_AXIS]<size[X_AXIS]; ind[X_AXIS]++)
            tout << this->t[timestep]->cells[ ind ].value << " ";
        tout << std::endl;
    }
    tout << std::endl << std::flush;

    /* Wait a bit till everything is flushed out, to not get scrambled output */
    double t0 = MPI_Wtime();
    while( MPI_Wtime() - t0 < 0.1*( double(rand())/double(RAND_MAX) ) ) {}
}
#endif

template<int T_DIM, int T_GUARDSIZE>
typename SimulationBox<T_DIM, T_GUARDSIZE>::VecD SimulationBox<T_DIM,T_GUARDSIZE>::getGlobalPosition ( const IteratorType it ) const {
    return abspos + ( VecD(it.icell) * cellsize );
}

template<int T_DIM, int T_GUARDSIZE>
void SimulationBox<T_DIM,T_GUARDSIZE>::CopyCurrentToPriorTimestep( void ) {
    // t[ntimesteps-1] is the oldest position => shift t[i] to t[i+1]
    TimeData * tmp = t[ntimesteps-1];
    for (int j = ntimesteps-1; j > 0; j--) {
        // this only swaps pointers ( but it also works with copy 
        // assignment operators defined in the data structure )
        t[j] = t[j-1];
        t[j] = t[j-1];
    }
    for ( int pos=0; pos < this->t[0]->cells.getSize().product(); pos++ )
        tmp->cells[pos] = t[0]->cells[pos];
    t[0] = tmp;
}


} // namespace SimulationBox
