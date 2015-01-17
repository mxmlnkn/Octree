#pragma once

#include "SimulationBox.h"

namespace SimulationBox {

template<int T_DIM>
bool InArea( const Vec<int,T_DIM> & pos, const int & area, const Vec<int,T_DIM> & localcells, int guardsize ) {
    typedef Vec<int,T_DIM> VecI;
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

/********************************* Constructor ********************************/
template<int T_DIM, typename T_CELLDATA>
SimulationBox<T_DIM,T_CELLDATA>::SimulationBox(
    VecD abspos, VecI localcells, VecD cellsize, int guardsize, int bufferpages
) : ntimesteps(bufferpages), abspos(abspos), localcells(localcells), cellsize(cellsize),
    guardsize(guardsize)
{
    t = (TimeData**) malloc( ntimesteps * sizeof(TimeData*) );
    for (int i=0; i<this->ntimesteps; i++) {
        t[i] = new TimeData;
        t[i]->cells = CellMatrix( localcells + VecI(2*this->guardsize) );
    }
}

/********************************* Destructor *********************************/
template<int T_DIM, typename T_CELLDATA>
SimulationBox<T_DIM,T_CELLDATA>::~SimulationBox(void) {
    if ( t != NULL ) {
        for (int i=0; i<this->ntimesteps; i++)
            if ( this->t[i] != NULL )
                delete this->t[i];
        free(t);
    }
}

template<int T_DIM, typename T_CELLDATA>
bool SimulationBox<T_DIM,T_CELLDATA>::inArea( const VecI & pos, const int & area ) const {
    return InArea( pos, area, this->localcells, this->guardsize );
}

/******************************** getIterator *********************************/
template<int T_DIM, typename T_CELLDATA>
typename SimulationBox<T_DIM,T_CELLDATA>::IteratorType SimulationBox<T_DIM,T_CELLDATA>::getIterator( const int area ) const {
    return IteratorType( area, localcells + VecI(2*guardsize), guardsize );
}

#if DEBUG_SIMBOX >= 1
template<int T_DIM, typename T_CELLDATA>
void SimulationBox<T_DIM,T_CELLDATA>::PrintValues( int timestep ) const {
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

/****************************** getGlobalPosition *****************************/
template<int T_DIM, typename T_CELLDATA>
typename SimulationBox<T_DIM,T_CELLDATA>::VecD SimulationBox<T_DIM,T_CELLDATA>::getGlobalPosition ( const IteratorType it ) const {
    return abspos + VecD(it.icell)*cellsize;
}

/************************* copyCurrentToPriorTimestep ************************/
template<int T_DIM, typename T_CELLDATA>
void SimulationBox<T_DIM,T_CELLDATA>::copyCurrentToPriorTimestep( void ) {
    /* t[ntimesteps-1] is the oldest position => shift t[i] to t[i+1] */
    TimeData * tmp = t[ntimesteps-1];
    for (int j = ntimesteps-1; j > 0; j--)
        t[j] = t[j-1];
    t[0] = tmp;
    /* copy all matrix elements to obsolete t[0] data */
    for ( int pos=0; pos < this->t[1]->cells.getSize().product(); pos++ )
        t[0]->cells[pos] = t[1]->cells[pos];
}

} // namespace SimulationBox
