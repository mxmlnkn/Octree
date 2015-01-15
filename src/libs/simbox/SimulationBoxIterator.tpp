#pragma once

#include "SimulationBoxIterator.h"

namespace SimulationBox {

/******************************** Constructor *********************************/
template<int T_DIM>
Iterator<T_DIM>::Iterator( const int area, const VecI ncells, const int guardsize )
: guardsize(guardsize), area(area), ncells(ncells) {
    assert( area != 0 ); // would be trivial iterator
    assert( area <= CORE+BORDER+GUARD );
    assert( area != CORE+GUARD ); // more difficult to implement,
                                  // and makes no sense
    this->core   = area & CORE;
    this->border = area & BORDER;
    this->guard  = area & GUARD;
    icell = 0;
}

/******************************** Denstructor *********************************/
template<int T_DIM>
Iterator<T_DIM>::~Iterator( void ) {};

/****************************** Copy Constructor ******************************/
template<int T_DIM>
Iterator<T_DIM>::Iterator( const Iterator & src ) {
    this->core   = src.core  ;
    this->border = src.border;
    this->guard  = src.guard ;
    this->area   = src.area ;
    this->icell  = src.icell ;
    this->ncells = src.ncells;
}

/**********************************************************************
 * If we have 2 slates of 3x4 length and we wanted the index          *
 * [i,j,k] = [2,3,2] , then the matrix lies in the memory like:       *
 *   [oooo|oooo|oooo] [oooo|oooo|oxoo]                                *
 * This means to address that element (above marked x) we need to     *
 * calculate:                                                         *
 *   (i=2)*[ (nj=3) * (nk=4) ] + (j=3)*[ (nk=4) ] + (k=2)*[ 1 ]       *
 * That argument in [] will be named 'prevrange'                      *
 **********************************************************************/
template<int T_DIM>
int Iterator<T_DIM>::getCellIndex( void ) const {
    int index     = 0;
    int prevrange = 1;
    for (int i=T_DIM-1; i>=0; i--) {
        index     += icell[i] * prevrange;
        prevrange *= ncells[i];
    }
    return index;
}

template<int T_DIM>
Iterator<T_DIM> Iterator<T_DIM>::begin( void ) const {
    Iterator tmp( *this );

    int upperleftcorner;
    if (guard)
        upperleftcorner = 0;
    else if (border)
        upperleftcorner = guardsize; // not guardsize-1 !
    else if (core)
        upperleftcorner = 2*guardsize;
    // Above is equivalent to cryptic: upperleftcorner = area / 2;

    tmp.icell = VecI( upperleftcorner );
    return tmp;
}

/* return VecI(-1) which is a end marker. Will have to set it to that *
 * in operator++ for example                                          */
template<int T_DIM>
Iterator<T_DIM> Iterator<T_DIM>::end( void ) const {
    Iterator tmp( *this );
    tmp.icell  = VecI(-1);
    return tmp;
}

/* Set the iterator to some index */
template<int T_DIM>
Iterator<T_DIM>& Iterator<T_DIM>::operator<<( const int index[T_DIM] ) {
    VecI ind = VecI(index);
    assert( InArea( ind, area, (ncells-VecI(2*guardsize)), guardsize ) );
    *this = this->begin();
    this->icell += ind;
    return *this;
}

template<int T_DIM>
Iterator<T_DIM> & Iterator<T_DIM>::operator++( void ) { // only prefix, because it's faster!
    do {
        int linindex = ConvertVectorToLinearIndex( icell, ncells ) + 1;
        if ( !( linindex < ncells.product() ) ) {
            this->icell = VecI(-1);
            break;
        }
        this->icell = ConvertLinearToVectorIndex( linindex, ncells );
    } while( !InArea( this->icell, area, (ncells-VecI(2*guardsize)), guardsize ) );
    return *this;
}

template<int T_DIM>
bool Iterator<T_DIM>::operator==( const Iterator & it ) {
    bool alland = true;
    alland = alland & ( this->core   == it.core  );
    alland = alland & ( this->border == it.border );
    alland = alland & ( this->guard  == it.guard  );
    alland = alland & ( this->area   == it.area  );

    alland = alland & ( this->icell  == it.icell  );
    alland = alland & ( this->ncells == it.ncells );

    return alland;
}

template<int T_DIM>
bool Iterator<T_DIM>::operator!=( const Iterator & it ) {
    return !( (*this) == it );
}

template<int T_DIM>
void Iterator<T_DIM>::operator=( const Iterator & src ) {
    this->core   = src.core  ;
    this->border = src.border;
    this->guard  = src.guard ;
    this->area   = src.area ;
    this->icell  = src.icell ;
    this->ncells = src.ncells;
}

}