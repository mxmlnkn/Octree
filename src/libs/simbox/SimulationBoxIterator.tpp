#pragma once

#include "SimulationBoxIterator.h"

namespace SimulationBox {

/******************************** Constructor *********************************/
template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA>::Iterator
( const int p_area, const VecI p_ncells, const int p_guardsize,
  CellMatrix & p_srcmat )
: area(p_area), core(p_area & CORE), border(p_area & BORDER), guard (p_area &
  GUARD), srcmat(&p_srcmat), icell(0), ncells(p_ncells), guardsize(p_guardsize)
{
    assert( area != 0 ); // would be trivial iterator
    assert( area <= CORE+BORDER+GUARD );
    /* more difficult to implement, and makes no sense */
    assert( area != CORE+GUARD );
}

template<int T_DIM, typename T_CELLDATA> Iterator<T_DIM,T_CELLDATA>::Iterator( void )
: area(-1), core(0), border(0), guard (0), srcmat(NULL), icell(0), ncells(0), guardsize(0)
{}

/******************************** Destructor **********************************/
template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA>::~Iterator( void )
{}

/****************************** Copy Constructor ******************************/
template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA>::Iterator( const Iterator & src )
 : area(src.area), core(src.core), border(src.border), guard(src.guard),
   srcmat(src.srcmat), icell(src.icell), ncells(src.ncells), guardsize(src.guardsize)
{}

/**********************************************************************
 * If we have 2 slates of 3x4 length and we wanted the index          *
 * [i,j,k] = [2,3,2] , then the matrix lies in the memory like:       *
 *   [oooo|oooo|oooo] [oooo|oooo|oxoo]                                *
 * This means to address that element (above marked x) we need to     *
 * calculate:                                                         *
 *   (i=2)*[ (nj=3) * (nk=4) ] + (j=3)*[ (nk=4) ] + (k=2)*[ 1 ]       *
 * That argument in [] will be named 'prevrange'                      *
 **********************************************************************/
template<int T_DIM, typename T_CELLDATA>
int Iterator<T_DIM,T_CELLDATA>::getCellIndex( void ) const {
    int index     = 0;
    int prevrange = 1;
    for (int i=T_DIM-1; i>=0; i--) {
        index     += icell[i] * prevrange;
        prevrange *= ncells[i];
    }
    return index;
}

template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA> Iterator<T_DIM,T_CELLDATA>::begin( void ) const {
    Iterator tmp( *this );

    int upperleftcorner(0);
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
template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA> Iterator<T_DIM,T_CELLDATA>::end( void ) const {
    Iterator tmp( *this );
    tmp.icell  = VecI(-1);
    return tmp;
}

/* Set the iterator to some index */
template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA>& Iterator<T_DIM,T_CELLDATA>::operator<<( const int index[T_DIM] ) {
    VecI ind = VecI(index);
    assert( InArea( ind, area, (ncells-VecI(2*guardsize)), guardsize ) );
    *this = this->begin();
    this->icell += ind;
    return *this;
}

template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA> & Iterator<T_DIM,T_CELLDATA>::operator++( void ) { // only prefix, because it's faster!
    do {
        int linindex = ConvertVectorToLinearIndex( icell, ncells ) + 1;
        if ( !( linindex < ncells.product() ) ) {
            this->icell = VecI(-1);
            break;
        }
        this->icell = ConvertLinearToVectorIndex( linindex, ncells );
    } while( !InArea( icell, area, (ncells-VecI(2*guardsize)), guardsize ) );
    return *this;
}

template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA> Iterator<T_DIM,T_CELLDATA>::operator+( VecI offset ) {
    Iterator shifted = *this;
    shifted.icell += offset;
    VecI localcells = ncells - VecI(2*guardsize);
    if ( InArea( shifted.icell, CORE, localcells, guardsize ) ) {
        shifted.area   = CORE;
        shifted.core   = true;
        shifted.border = false;
        shifted.guard  = false;
    } else if ( InArea( shifted.icell, BORDER, localcells, guardsize ) ) {
        shifted.area   = BORDER;
        shifted.core   = false;
        shifted.border = true;
        shifted.guard  = false;
    } else if ( InArea( shifted.icell, GUARD, localcells, guardsize ) ) {
        shifted.area   = GUARD;
        shifted.core   = false;
        shifted.border = false;
        shifted.guard  = true;
    }
    assert( InArea( shifted.icell, shifted.area, localcells, guardsize ) );
    return shifted;
}

template<int T_DIM, typename T_CELLDATA>
bool Iterator<T_DIM,T_CELLDATA>::operator==( const Iterator & it ) {
    bool alland = true;
    alland = alland & ( this->core   == it.core  );
    alland = alland & ( this->border == it.border );
    alland = alland & ( this->guard  == it.guard  );
    alland = alland & ( this->area   == it.area  );

    alland = alland & ( this->icell  == it.icell  );
    alland = alland & ( this->ncells == it.ncells );
    alland = alland & ( this->srcmat == it.srcmat );

    return alland;
}

template<int T_DIM, typename T_CELLDATA>
bool Iterator<T_DIM,T_CELLDATA>::operator!=( const Iterator & it ) {
    return !( (*this) == it );
}

template<int T_DIM, typename T_CELLDATA>
Iterator<T_DIM,T_CELLDATA> & Iterator<T_DIM,T_CELLDATA>::operator=
( const Iterator & src ) {
    this->core      = src.core  ;
    this->border    = src.border;
    this->guard     = src.guard ;
    this->area      = src.area ;
    this->srcmat    = src.srcmat;
    this->icell     = src.icell ;
    this->ncells    = src.ncells;
    this->guardsize = src.guardsize;
    return *this;
}

template<int T_DIM, typename T_CELLDATA>
T_CELLDATA& Iterator<T_DIM,T_CELLDATA>::operator*( void ) const {
    return (*srcmat)[icell];
}

template<int T_DIM, typename T_CELLDATA>
T_CELLDATA* Iterator<T_DIM,T_CELLDATA>::operator->( void ) const {
    return &((*srcmat)[icell]);
}

}