#pragma once

#include <cassert>
#include <iostream>
#include <cstring>   // memcpy
#include "TVector.h"

#define DEBUG_BASEMATRIX 2

/**************************************************************************
 * If we have 2 slates of 3x4 length and we wanted the index  [i,j,k] =   *
 * [1,2,1] (begin with 0!), then the matrix lies in the memory like:      *
 *   [oooo|oooo|oooo] [oooo|oooo|oxoo]                                    *
 * This means to address that element (above marked x) we need to         *
 * calculate:                                                             *
 *   21 = (i=1)*[ (nj=3) * (nk=4) ] + (j=2)*[ (nk=4) ] + (k=1)*[ 1 ]      *
 * That argument in [] will be named 'prevrange'                          *
 **************************************************************************/
template<int T_DIM>
int ConvertVectorToLinearIndex( const Vec<int,T_DIM> & pos, const Vec<int,T_DIM> & size );

/**************************************************************************
 * If we write the above formula in a different way we can derive the     *
 * reverse operation, meaning getting the indices with modulo and div:    *
 *   21 = (k=1) + (nk=4)*[ (j=2) + (nj=3)*[(i=1)] ]                       *
 * So mod (nk=4) will give us (k=2), then div (nk=4) will give us the     *
 * number we have to consider recursively:                                *
 *   k   = 21 mod (nk=4) = 1                                              *
 *   tmp = 21  /  (nk=4) = 5                                              *
 *   j   = 5  mod (nj=3) = 2                                              *
 *   tmp = 5   /  (nj=3) = 1                                              *
 *   i   = 1  mod (ni=2) = 1                                              *
 * This works, because every new summand has the previous factor in       *
 * prevrange, so that summand mod factor = 0                              *
 **************************************************************************/
template<int T_DIM>
Vec<int,T_DIM> ConvertLinearToVectorIndex( const int & linindex, const Vec<int,T_DIM> & size );

template<typename T_DTYPE, int T_DIM>
class BaseMatrix {
public:
    Vec<int,T_DIM> size;
    T_DTYPE* data;

    static const int dim = T_DIM;
    typedef Vec<int,T_DIM> VecI;
    typedef T_DTYPE Datatype;

	/* should only be used in conjunction with following copy constructor,    *
	 * because this initializes a matrix of size 0 with no way to adjust the  *
	 * the size thereafter                                                    */
    BaseMatrix( void );
    BaseMatrix( const BaseMatrix & m );
    ~BaseMatrix( void );
    BaseMatrix& operator= (const BaseMatrix & m);

    BaseMatrix( VecI size );
    BaseMatrix& operator= (const T_DTYPE a);
    int getLinearIndex( const VecI & pos ) const;
    VecI getVectorIndex( const int & linindex ) const;
    T_DTYPE operator[] ( const int i ) const;
    T_DTYPE & operator[] ( const int i );
    T_DTYPE operator[] ( const VecI pos ) const;
    T_DTYPE & operator[] ( const VecI pos );
    VecI getSize( void ) const;
    BaseMatrix getPartialMatrix( const VecI & pos, const VecI & size ) const;
    void insertMatrix( const VecI & pos, const BaseMatrix & m );
};

#include "TBaseMatrix.tpp"
