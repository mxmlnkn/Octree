#pragma once

#include <cassert>
#include <iostream>
#include <cstring> // memcpy
#include "TVector.h"

using namespace std;

#pragma once

#define DEBUG_BASEMATRIX 2

//namespace BaseMatrix

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
int ConvertVectorToLinearIndex( const Vec<int,T_DIM> & pos, const Vec<int,T_DIM> & size ) {
    assert( pos < size );
    assert( pos >= (Vec<int,T_DIM>)(0) );
    int index     = 0;
    int prevrange = 1;
    for (int i=T_DIM-1; i>=0; i--) {
        index     += pos[i] * prevrange;
        prevrange *= size[i];
    }
    assert( index < size.product() );
    return index;
}

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
Vec<int,T_DIM> ConvertLinearToVectorIndex( const int & linindex, const Vec<int,T_DIM> & size ) {
    assert( linindex < size.product() );
    Vec<int,T_DIM> index;
    int tmp = linindex;
    for (int i=T_DIM-1; i>=0; i--) {
        index[i]  = tmp % size[i];
        tmp       = tmp / size[i];
    }
    assert( tmp == 0 );
    assert( index <  size );
    assert( index >= (Vec<int,T_DIM>)(0) );
    return index;
}
    

template<typename T_DTYPE, int T_DIM>
class BaseMatrix {
private:
    Vec<int,T_DIM> size;

public:
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



/******************************* Magic Methods ********************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( void ) {
    this->size = 0;
    this->data = NULL;
}

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( const BaseMatrix & m ) {
    this->size = m.size;
    this->data = new T_DTYPE[size.product()];
    memcpy( this->data, m.data, sizeof(T_DTYPE) * m.size.product() );
}

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::~BaseMatrix( void ) {
    if ( this->data != NULL )
        delete[] this->data;
}

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>& BaseMatrix<T_DTYPE,T_DIM>::operator= (const BaseMatrix & m) {
    this->~BaseMatrix();
    this->size = m.size;
    this->data = new T_DTYPE[size.product()];

    memcpy( this->data, m.data, sizeof(T_DTYPE) * m.size.product() );
    return *this;
}

/***************************** Other Constructors *****************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( VecI size ) {
    this->size = size;
    this->data = new T_DTYPE[size.product()];
}

/************************** Assignment Operators **************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>& BaseMatrix<T_DTYPE,T_DIM>::operator= (const T_DTYPE a) {
    for (int i=0; i < size.product(); i++)
        this->data[i] = a;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
int BaseMatrix<T_DTYPE,T_DIM>::getLinearIndex( const VecI & pos ) const {
    return ConvertVectorToLinearIndex( pos, this->size );
}

template<typename T_DTYPE, int T_DIM>
Vec<int,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getVectorIndex( const int & linindex ) const {
    return ConvertLinearToVectorIndex( linindex, this->size );
}

/**************************** Access Operators ****************************/
template<typename T_DTYPE, int T_DIM>
T_DTYPE BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const int i ) const {
    assert( i < size.product() );
    return data[i];
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE & BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const int i ) {
    assert( i < size.product() );
    return data[i];
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const VecI pos ) const {
    assert( pos < size );
    return data[ getLinearIndex(pos) ];
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE & BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const VecI pos ) {
    assert( pos < size );
    return data[ getLinearIndex(pos) ];
}


/***************************** Other Functions ****************************/
template<typename T_DTYPE, int T_DIM>
Vec<int,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getSize( void ) const {
    return this->size;
}

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getPartialMatrix( const VecI & pos, const VecI & size ) const {
    assert( pos+size <= this->size );
    BaseMatrix<T_DTYPE,T_DIM> tmp( size );
    VecI ind(0);

    for ( int i=0; i<size.product(); i++ )
        tmp[i] = (*this)[ pos + tmp.getVectorIndex(i) ];

    return tmp;
}

template<typename T_DTYPE, int T_DIM>
void BaseMatrix<T_DTYPE,T_DIM>::insertMatrix( const VecI & pos, const BaseMatrix & m ) {
    assert( pos+m.getSize() <= this->size );
    for ( int i=0; i < m.getSize().product(); i++ ) {
        VecI index = pos + m.getVectorIndex(i);
        assert( index < this->size );
        (*this)[ index ] = m[i];
    }
}


/* Enables cout << Vec<int,2>(1); This alos works with fstream and therefore with cout */
template<typename T_DTYPE, int T_DIM>
ostream& operator<<( ostream& out, const BaseMatrix<T_DTYPE,T_DIM>& m ) {
    typedef typename BaseMatrix<T_DTYPE,T_DIM>::VecI VecI; // xD ... typisch picongpu
    VecI size = m.getSize();
    out << "This " << size << "Matrix:" << endl;
    VecI ind(0);
    for ( ind[1]=0; ind[1]<size[1]; ind[1]++) {
        for ( ind[0]=0; ind[0]<size[0]; ind[0]++)
            out << m[ ind ] << " ";
        out << endl;
    }
    out << endl;
    return out;
}
