#pragma once

#include <cassert>
#include <iostream>
#include <cstring> // memcpy
#include "Vector.h"

using namespace std;

#pragma once

#define DEBUG_BASEMATRIX 2


template<typename T_DTYPE, int T_DIMENSION>
class BaseMatrix {
public:
    static const int dim = T_DIMENSION;
    typedef Vec<int,T_DIMENSION> VecI;
    typedef T_DTYPE Datatype;

    T_DTYPE* data;
    VecI size;

    /****************************** Constructors ******************************/
    BaseMatrix( void ) {
        #if DEBUG_BASEMATRIX >= 2
            cout << "Standard Constructor called" << endl << flush;
        #endif
        this->size = 0;
        this->data = NULL;
    }
    
    BaseMatrix( VecI size ) { 
        #if DEBUG_BASEMATRIX >= 2
            cout << "Constructor of size: "; size.Print(); cout << endl;
        #endif
        this->size = size;
        this->data = (T_DTYPE*) malloc( sizeof(T_DTYPE) * size.product() );
    }
    
    /* Copy Constructor. It's a constructor, so data and size are assumed to  *
     * be absolutely random!                                                  */
    BaseMatrix( const BaseMatrix & m ) { 
        #if DEBUG_BASEMATRIX >= 2
            cout << "Copy Constructor called with Matrix of size: "; m.size.Print(); cout << endl;
        #endif
        this->size = m.size;
        this->data = (T_DTYPE*) malloc( sizeof(T_DTYPE) * size.product() );
        memcpy( this->data, m.data, sizeof(T_DTYPE) * m.size.product() );
    }

    ~BaseMatrix( void ) {
        if ( this->data != NULL )
            free( this->data );
    }

    /************************** Assignment Operators **************************/
    BaseMatrix& operator= (const T_DTYPE a) {
        for (int i=0; i < size.product(); i++)
            this->data[i] = a;
        return *this;
    }

    BaseMatrix& operator= (const BaseMatrix & m) {
        #if DEBUG_BASEMATRIX >= 2
            cout << "Assignment of size: "; m.size.Print(); cout << endl << flush;
        #endif
        this->~BaseMatrix();
        this->size = m.size;
        this->data = (T_DTYPE*) malloc( sizeof(T_DTYPE) * size.product() );

        memcpy( this->data, m.data, sizeof(T_DTYPE) * m.size.product() );
        return *this;
    }

    /**************************************************************************
     * If we have 2 slates of 3x4 length and we wanted the index  [i,j,k] =   *
     * [1,2,1] (begin with 0!), then the matrix lies in the memory like:      *
     *   [oooo|oooo|oooo] [oooo|oooo|oxoo]                                    *
     * This means to address that element (above marked x) we need to         *
     * calculate:                                                             *
     *   21 = (i=1)*[ (nj=3) * (nk=4) ] + (j=2)*[ (nk=4) ] + (k=1)*[ 1 ]      *
     * That argument in [] will be named 'prevrange'                          *
     **************************************************************************/
    int getLinearIndex( const VecI & pos ) const {
        assert( pos < size );
        assert( pos >= VecI(0) );
        int index     = 0;
        int prevrange = 1;
        for (int i=T_DIMENSION-1; i>=0; i--) {
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
    VecI getVectorIndex( const int & linindex ) const {
        assert( linindex < this->size.product() );
        VecI index;
        int tmp = linindex;
        for (int i=T_DIMENSION-1; i>=0; i--) {
            index[i]  = tmp % this->size[i];
            tmp       = tmp / this->size[i];
        }
        assert( tmp == 0 );
        assert( index < this->size );
        assert( index >= VecI(0) );
        return index;
    }

    /**************************** Access Operators ****************************/
    T_DTYPE operator[] ( const int i ) const {
        assert( i < size.product() );
        return data[i];
    }
    
    T_DTYPE & operator[] ( const int i ) {
        assert( i < size.product() );
        return data[i];
    }
    
    T_DTYPE operator[] ( const VecI pos ) const {        
        #if DEBUG_BASEMATRIX >= 1
        if ( ! (pos < size) ) {
            cout << "Incorrect Access try to: ";
            pos.Print();
            cout << " in ";
            this->size.Print();
            cout << "-Matrix";
            cout << endl;
        }
        #endif
        assert( pos < size );
        return data[ getLinearIndex(pos) ];
    }

    T_DTYPE & operator[] ( const VecI pos ) {
        assert( pos < size );
        return data[ getLinearIndex(pos) ];
    }


    /***************************** Other Functions ****************************/
    VecI getSize( void ) const {
        return this->size;
    }
    
    BaseMatrix getPartialMatrix( const VecI & pos, const VecI & size ) const {
        #if DEBUG_BASEMATRIX >= 2
            cout << "pos: "; pos.Print();
            cout << " size: "; size.Print();
            cout << " this->size: "; this->size.Print();
            cout << endl;
        #endif
        assert( pos+size <= this->size );
        BaseMatrix<T_DTYPE,T_DIMENSION> tmp( size );
        VecI ind(0);
        
        for ( int i=0; i<size.product(); i++ )
            tmp[i] = (*this)[ pos + tmp.getVectorIndex(i) ];

        return tmp;
    }
    
    void insertMatrix( const VecI & pos, const BaseMatrix & m ) {
#if DEBUG_BASEMATRIX >= 2
        cout << "(pos="; pos.Print(); cout << ") + (m.size=";
        m.size.Print(); cout << ") <= (this->size="; this->size.Print();
        cout << "?" << endl << flush;
#endif
        assert( pos+m.getSize() <= this->size );
        for ( int i=0; i < m.getSize().product(); i++ ) {
            VecI index = pos + m.getVectorIndex(i);
            assert( index < this->size );
            (*this)[ index ] = m[i];
        }
    }
    
#if DEBUG_BASEMATRIX >= 1
    void Print( void ) const {
        VecI size = this->getSize();
        cout << "This "; size.Print();
        cout << "Matrix:" << endl;
        VecI ind(0);
        for ( ind[1]=0; ind[1]<size[1]; ind[1]++) {
            for ( ind[0]=0; ind[0]<size[0]; ind[0]++)
                cout << (*this)[ ind ] << " ";
            cout << endl;
        }
        cout << endl;
    }
#endif
};

