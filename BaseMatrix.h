#pragma once

#include <cassert>
#include <iostream>
#include "Vector.h"

using namespace std;

#pragma once

#define DEBUG_BASEMATRIX 1


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
        this->size = 0;
        data = NULL;
    }
    BaseMatrix( const VecI size ) {
        this->size = size;
        data = (T_DTYPE*) malloc( sizeof(T_DTYPE) * size.product() );
        #if DEBUG_BASEMATRIX >= 2
            cout << "Constructor of size: "; size.Print(); cout << endl;
        #endif
    }

    ~BaseMatrix( void ) {
        if (data != NULL)
            free( data );
    }

//private:
    /**************************************************************************
     * If we have 2 slates of 3x4 length and we wanted the index  [i,j,k] =   *
     * [1,2,1] (begin with 0!), then the matrix lies in the memory like:      *
     *   [oooo|oooo|oooo] [oooo|oooo|oxoo]                                    *
     * This means to address that element (above marked x) we need to         *
     * calculate:                                                             *
     *   21 = (i=1)*[ (nj=3) * (nk=4) ] + (j=2)*[ (nk=4) ] + (k=1)*[ 1 ]      *
     * That argument in [] will be named 'prevrange'                          *
     **************************************************************************/
    int getLinearIndex( VecI pos ) {
        #if DEBUG_BASEMATRIX >= 2
            cout << "pos: "; pos.Print();
            cout << " size: "; size.Print();
            cout << endl;
        #endif
        assert( pos < size );
        int index     = 0;
        int prevrange = 1;
        for (int i=T_DIMENSION-1; i>=0; i--) {
            index     += pos[i] * prevrange;
            prevrange *= size[i];
        }
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
    VecI getVectorIndex( const int linindex ) {
        assert( linindex < size.product() );
        VecI index;
        int tmp = linindex;
        for (int i=T_DIMENSION-1; i>=0; i--) {
            index[i]  = tmp % size[i];
            tmp       = tmp / size[i];
        }
        assert( tmp == 0 );
        return index;
    }

public:
    /**************************** Access Operators ****************************/
    T_DTYPE & operator[] ( const VecI pos ) const {
        return data[ getLinearIndex(pos) ];
    }

    T_DTYPE & operator[] ( const VecI pos ) {
        return data[ getLinearIndex(pos) ];
    }

    T_DTYPE & operator[] ( const int i ) {
        assert( i < size.product() );
        return data[i];
    }

    /************************** Assignment Operators **************************/
    /*BaseMatrix& operator= (const T_DTYPE a) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] = a;
        return *this;
    }*/

    BaseMatrix& operator= (const BaseMatrix & v) {
        this->~BaseMatrix();
        this->size = v.size;
        data = (T_DTYPE*) malloc( sizeof(T_DTYPE) * size.product() );
        #if DEBUG_BASEMATRIX >= 2
            cout << "Assignment of size: "; v.size.Print(); cout << endl;
        #endif

        memcpy( this->data, v.data, sizeof(T_DTYPE) * size.product() );
        return *this;
    }

    BaseMatrix getPartialMatrix( VecI pos, VecI size ) {
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

#if 1 == 0
        int currentaxis = T_DIMENSION-1;
        while(true) {
            /* Actual operation we want to do! Rest is only another way to    *
             * write arbitrary for-loops over all indices ... -> reverse      *
             * getLinearIndex would have been easier -.-                      */
            tmp[ind] = this->data[pos+ind];
            ind[currentaxis]++;
            while( (ind[currentaxis] >= size[currentaxis]) and (currentaxis >= 0) ) {
                ind[currentaxis] = 0;
                curentaxis--;
                ind[currentaxis]++;
            }
            /* Overflow => all elements cycled through */
            if (currentaxis < 0)
                break;
        }
#endif

        return tmp;
    }
    
    void insertMatrix( VecI pos, BaseMatrix m ) {
        assert( pos+m.size <= this->size );
        for ( int i=0; i<m.size.product(); i++ )
            (*this)[ pos + m.getVectorIndex(i) ] = m[i];
    }
};

