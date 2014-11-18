#include <cassert>
#include <iostream>

using namespace std;

#pragma once

#define DEBUG_VECTOR 1


template<typename T_DTYPE, int T_DIMENSION>
class Vec {
public:
    static const int dim = T_DIMENSION;
    T_DTYPE data[T_DIMENSION];
    /**************************************************************************
        Vec(void);
        Vec(const T_ETYPE v[T_DIMENSION] );
        Vec(const Vec<T_ETYPE,T_DIMENSION> & v);
        Vec(const T_DTYPE a);

        T_DTYPE operator[] (const int i) const;
        T_DTYPE & operator[] (const int i);

        operator const T_DTYPE*() const;

        Vec& operator= (const T_DTYPE a);
        Vec& operator= (const Vec<T_ETYPE, T_DIMENSION> & v);
        Vec& operator+= (const Vec<T_ETYPE, T_DIMENSION> & v);
        Vec& operator*= (const Vec<T_ETYPE,T_DIMENSION> & v);
        inline Vec& operator/= (const Vec<T_ETYPE, T_DIMENSION> & v);
        Vec operator/ (const Vec<T_ETYPE, T_DIMENSION> & v) const;

        bool operator== (const Vec & v) const;
        bool operator< (const Vec & v) const;
        bool operator> (const Vec & v) const;
        bool operator<= (const Vec & v) const;
        bool operator>= (const Vec & v) const;

        double norm2() const;
        T_DTYPE product( void ) const;
        T_DTYPE sum( void ) const;

        inline Vec operator+ (const Vec<T_ETYPE, T_DIMENSION> & v) const;
        inline Vec& operator-= (const Vec<T_ETYPE, T_DIMENSION> & v);
        inline Vec operator- (const Vec<T_ETYPE, T_DIMENSION> & v) const;
        inline Vec operator* (const T_ETYPE a) const;
        inline Vec operator* (const Vec<T_ETYPE,T_DIMENSION> & v) const;
        Vec& operator*= (const T_ETYPE a);
        inline Vec& operator/= (const T_ETYPE a);
        inline Vec operator/ (const T_ETYPE a) const;
        inline bool operator!= (const Vec & v) const;
        double norm() const;
        void Print( void ) const;
     **************************************************************************/

    /****************************** Constructors ******************************/
    Vec(void) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] = 0;
    }

    template<typename T_ETYPE>
    Vec(const T_ETYPE v[T_DIMENSION] ) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] = v[i];
    }

    template<typename T_ETYPE>
    Vec(const Vec<T_ETYPE,T_DIMENSION> & v) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] = v[i];
    }

    Vec(const T_DTYPE a) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] = a;
    }

    /**************************** Access Operators ****************************/
    T_DTYPE operator[] (const int i) const {
        assert( i>= 0 and i < T_DIMENSION );
        return data[i];
    }

    T_DTYPE & operator[] (const int i) {
        assert( i>= 0 and i < T_DIMENSION );
        return data[i];
    }

    /************************** Conversion Operators **************************/
    operator const T_DTYPE*() const {
        return this->data;
    }

    /************************** Assignment Operators **************************/
    Vec& operator= (const T_DTYPE a) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] = a;
        return *this;
    }

    template<typename T_ETYPE>
    Vec& operator= (const Vec<T_ETYPE, T_DIMENSION> & v) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] = v[i];
        return *this;
    }

    template<typename T_ETYPE>
    Vec& operator+= (const Vec<T_ETYPE, T_DIMENSION> & v) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] += v[i];
        return *this;
    }

    /* Elementwise Multiplication with another Vec of double or int */
    template<typename T_ETYPE>
    Vec& operator*= (const Vec<T_ETYPE,T_DIMENSION> & v) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] *= v[i];
        return *this;
    }

    /* Elementwise Division with another Vec of double or int */
    template<typename T_ETYPE>
    inline Vec& operator/= (const Vec<T_ETYPE, T_DIMENSION> & v) {
        for (int i=0; i<T_DIMENSION; i++)
            this->data[i] /= v[i];
        return *this;
    }

    /* Elementwise Division with another Vec of double or int */
    template<typename T_ETYPE>
    Vec operator/ (const Vec<T_ETYPE, T_DIMENSION> & v) const {
        Vec tmp(*this);
        tmp /= v;
        return tmp;
    }

    /************************** Comparison Operators **************************/
    bool operator== (const Vec & v) const {
        bool alland = true;
        for (int i=0; i<T_DIMENSION; i++)
            alland = alland & ( this->data[i] == v[i] );
        return alland;
    }

    bool operator< (const Vec & v) const {
        bool alland = true;
        for (int i=0; i<T_DIMENSION; i++)
            alland = alland & ( this->data[i] < v[i] );
        return alland;
    }

    bool operator> (const Vec & v) const {
        bool alland = true;
        for (int i=0; i<T_DIMENSION; i++)
            alland = alland & ( this->data[i] > v[i] );
        return alland;
    }

    bool operator<= (const Vec & v) const {
        bool alland = true;
        for (int i=0; i<T_DIMENSION; i++)
            alland = alland & ( this->data[i] <= v[i] );
        return alland;
    }

    bool operator>= (const Vec & v) const {
        bool alland = true;
        for (int i=0; i<T_DIMENSION; i++)
            alland = alland & ( this->data[i] >= v[i] );
        return alland;
    }


    /********************** Horizontal Vector Operations **********************/

    /* Scalarproduct */
    /*T_DTYPE operator* (const Vec & v) const {
        T_DTYPE sum = 0;
        for (int i=0; i<T_DIMENSION; i++)
            sum += (this->data[i]) * v[i];
        return sum;
    }*/

    double norm2() const {
        double sum = 0;
        for (int i=0; i<T_DIMENSION; i++)
            sum += (this->data[i]) * (this->data[i]);
        return sum;
    }

    T_DTYPE product( void ) const {
        T_DTYPE tmp = 1;
        for (int i=0; i<T_DIMENSION; i++)
            tmp *= this->data[i];
        return tmp;
    }

    T_DTYPE sum( void ) const {
        T_DTYPE tmp = 0;
        for (int i=0; i<T_DIMENSION; i++)
            tmp += this->data[i];
        return tmp;
    }

    /*Vec cross(const Vec & b) const {
        Vec res(0,0,0);
        res.x += this->y * b.z - this->z * b.y;
        res.y += this->z * b.x - this->x * b.z;
        res.z += this->x * b.y - this->y * b.x;
        return res;
    } */

    /************************ Derived Declarations ********************************/
    template<typename T_ETYPE>
    inline Vec operator+ (const Vec<T_ETYPE, T_DIMENSION> & v) const {
        Vec res = *this;
        res += v;
        return res;
    }

    template<typename T_ETYPE>
    inline Vec& operator-= (const Vec<T_ETYPE, T_DIMENSION> & v) {
        return (*this) += v * (-1);
    }

    template<typename T_ETYPE>
    inline Vec operator- (const Vec<T_ETYPE, T_DIMENSION> & v) const {
        return (*this)+( v*(-1) );
    }

    /* Elementwise Multiplication, derived from *= */
    template<typename T_ETYPE>
    inline Vec operator* (const T_ETYPE a) const {
        Vec res = *this;
        res *= a;
        return res;
    }
    /* Broadcasting Multiplication, derived from *= */
    template<typename T_ETYPE>
    inline Vec operator* (const Vec<T_ETYPE,T_DIMENSION> & v) const {
        Vec res = *this;
        res *= v;
        return res;
    }

    /* Broadcasting Multiplication, derived from elementwise Mul */
    template<typename T_ETYPE>
    Vec& operator*= (const T_ETYPE a) {
        (*this) *= Vec(a);
        return *this;
    }
    /* Broadcasting division from double or int, derived from elementwise div */
    template<typename T_ETYPE>
    inline Vec& operator/= (const T_ETYPE a) {
        (*this) /= a;
        return *this;
    }
    /* Broadcasting division from double or int, derived from elementwise mul */
    template<typename T_ETYPE>
    inline Vec operator/ (const T_ETYPE a) const {
        return (*this)*(Vec(1.0)/a);
    }


    inline bool operator!= (const Vec & v) const {
        return !((*this) == v);
    }

    double norm() const {
        double res = sqrt((*this).norm2());
        return res;
    }

    /*************************** Debug methods ********************************/
#if DEBUG_VECTOR >= 1
    void Print( void ) const {
        cout << "(";
        for (int i=0; i<T_DIMENSION-1; i++)
            cout << this->data[i] << ",";
        cout << this->data[T_DIMENSION-1] << ")";
    }
#endif
};

template<typename T, typename T_DTYPE, int T_DIMENSION, typename T_ETYPE>
Vec<T_DTYPE,T_DIMENSION> operator*(T const& scalar, Vec<T_DTYPE,T_DIMENSION> rhs)
{
    // scalar multiplication is commutative: s M = M s
    return rhs *= scalar; // calls rhs.operator*=(scalar);
}