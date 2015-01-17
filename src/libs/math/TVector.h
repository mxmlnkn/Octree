#include <cassert>
#include <iostream>

#pragma once

#define DEBUG_VECTOR 1


template<typename T_DTYPE, int T_DIM>
class Vec {
public:
    const int dim = T_DIM;
    const int & size = dim;
    T_DTYPE data[T_DIM];
    
    Vec(void);
    Vec(const Vec & v);
    ~Vec(void);
    Vec& operator=(const Vec & v);

    template<typename T_ETYPE> Vec(const T_ETYPE v[T_DIM] );
    template<typename T_ETYPE> Vec(const Vec<T_ETYPE,T_DIM> & v);
    Vec(const T_DTYPE a);

    T_DTYPE operator[] (const int i) const;
    T_DTYPE & operator[] (const int i);

    operator const T_DTYPE*() const;

    Vec& operator= (const T_DTYPE a);
    template<typename T_ETYPE> Vec& operator= (const Vec<T_ETYPE, T_DIM> & v);
    template<typename T_ETYPE> Vec& operator+=(const Vec<T_ETYPE, T_DIM> & v);
    template<typename T_ETYPE> Vec& operator*=(const Vec<T_ETYPE, T_DIM> & v);
    template<typename T_ETYPE> Vec& operator/=(const Vec<T_ETYPE, T_DIM> & v);
    template<typename T_ETYPE> Vec  operator/ (const Vec<T_ETYPE, T_DIM> & v) const;

    bool operator==(const Vec & v) const;
    bool operator< (const Vec & v) const;
    bool operator> (const Vec & v) const;
    bool operator<=(const Vec & v) const;
    bool operator>=(const Vec & v) const;
    Vec<bool,T_DIM> GreaterThan(const Vec & v) const;
    Vec<bool,T_DIM> SmallerThan(const Vec & v) const;

    double norm2( void ) const;
    T_DTYPE product( void ) const;
    T_DTYPE sum( void ) const;
    T_DTYPE min( void ) const;
    T_DTYPE max( void ) const;
    T_DTYPE mean( void ) const;

    template<typename T_ETYPE> Vec  operator+ (const Vec<T_ETYPE, T_DIM> & v) const;
    template<typename T_ETYPE> Vec& operator-=(const Vec<T_ETYPE, T_DIM> & v);
    template<typename T_ETYPE> Vec  operator- (const Vec<T_ETYPE, T_DIM> & v) const;
    template<typename T_ETYPE> Vec  operator* (const Vec<T_ETYPE, T_DIM> & v) const;
    
    template<typename T_ETYPE> Vec& operator+=(const T_ETYPE a);
    template<typename T_ETYPE> Vec& operator-=(const T_ETYPE a);
    template<typename T_ETYPE> Vec& operator*=(const T_ETYPE a);
    template<typename T_ETYPE> Vec& operator/=(const T_ETYPE a);
    template<typename T_ETYPE> Vec  operator+ (const T_ETYPE a) const;
    template<typename T_ETYPE> Vec  operator- (const T_ETYPE a) const;
    template<typename T_ETYPE> Vec  operator* (const T_ETYPE a) const;
    template<typename T_ETYPE> Vec  operator/ (const T_ETYPE a) const;

    bool operator!= (const Vec & v) const;
    double norm() const;
};


template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator+( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator-( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator*( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator/( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );

template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out, const Vec<T_DTYPE,T_DIM> v );

#include "TVector.tpp"