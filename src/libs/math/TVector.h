#include <cassert>
#include <iostream>
#include <cmath>     // std::abs

#pragma once

#define DEBUG_VECTOR 1


template<typename T_DTYPE, int T_DIM>
class Vec {
public:
    const int dim = T_DIM;
    const int & size = dim;
    T_DTYPE data[T_DIM];

    inline Vec(void);
    inline Vec(const Vec & v);
    inline ~Vec(void);
    inline Vec& operator=(const Vec & v);

    template<typename T_ETYPE> inline Vec(const T_ETYPE v[T_DIM]);
    template<typename T_ETYPE> inline Vec(const T_ETYPE a0, const T_ETYPE a1);
    template<typename T_ETYPE> inline Vec(const T_ETYPE a0, const T_ETYPE a1, const T_ETYPE a2);
    template<typename T_ETYPE> inline Vec(const Vec<T_ETYPE,T_DIM> & v);
    inline Vec(const T_DTYPE a);

    inline T_DTYPE operator[] (const int i) const;
    inline T_DTYPE & operator[] (const int i);

    inline operator const T_DTYPE*() const;

    inline Vec& operator= (const T_DTYPE a);
    inline Vec& operator+=(const Vec & v);
    inline Vec& operator*=(const Vec & v);
    inline Vec& operator/=(const Vec & v);
    inline Vec  operator/ (const Vec & v) const;

    inline bool operator==(const Vec & v) const;
    inline bool operator< (const Vec & v) const;
    inline bool operator> (const Vec & v) const;
    inline bool operator<=(const Vec & v) const;
    inline bool operator>=(const Vec & v) const;
    inline Vec<bool,T_DIM> GreaterThan(const Vec & v) const;
    inline Vec<bool,T_DIM> GreaterOrEqualThan(const Vec & v) const;
    inline Vec<bool,T_DIM> SmallerThan(const Vec & v) const;
    inline Vec min(const Vec & v) const;
    inline Vec max(const Vec & v) const;

    inline T_DTYPE scp(const Vec & v) const;
    inline Vec abs(void) const;
    inline double norm2( void ) const;
    inline T_DTYPE product( void ) const;
    inline T_DTYPE sum( void ) const;
    inline T_DTYPE min( void ) const;
    inline T_DTYPE max( void ) const;
    inline T_DTYPE mean( void ) const;
    inline Vec cross(const Vec & b) const;

    inline Vec  operator+ (const Vec & v) const;
    inline Vec& operator-=(const Vec & v);
    inline Vec  operator- (const Vec & v) const;
    inline Vec  operator* (const Vec & v) const;

    inline Vec& operator+=(const T_DTYPE a);
    inline Vec& operator-=(const T_DTYPE a);
    inline Vec& operator*=(const T_DTYPE a);
    inline Vec& operator/=(const T_DTYPE a);
    inline Vec  operator+ (const T_DTYPE a) const;
    inline Vec  operator- (const T_DTYPE a) const;
    inline Vec  operator* (const T_DTYPE a) const;
    inline Vec  operator/ (const T_DTYPE a) const;

    int getSize(void) const;

    inline bool operator!= (const Vec & v) const;
    inline double norm() const;
};


template<typename T, typename T_DTYPE, int T_DIM>
inline Vec<T_DTYPE,T_DIM> operator+
( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );

template<typename T, typename T_DTYPE, int T_DIM>
inline Vec<T_DTYPE,T_DIM> operator-
( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );

template<typename T, typename T_DTYPE, int T_DIM>
inline Vec<T_DTYPE,T_DIM> operator*
( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );

template<typename T, typename T_DTYPE, int T_DIM>
inline Vec<T_DTYPE,T_DIM> operator/
( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside );

template<typename T_DTYPE, int T_DIM>
inline std::ostream& operator<<( std::ostream& out, const Vec<T_DTYPE,T_DIM> v );

#include "TVector.tpp"
