#pragma once

#include "TVector.h"


/******************************* Magic Methods ********************************/

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>::Vec( void ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = 0;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>::Vec( const Vec & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = v.data[i];
        // have to use .data instead of v::operator[], which would result in "error: invalid use of template-name ‘Vec’ without an argument list"
        // Also this Argument is the only exception where we don't need to specify the template arguments for Vec !!! "const Vec &" suffices instead of "const Vec<T_DTYPE,T_DIM> &", which would even be wrong !!!
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>::~Vec( void ) {
    return;
}

/********************************* Assignment *********************************/
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator=( const Vec & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = v.data[i];
    return *this;
}

/******************************** Constructors ********************************/

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>::Vec( const T_ETYPE v[T_DIM] ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = (T_DTYPE) v[i];
}

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>::Vec( const T_ETYPE a0, const T_ETYPE a1 ) {
    assert( T_DIM <= 3 );
    this->data[0] = (T_DTYPE) a0;
    if ( T_DIM >= 2 )
        this->data[1] = (T_DTYPE) a1;
    for ( int i=2; i<T_DIM; i++ )
        data[i] = 0;
}

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>::Vec( const T_ETYPE a0, const T_ETYPE a1, const T_ETYPE a2 ) {
    assert( T_DIM <= 3 );
    this->data[0] = (T_DTYPE) a0;
    if ( T_DIM >= 2 )
        this->data[1] = (T_DTYPE) a1;
    if ( T_DIM >= 3 )
        this->data[2] = (T_DTYPE) a2;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>::Vec(const T_DTYPE a) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = a;
}

/****************************** Copy Constructor ******************************/
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>::Vec( const Vec<T_ETYPE,T_DIM> & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = (T_DTYPE) v.data[i];
}

/****************************** Access Operators ******************************/
template<typename T_DTYPE, int T_DIM>
T_DTYPE Vec<T_DTYPE,T_DIM>::operator[] (const int i) const {
    assert( i>= 0 and i < T_DIM );
    return data[i];
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE & Vec<T_DTYPE,T_DIM>::operator[] (const int i) {
    assert( i>= 0 and i < T_DIM );
    return data[i];
}

/************************** Conversion Operators **************************/
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>::operator const T_DTYPE*() const {
    return this->data;
}

/************************** Assignment Operators **************************/
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator=( const T_DTYPE a ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = a;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator+=( const Vec & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] += v[i];
    return *this;
}

/* Elementwise Multiplication with another Vec of double or int */
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator*=( const Vec & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] *= v[i];
    return *this;
}

/* Elementwise Division with another Vec of double or int */
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator/=( const Vec & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] /= v[i];
    return *this;
}

/* Elementwise Division with another Vec of double or int */
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator/( const Vec & v ) const {
    Vec tmp(*this);
    tmp /= v;
    return tmp;
}

/************************** Comparison Operators **************************/
template<typename T_DTYPE, int T_DIM>
bool Vec<T_DTYPE,T_DIM>::operator==( const Vec & v ) const {
    bool alland = true;
    for (int i=0; i<T_DIM; i++)
        alland = alland & ( this->data[i] == v[i] );
    return alland;
}

template<typename T_DTYPE, int T_DIM>
bool Vec<T_DTYPE,T_DIM>::operator<( const Vec & v ) const {
    bool alland = true;
    for (int i=0; i<T_DIM; i++)
        alland = alland & ( this->data[i] < v[i] );
    return alland;
}

template<typename T_DTYPE, int T_DIM>
bool Vec<T_DTYPE,T_DIM>::operator>( const Vec & v ) const {
    bool alland = true;
    for (int i=0; i<T_DIM; i++)
        alland = alland & ( this->data[i] > v[i] );
    return alland;
}

template<typename T_DTYPE, int T_DIM>
bool Vec<T_DTYPE,T_DIM>::operator<=( const Vec & v ) const {
    bool alland = true;
    for (int i=0; i<T_DIM; i++)
        alland = alland & ( this->data[i] <= v[i] );
    return alland;
}

template<typename T_DTYPE, int T_DIM>
bool Vec<T_DTYPE,T_DIM>::operator>=( const Vec & v ) const {
    bool alland = true;
    for (int i=0; i<T_DIM; i++)
        alland = alland & ( this->data[i] >= v[i] );
    return alland;
}

template<typename T_DTYPE, int T_DIM>
Vec<bool,T_DIM> Vec<T_DTYPE,T_DIM>::GreaterThan( const Vec & v ) const {
    Vec<bool,T_DIM> res;
    for (int i=0; i<T_DIM; i++)
        res[i] = this->data[i] > v[i];
    return res;
}

template<typename T_DTYPE, int T_DIM>
Vec<bool,T_DIM> Vec<T_DTYPE,T_DIM>::GreaterOrEqualThan( const Vec & v ) const {
    Vec<bool,T_DIM> res;
    for (int i=0; i<T_DIM; i++)
        res[i] = this->data[i] >= v[i];
    return res;
}

template<typename T_DTYPE, int T_DIM>
Vec<bool,T_DIM> Vec<T_DTYPE,T_DIM>::SmallerThan( const Vec & v ) const {
    Vec<bool,T_DIM> res;
    for (int i=0; i<T_DIM; i++)
        res[i] = this->data[i] < v[i];
    return res;
}


/********************** Horizontal Vector Operations **********************/

/* Scalarproduct */
/*T_DTYPE operator* (const Vec & v) const {
    T_DTYPE sum = 0;
    for (int i=0; i<T_DIM; i++)
        sum += (this->data[i]) * v[i];
    return sum;
}*/

template<typename T_DTYPE, int T_DIM>
double Vec<T_DTYPE,T_DIM>::norm2( void ) const {
    double tmp = 0;
    for (int i=0; i<T_DIM; i++)
        tmp += (this->data[i]) * (this->data[i]);
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE Vec<T_DTYPE,T_DIM>::product( void ) const {
    T_DTYPE tmp = 1;
    for (int i=0; i<T_DIM; i++)
        tmp *= this->data[i];
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE Vec<T_DTYPE,T_DIM>::sum( void ) const {
    T_DTYPE tmp = 0;
    for (int i=0; i<T_DIM; i++)
        tmp += this->data[i];
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE Vec<T_DTYPE,T_DIM>::min( void ) const {
    T_DTYPE tmp = this->data[0];
    for (int i=0; i<T_DIM; i++)
        if ( this->data[i] < tmp )
            tmp = this->data[i];
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE Vec<T_DTYPE,T_DIM>::max( void ) const {
    T_DTYPE tmp = this->data[0];
    for (int i=0; i<T_DIM; i++)
        if ( this->data[i] > tmp )
            tmp = this->data[i];
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE Vec<T_DTYPE,T_DIM>::mean( void ) const {
    return this->sum() / double(T_DIM);
}

/*template<typename T_DTYPE> // only available for 3D
Vec<T_DTYPE,3> Vec<T_DTYPE,3>::cross(const Vec<T_DTYPE,3> & b) const {
    Vec<T_DTYPE,3> res(0);
    res.data[0] += this->data[1] * b.data[2] - this->data[2] * b.data[1];
    res.data[1] += this->data[2] * b.data[0] - this->data[0] * b.data[2];
    res.data[2] += this->data[0] * b.data[1] - this->data[1] * b.data[0];
    return res;
}*/

/************************ Derived Declarations ********************************/
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator+ (const Vec & v) const {
    Vec res = *this;
    res += v;
    return res;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator-= (const Vec & v) {
    return (*this) += v * (-1);
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator- (const Vec & v) const {
    return (*this)+( v*(-1) );
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator*( const Vec & v ) const {
    Vec res = *this;
    res *= v;
    return res;
}

template<typename T_DTYPE, int T_DIM>
bool Vec<T_DTYPE,T_DIM>::operator!=( const Vec & v ) const {
    return !((*this) == v);
}

template<typename T_DTYPE, int T_DIM>
double Vec<T_DTYPE,T_DIM>::norm() const {
    double res = sqrt((*this).norm2());
    return res;
}

/***************** Broadcasting Operators to save memory **********************/
// this could also be done with , but we want to save memory !
template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator+= (const T_DTYPE a) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] += a;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator-= (const T_DTYPE a) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] -= a;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator*= (const T_DTYPE a) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] *= a;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator/= (const T_DTYPE a) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] /= a;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator+ (const T_DTYPE a) const {
    Vec tmp( *this );
    tmp += a;
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator- (const T_DTYPE a) const {
    Vec tmp( *this );
    tmp -= a;
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator* (const T_DTYPE a) const {
    Vec tmp( *this );
    tmp *= a;
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator/ (const T_DTYPE a) const {
    Vec tmp( *this );
    tmp /= a;
    return tmp;
}


/************************** global overloading ********************************/

/* overload global operator* to allow 3*Vec additionally to Vec*3 */
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator+( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside ) {
    return righthandside + (T_DTYPE) scalar;
}
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator-( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside ) {
    return Vec<T_DTYPE,T_DIM>(scalar) - righthandside;
}
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator*( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside ) {
    return righthandside * (T_DTYPE) scalar;
}
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator/( const T scalar, const Vec<T_DTYPE,T_DIM> & righthandside ) {
    return Vec<T_DTYPE,T_DIM>(scalar) / righthandside;
}

/* Enables cout << Vec<int,2>(1); This also works with fstream and therefore with tout */
template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out, const Vec<T_DTYPE,T_DIM> v ) {
    out << "(";
    for (int i=0; i<T_DIM-1; i++)
        out << v[i] << ",";
    out << v[T_DIM-1] << ")";
    return out;
}
