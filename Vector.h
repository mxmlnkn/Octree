#include <cassert>
#include <iostream>

using namespace std;

#pragma once

#define DEBUG_VECTOR 1


template<typename T_DTYPE, int T_DIM>
class Vec {
public:
    static const int dim = T_DIM;
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
    template<typename T_ETYPE> Vec& operator*=(const Vec<T_ETYPE,T_DIM> & v);
    template<typename T_ETYPE> Vec& operator/=(const Vec<T_ETYPE, T_DIM> & v);
    template<typename T_ETYPE> Vec  operator/ (const Vec<T_ETYPE, T_DIM> & v) const;

    bool operator== (const Vec & v) const;
    bool operator<  (const Vec & v) const;
    bool operator>  (const Vec & v) const;
    bool operator<= (const Vec & v) const;
    bool operator>= (const Vec & v) const;

    double norm2( void ) const;
    T_DTYPE product( void ) const;
    T_DTYPE sum( void ) const;

    template<typename T_ETYPE> Vec  operator+ (const Vec<T_ETYPE, T_DIM> & v) const;
    template<typename T_ETYPE> Vec& operator-=(const Vec<T_ETYPE, T_DIM> & v);
    template<typename T_ETYPE> Vec  operator- (const Vec<T_ETYPE, T_DIM> & v) const;
    template<typename T_ETYPE> Vec  operator* (const T_ETYPE a) const;
    template<typename T_ETYPE> Vec  operator* (const Vec<T_ETYPE,T_DIM> & v) const;
    template<typename T_ETYPE> Vec& operator*=(const T_ETYPE a);
    template<typename T_ETYPE> Vec& operator/=(const T_ETYPE a);
    template<typename T_ETYPE> Vec  operator/ (const T_ETYPE a) const;

    bool operator!= (const Vec & v) const;
    double norm() const;
    void Print( void ) const;
};





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

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator=( const Vec<T_DTYPE,T_DIM> & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = v.data[i];
    return *this;
}

/****************************** Constructors ******************************/

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>::Vec( const T_ETYPE v[T_DIM] ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = v[i];
}

template<typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM>::Vec(const T_DTYPE a) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = a;
}

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>::Vec( const Vec<T_ETYPE,T_DIM> & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = v.data[i]; // implicit conversion from T_ETYPE to T_DTYPE
}

/**************************** Access Operators ****************************/
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
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator=( const Vec<T_ETYPE, T_DIM> & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] = v[i];
    return *this;
}

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator+=( const Vec<T_ETYPE, T_DIM> & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] += v[i];
    return *this;
}

/* Elementwise Multiplication with another Vec of double or int */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator*=( const Vec<T_ETYPE,T_DIM> & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] *= v[i];
    return *this;
}

/* Elementwise Division with another Vec of double or int */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator/=( const Vec<T_ETYPE, T_DIM> & v ) {
    for (int i=0; i<T_DIM; i++)
        this->data[i] /= v[i];
    return *this;
}

/* Elementwise Division with another Vec of double or int */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator/( const Vec<T_ETYPE, T_DIM> & v ) const {
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
    double sum = 0;
    for (int i=0; i<T_DIM; i++)
        sum += (this->data[i]) * (this->data[i]);
    return sum;
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
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator+ (const Vec<T_ETYPE, T_DIM> & v) const {
    Vec res = *this;
    res += v;
    return res;
}

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator-= (const Vec<T_ETYPE, T_DIM> & v) {
    return (*this) += v * (-1);
}

template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator- (const Vec<T_ETYPE, T_DIM> & v) const {
    return (*this)+( v*(-1) );
}

/* Elementwise Multiplication, derived from *= */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator*( const T_ETYPE a ) const {
    Vec res = *this;
    res *= a;
    return res;
}

/* Broadcasting Multiplication, derived from *= */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator*( const Vec<T_ETYPE,T_DIM> & v ) const {
    Vec res = *this;
    res *= v;
    return res;
}

/* Broadcasting Multiplication, derived from elementwise Mul */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator*=( const T_ETYPE a ) {
    (*this) *= Vec(a);
    return *this;
}
/* Broadcasting division from double or int, derived from elementwise div */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM>& Vec<T_DTYPE,T_DIM>::operator/=( const T_ETYPE a ) {
    (*this) /= a;
    return *this;
}
/* Broadcasting division from double or int, derived from elementwise mul */
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
Vec<T_DTYPE,T_DIM> Vec<T_DTYPE,T_DIM>::operator/( const T_ETYPE a ) const {
    return (*this)*(Vec(1.0)/a);
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

/***************************** Debug methods **********************************/
#if DEBUG_VECTOR >= 1
template<typename T_DTYPE, int T_DIM>
void Vec<T_DTYPE,T_DIM>::Print( void ) const {
    cout << "(";
    for (int i=0; i<T_DIM-1; i++)
        cout << this->data[i] << ",";
    cout << this->data[T_DIM-1] << ")";
}
#endif

/* overload global operator* to allow 3*Vec additionally to Vec*3 */
/*template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator+( const T & scalar, const Vec<T_DTYPE,T_DIM> & righthandside )
{
    // scalar multiplication is commutative: s M = M s
    return righthandside + scalar;
}
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator-( const T & scalar, const Vec<T_DTYPE,T_DIM> & righthandside ) {
    return righthandside - scalar;
}
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator*( const T & scalar, const Vec<T_DTYPE,T_DIM> & righthandside ) {
    return righthandside * scalar;
}
template<typename T, typename T_DTYPE, int T_DIM>
Vec<T_DTYPE,T_DIM> operator/( const T & scalar, const Vec<T_DTYPE,T_DIM> & righthandside ) {
    return righthandside / scalar;
} */
/*
ostream& operator <<(ostream& osObject, const storageRentals& rentals)
{

  for (int count = 0; count < 8; count++) {
      osObject << "Unit: " << count + 1 << "    " << rentals.stoUnits[count] << endl;
  }
  return osObject;
}

If the stoUnits member is private you need to make the stream function a friend of your storage class.

friend ostream& operator<<(ostream& osObject, const storageRentals& rentals);*/