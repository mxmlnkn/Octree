#pragma once

#include <fstream>
#include <algorithm> // std::min
#include <cstdlib>   // malloc
#include "math/TBaseMatrix.h"
#include "math/TVector.h"

template<typename T_DTYPE>
class MathMatrix : public BaseMatrix<T_DTYPE,2> {
public:
    typedef Vec<int,2> VecI;

    MathMatrix( void );
    MathMatrix( VecI psize );
    MathMatrix( int m, int n );
    MathMatrix( const MathMatrix & m );
    //Matrix(char* filename);                     //Automatically read from file
    ~MathMatrix( void );
    MathMatrix & operator=(const MathMatrix m); // if I make this a call by reference, then it won't be recognized as a specialization of the next assignment and thereby breaking everything, but copy-arguments could be memory intensive
    template<typename T_ETYPE> MathMatrix & operator=(T_ETYPE a); // broadcast value to all elements

    /* Conversion Operators */
    // template<int T_DIM> operator Vec<double,T_DIM>() const;
    // template<int T_DIM> operator Vec<double,T_DIM>() const;
    // template<int T_DIM> Matrix& operator=(const Vec<double,T_DIM> val);
    // template<int T_DIM> Matrix(const Vec<double,T_DIM> v);

    MathMatrix & operator+=(const MathMatrix &mat);    // Add up to matrices
    MathMatrix & operator*=(const MathMatrix &mat);    // Matrix Multiplication
    MathMatrix & operator-=(const MathMatrix &mat);    // Subtraction
    template<typename T_ETYPE>
    inline MathMatrix & operator/=(T_ETYPE a);         // scalar Division based on Multiplication
    template<typename T_ETYPE>
    MathMatrix & operator*=(T_ETYPE a);                // scalar Multiplication

    MathMatrix operator+(const MathMatrix &mat) const; // Add up to matrices
    MathMatrix operator*(const MathMatrix &mat) const; // Matrix Multiplication
    MathMatrix operator-(const MathMatrix &mat) const; // Subtraction
    template<typename T_ETYPE>
    inline MathMatrix operator/(T_ETYPE a) const;       // scalar Division based on Multiplication
    template<typename T_ETYPE>
    MathMatrix operator*(T_ETYPE a) const;              // scalar Multiplication
    bool operator==(const MathMatrix &mat) const;
    inline bool operator!=(const MathMatrix &mat) const;

    // /* Returns Norm of Matrix if it is a Vector, else error (returns -1) */
    /* most of these functions only make sense if T_DTYPE is double! -> make another subclass? */
    // double norm(void) const;
    // double (minor)(int row, int col) const;             // Counting from 1. clash with minor-macro -.-
    // double det(void) const;
    // MathMatrix invert(void) const;
    // MathMatrix adjugate(void) const;
    // /* Returns matrix with only positive elements (aij -> |aij|) */
    // MathMatrix abs(void) const;
    // int rank(void) const;
    // MathMatrix rowEchelon(void) const;
    MathMatrix transpose(void) const;
    T_DTYPE trace(void) const;

    inline bool isSquare(void) const;
    inline bool isVector(void) const;

    T_DTYPE & operator()( int i, int j );     // Get reference to element
    T_DTYPE   operator()( int i, int j ) const;     // Just read element

    inline int getVectorDim(void) const;
    inline int getSquareDim(void) const;

    MathMatrix & setDiagonal(T_DTYPE a=1, int k = 0);
    MathMatrix & setDiagonal(T_DTYPE a[], int k = 0);
    template<typename T_ETYPE> inline MathMatrix & setAll(T_DTYPE a);
    MathMatrix & setSize(int m, int n);

    // deletes i-th Row counting from 0 !!!!!!! was changed from counting from 1 => adjust all which use this: !!!minor!!!
    MathMatrix & delRow(int irow = 0, int nrows = 1);
    MathMatrix & delCol(int icol = 0, int ncols = 1);
    //MathMatrix augment(const MathMatrix &mat) const;    // merges two matrices (if number of rows is equivalent)

    // int Load(char* filename);
    // void Save(char* filename);

    // explicit operator Vec(); // convert to vector class
};

template<typename T_DTYPE, typename T_ETYPE>
MathMatrix<T_DTYPE> operator*( const T_ETYPE a, const MathMatrix<T_DTYPE> & rhs );

#include "Matrix.cpp"
