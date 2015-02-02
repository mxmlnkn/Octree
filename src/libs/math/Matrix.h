#pragma once

#include <fstream>
#include <algorithm> // std::min
#include "math/TVector.h"

struct MDIM {
    int m,n;
};

class Matrix {
    private:
        int m,n;
        double **data;
        void Setup(int m, int n);                   //Allocates Memory for Matrix.data and stores m,n
        void Clear(void);
    
    public:
        Matrix(const Matrix &mat);                  //Copy method/Constructor
        Matrix(int m, int n);                       //Create new 2d-array with (m,n) being the dimension/Constructor
        Matrix(double x);
        Matrix(char* filename);                     //Automatically read from file
        template<int T_DIM> Matrix(const Vec<double,T_DIM> v);
        ~Matrix(void);                              //delete matrix from heap
        
        template<int T_DIM> operator Vec<double,T_DIM>() const;

        Matrix& operator=(const Matrix &mat);
        Matrix& operator=(const double val);
        template<int T_DIM> Matrix& operator=(const Vec<double,T_DIM> val);
        double& operator()(int i, int j) const;     //Get ptr to element m,n
        double& operator[](int x) const;            //Get ptr to element m,1 or 1,n if n==1 or m==1
        // make indexing with negative indices possible !!! or subspaces (slices,...)
        Matrix operator+(const Matrix &mat);        //Add up to matrices
        Matrix operator*(const Matrix &mat);        //Matrix Multiplication
        Matrix operator-(const Matrix &mat);        //Subtraction
        inline Matrix operator/(double divisor);    //scalar Division based on Multiplication
        Matrix operator*(double factor);            //scalar Multiplication
        bool operator==(const Matrix &mat);

        Matrix operator* (const double a) const;
        
        double Minor(int row, int col) const;       //Counting from 1
        double Det(void) const;
        Matrix Invert(void) const;
        Matrix Adjugate(void) const;
        Matrix Transpose(void) const;
        /* Returns Norm of Matrix if it is a Vector, else error (returns -1) */
        double Norm(void) const;
        /* Returns matrix with only positive elements (aij -> |aij|) */
        Matrix Abs(void) const;     
        int Rank(void) const;
        double Trace(void) const;
        Matrix RowEchelon(void) const;
        inline bool IsSquare(void) const;
        inline bool IsVector(void) const;
        MDIM GetDim(void) const;
        int GetVectorDim(void) const;
        int GetSquareDim(void) const;

        void SetDiagonal(double val=1);
        void SetDiagonal(double val[]);
        void SetAll(double val);
        
        void DelRow(int row, int amount=1);         //Deletes i-th Row counting from 1
        void DelCol(int col, int amount=1);
        Matrix Augment(const Matrix &mat) const;    //conjoines object with mat (if number of rows is equivalent)

        int Load(char* filename);
        void Save(char* filename);
        void Print(void) const;
};

template<typename T>
Matrix operator*( const T scalar, const Matrix & rhs );

#include "Matrix.cpp"
