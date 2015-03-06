#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "math/Matrix.h"
#include <cctype> // isalpha, isalnum, isdigit, islower, isspace, isupper

#define POLYNOMIAL_DEBUG 100

/* the template parameter must be explicitely given and is very important! e.g.
 * MathMatrix<int> or MathMatrix<MathMatrix<int>> and so on. varnames should
 * contain at least as much variables as nested levels here */
template<typename T_POLYNOM>
class Polynomial {
public:
    T_POLYNOM data;

private:
    int np; /* number of powers */
    int n0; /* index where power 0 is being stored */
    std::vector<std::string> varnames;

public:
    Polynomial( std::vector<std::string> pvarnames, std::string pol = "0" );
    Polynomial( std::string pvarnames = "x", std::string pol = "0" );
    int checkVariableNames( void );

    /* Counts recursion MathMatrix */
    int countLevels(void) const;
    /* Print the polynomial representation of given v, e.g. {0,1,2} -> x+2x**2 */
    std::string toString(void) const;
    Polynomial<T_POLYNOM> & fromString( std::string sr );
    Polynomial<T_POLYNOM> & fromOneVarString( std::string sr );
    /* convert vector to multiplication operator i.e. matrix */
    Polynomial<T_POLYNOM> toMultiplyer(void) const;
    Polynomial<T_POLYNOM> & setSize( int pnp );
    Polynomial<T_POLYNOM> & setZero( int pn0 );

private:
    template<typename T_ELEM>
    MathMatrix<MathMatrix<T_ELEM>> & setSize( MathMatrix<MathMatrix<T_ELEM>> & pol, int pnp ) const;
    template<typename T_ELEM>
    MathMatrix<T_ELEM> & setSize( MathMatrix<T_ELEM> & pol, int pnp ) const;

    template<typename T_ELEM>
    int countLevels( MathMatrix<T_ELEM> v, int level = 0 ) const;
    /* Last specialized call for recursion toPolynomial goes here */
    template<typename T_ELEM>
    int countLevels( T_ELEM v, int level = 0 ) const;

    template<typename T_ELEM>
    std::string toString( const MathMatrix<T_ELEM> & v, int level = 0) const;
    /* Last specialized call for recursion toPolynomial goes here */
    template<typename T_ELEM>
    std::string toString( const T_ELEM & v, int level = 0) const;

    template<typename T_ELEM>
    MathMatrix<T_ELEM> & fromOneVarString( std::string str, MathMatrix<T_ELEM> & result, int level = 0 ) const;
    template<typename T_ELEM>
    MathMatrix<MathMatrix<T_ELEM>> & fromOneVarString ( std::string str, MathMatrix<MathMatrix<T_ELEM>> & result, int level = 0 ) const;

    template<typename T_ELEM>
    MathMatrix<T_ELEM> & fromString ( std::string sr, MathMatrix<T_ELEM> & result ) const;

    template<typename T_ELEM>
    MathMatrix<MathMatrix<T_ELEM>> toMultiplyer( const MathMatrix<MathMatrix<T_ELEM>> & v ) const;
    template<typename T_ELEM>
    MathMatrix<T_ELEM> toMultiplyer( const MathMatrix<T_ELEM> & v ) const;
};

typedef Polynomial< MathMatrix<int> > Pol1Vars;
typedef Polynomial< MathMatrix<MathMatrix<int>> > Pol2Vars;
typedef Polynomial< MathMatrix<MathMatrix<MathMatrix<int>>> > Pol3Vars;

template<typename T_POLYNOM>
std::ostream& operator<<(std::ostream& out, const Polynomial<T_POLYNOM>& data);

#include "Polynomial.tpp"
