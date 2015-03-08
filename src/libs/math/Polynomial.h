#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "math/Matrix.h"
#include <cctype> // isalpha, isalnum, isdigit, islower, isspace, isupper

#define POLYNOMIAL_DEBUG 98
#define MAX_POWERS 5

/* the template parameter must be explicitely given and is very important! e.g.
 * MathMatrix<int> or MathMatrix<MathMatrix<int>> and so on. varnames should
 * contain at least as much variables as nested levels here */
template<typename T_COEFF>
class Polynomial {
public:
    T_COEFF data;
    std::vector<std::string> varnames;

private:
    int np; /* number of powers */
    int n0; /* index where power 0 is being stored */

    std::vector<std::string> splitByDelimiter( std::string src, char delimiter = -1 ) const;

public:
    /* first variable name in list or string will be the most outward variable*/
    Polynomial( std::vector<std::string> pvarnames, std::string pol = "0", int pn0 = 0 );
    Polynomial( std::string pvarnames = "x", std::string pol = "0", int pn0 = 0 );
    Polynomial( const T_COEFF &, int pn0 = 0 );
    Polynomial( std::vector<std::string> pvarnames, const T_COEFF & src, int pn0 = 0 );
    int checkVariableNames( void );

    /* Counts recursion MathMatrix */
    int countLevels(void) const;
    /* Print the polynomial representation of given v, e.g. {0,1,2} -> x+2x**2 */
    std::string toString(void) const;
    Polynomial<T_COEFF> & fromString( std::string sr );
    Polynomial<T_COEFF> & fromOneVarString( std::string sr );
    /* convert vector to multiplication operator i.e. matrix */
    T_COEFF toMultiplyer(void) const;
    Polynomial<T_COEFF> & setSize( int pnp );
    Polynomial<T_COEFF> & setZero( int pn0 );
    int getSize( void ) const;
    int getZero( void ) const;
    Polynomial<T_COEFF> & reorderVariableNames( std::string to );

    Polynomial<T_COEFF> & operator+=( double b );
    Polynomial<T_COEFF> & operator-=( double b );
    Polynomial<T_COEFF> & operator*=( double b );
    Polynomial<T_COEFF> & operator+=( std::string b );
    Polynomial<T_COEFF> & operator-=( std::string b );
    Polynomial<T_COEFF> & operator*=( std::string b );
    Polynomial<T_COEFF> & operator+=( Polynomial<T_COEFF> );
    Polynomial<T_COEFF> & operator-=( Polynomial<T_COEFF> );
    Polynomial<T_COEFF> & operator*=( Polynomial<T_COEFF> );
    template<typename T_ETYPE>
    Polynomial<T_COEFF> operator+( T_ETYPE b ) const;
    template<typename T_ETYPE>
    Polynomial<T_COEFF> operator-( T_ETYPE b ) const;
    template<typename T_ETYPE>
    Polynomial<T_COEFF> operator*( T_ETYPE b ) const;

    bool operator==( Polynomial<T_COEFF> b );

    void operator=( std::string );
    void operator=( const T_COEFF & );
    operator std::string(void);

    auto integrate( std::string variable, std::string from, std::string to ) const -> Polynomial<decltype(this->data[0])>;

private:
    template<typename T_ELEM>
    static MathMatrix<MathMatrix<T_ELEM>> & setSize
    ( MathMatrix<MathMatrix<T_ELEM>> & pol, int pnx, int pny = 1 );

    template<typename T_ELEM>
    static MathMatrix<T_ELEM> & setSize
    ( MathMatrix<T_ELEM> & pol, int pnx, int pny = 1 );

    template<typename T_ELEM>
    static int countLevels( MathMatrix<T_ELEM> v, int level = 0 );
    /* Last specialized call for recursion toPolynomial goes here */
    template<typename T_ELEM>
    static int countLevels( T_ELEM v, int level = 0 );

    template<typename T_ELEM>
    std::string toString( const MathMatrix<T_ELEM> & v, int level = 0) const;
    /* Last specialized call for recursion toPolynomial goes here */
    template<typename T_ELEM>
    std::string toString( const T_ELEM & v, int level = 0) const;

    /**************************** fromOneVarString ****************************
     * converts simple strings like 'x**2' or '2' into internal polynomial    *
     * format. Won't work with more complex strings like '3x'. Use fromString *
     * for that                                                               *
     **************************************************************************/
    template<typename T_ELEM>
    MathMatrix<T_ELEM> & fromOneVarString
    ( std::string str, MathMatrix<T_ELEM> & result, int level = 0 ) const;

    template<typename T_ELEM>
    MathMatrix<MathMatrix<T_ELEM>> & fromOneVarString
    ( std::string str, MathMatrix<MathMatrix<T_ELEM>> & result, int level = 0 ) const;

    /**************************** fromOneVarString ****************************
     * Parses almost arbitrary complex polynomials like 'x + ( x + 3y**2 )'   *
     * except e.g. (x+y)**3                                                   *
     **************************************************************************/
    template<typename T_ELEM>
    MathMatrix<T_ELEM> & fromString
    ( std::string sr, MathMatrix<T_ELEM> & result ) const;

    /****************************** to Multiplyer *****************************
     * returns a np x np matrix which can be multiplied onto                  *
     * this->data to calculate polynomial multiplication                      *
     **************************************************************************/
    template<typename T_ELEM>
    static MathMatrix<MathMatrix<T_ELEM>> toMultiplyer
    ( const MathMatrix<MathMatrix<T_ELEM>> & v, int n0 );

    template<typename T_ELEM>
    static MathMatrix<T_ELEM> toMultiplyer
    ( const MathMatrix<T_ELEM> & v, int n0 );

};

typedef Polynomial< MathMatrix<double> > Pol1Var;
typedef Polynomial< MathMatrix<MathMatrix<double>> > Pol2Vars;
typedef Polynomial< MathMatrix<MathMatrix<MathMatrix<double>>> > Pol3Vars;

template<typename T_COEFF>
std::ostream& operator<<(std::ostream& out, const Polynomial<T_COEFF>& data);

/* overwrites member operator*, so that this end in an infinite loop and crash */
/* template<typename T_COEFF, typename T_ETYPE>
Polynomial<T_COEFF> operator*( const T_ETYPE a, const Polynomial<T_COEFF> & rhs ); */

#include "Polynomial.tpp"
