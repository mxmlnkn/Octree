/*
 rm shapes.exe; g++ shapes.cpp -o shapes.exe -Wall -Wextra -Wchar-subscripts -Wcomment -Wconversion -Wformat -Wmissing-braces -Wmissing-noreturn -Wparentheses -Wreturn-type -Wshadow -Wsign-compare -Wstrict-aliasing -Wuninitialized -Wunknown-pragmas -Wunreachable-code -std=c++0x -I ../libs/ -g 2>&1 | grep 'error'; ./shapes.exe
*/

#include <iostream>
#include <algorithm> // count
#include <stack>
#include <vector>
#include <string>   // std::stoi, std::stod
#include <sstream>
#include "CompileTime.h"
#include "math/Matrix.h"

using namespace std;

//typedef BaseMatrix<int,2> Matrix;

/******************************************************************************
 * n  ... specifies the order of the shape:                                   *
 *            n=0: Point                                                      *
 *            n=1: Cloud in Cell (CIC)                                        *
 *            n=2: Triangular Shaped Cloud (TSC)                              *
 *            n=3: PSQ                                                        *
 * dx ... cell size, or a scaling factor, which is equal to the diameter of a *
 *        CIC macro particle in e.g. meters                                   *
 *   -> The shape has a cloud diameter = n * dx                               *
 * D  ... distance between the centers of the two radial sym. macro particles *
 *   -> for D >= 2n*dx the charges are not intersecting, therefore the        *
 *      resulting force is the normal Coulomb force! (see proof for radial    *
 *      symmetric charge distributions in my bachelor thesis)                 *
 * d  ... = D / dx is a unitless distance                                     *
 *   -> for d >= 2n the shapes are not intersecting                           *
 * xi ... = d / 2 is another unitless distance used                           *
 * j  ... = ceil(4n-2d) where 0 <= j <= 4n is an integer specifying one case  *
 *          differentiation when integrating.                                 *
 *   -> 2n - j/2 <= d <= 2n - j/2 + 1/2 where d = 2*xi = D/dx meaning for     *
 *         j=4n: 0 <= d <= 1, meaning both shapes overlap almost completely   *
 *         j=0 : 2n <= d <= 2n + 1/4, meaning the shapes are not intersecting *
 ******************************************************************************/

#if 1==0
MathMatrix<double> f1( int n, int j, int np ) {};
MathMatrix<double> f2( int n, int j, int np );
MathMatrix<double> f2Orange( int n, int j, int np ) {};
MathMatrix<double> f2Green ( int n, int j, int np ) {};
MathMatrix<double> f2Violet( int n, int j, int np ) {};
MathMatrix<double> f3( int n, int j, int np ) {};
MathMatrix<double> f3Brown ( int n, int j, int np ) {};
MathMatrix<double> f3Violet( int n, int j, int np ) {};

MathMatrix<double> f( int n, int j )
{
/* internal vector length i.e. polynomial space size. Meaning only powers     *
 * from x**0, x**0, ... x**(np-1) can be represented. If the formula actually *
 * includes higher power, then the result will be wrong. One can estimate the *
 * highest power from the analytical integral or you can test if higher np    *
 * change the result.                                                         *
 * Also the result has certain symmetries and simplifications and             *
 * redundancies which can be used to gauge whether it is correct              */
    const int np = 2*(n+2);
/* np-Vector containing all powers of xi from xi**0 to xi**(np-1)             */
    MathMatrix<double> xi(np,1);

    if ( j <= 2*n ) {
        return f1(n,j,np);
    } else if ( j <= 3*n ) {
        return f2(n,j,np);
    } else {
        return f3(n,j,np);
    }
}

MathMatrix<double> f2( int n, int j, int np ) {
   return  f2Orange(n,j,np) + f2Green(n,j,np) + f2Violet(n,j,np);
}
#endif

#if 1==0
/* recursively call -> may not be ... may be possible for operator << (call recursive function from operator<< */
template<typename T, int T_LEVELS>
void printSpaceMathMatrix( T mat, Vec<int,T_LEVELS> index )
{
    int size = mat.getSize()[0];
    for ( int i = 0; i < size; ++i )

    int size;
    size[0] = mat.getSize()[0];
    size[1] = mat[0].getSize()[0];

    for ( int i = 0; i < nestedlevels; ++i )
        std::cout << " " << spacevarnames[i] << "_{";

    BaseMatrix<int,1> index(nestedlevels);
    for ( int i = 0; i < size.product(); ++i )
    {
        /* recursively unpack mat which is Mathmatrix<Mathmatrix<...>...> *
         * using a local function i.e. struct lambda ...                  */
        struct_recunpack recunpack;
        std::cout << recunpack( mat, index );

        /* increment vector index */
        int closedBracktes = 0;
        index[nestedlevels-1] += 1;
        for ( int j = nestedlevels-1; j >= 0; --j ) {
            if ( index[j] >= size[j] ) {
                if ( j > 0 ) {
                    index[j] = 0;
                    index[j-1] += 1;
                } else
                    assert( i == size.product()-1 );
                closedBracktes++;
                std::cout << "}_" << spacevarnames[nestedlevels - closedBracktes]
                          << " ";
            }
        }
        if ( i != size.product()-1 ) {
            std::cout << ",";
            for ( int j = 0; j < closedBracktes; ++j )
                std::cout << " " << spacevarnames[nestedlevels - closedBracktes]
                          << "_{";
        }
    }
}
#endif

/* Print in another way */

/* Last specialized call for recursion toPolynomial goes here */
template<typename T>
std::string toPolynomial( T v, int n0, int maxlevel, int level = 0)
{
    /* first output doesn't have to have a plus sign */
    std::stringstream out;
    out << std::showpos << v << std::noshowpos;
    return out.str();
}

/* Print the polynomial representation of given v, e.g. {0,1,2} -> x+2x**2 */
template<typename T>
std::string toPolynomial( MathMatrix<T> v, int n0, int maxlevel, int level = 0)
{
    const char varnames[4][10] = { "d", "x", "x1" };
    std::stringstream out;

    bool firstprinted = false;
    for ( int i = 0; i < v.size[0]; ++i )
    {
        /* don't print summands whose coefficients are 0 */
        if ( v[i] == v[i]*0 )
            continue;

        std::string tmp = toPolynomial( v[i], n0, maxlevel, level+1 );
        /* count coefficients by counting + and - */
        int ncoeffs = 0;
        ncoeffs += std::count( tmp.begin(), tmp.end(), '+' );
        ncoeffs += std::count( tmp.begin(), tmp.end(), '-' );

        /* only print a sign, if needed i.e. if tmp does not begin with one   *
         * and if this is not the very first string to be printed             */
        if ( tmp[0] != '+' and tmp[0] != '-' and firstprinted )
            out << "+";
        firstprinted = true;

        /* if i-n0 > 0 and return of toPolynomial == "1", then don't print it*/
        if ( not ( (tmp.compare( std::string("+1") ) == 0) and (i-n0 > 0) ) )
            if ( ncoeffs > 1 ) {
                out << "(" << tmp << ")";
            } else {
                out << tmp;
            }

        /* print variable name: x**(n) if not x**0 */
        if ( i-n0 != 0 )
        {
            /* print negative powers as fractions */
            if ( i-n0 < 0 )
                out << "/";
            out << varnames[level];
            /* output power only if necessary */
            if ( i-n0 < 0 )
                out << "**" << -(i-n0);
            else if ( i-n0 != 1 )
                out << "**" << i-n0;
        }
    }

    return out.str();
}

template<typename T>
MathMatrix<T> & fromString
( std::string sr, std::vector<std::string> vars, MathMatrix<T> & result );

/* Last specialized call for recursion toMultiplyer goes here */
template<typename T>
T toMultiplyer( T v, int n0 ) {
    return v;
}

/* convert vector to multiplication operator i.e. matrix */
template<typename T>
MathMatrix<T> toMultiplyer( MathMatrix<T> v, int n0 )
{
    assert( v.isVector() );
    int np = v.getVectorDim();

    /* set up result matrix with correct dimension */
    MathMatrix<T> mul(np,np);
    for ( int i = 0; i < mul.getSize().product(); ++i )
       mul[i].setSize(np,np);
    mul = 0; /* setSize only allocates but does not set to 0 */

    for ( int i = 0; i < v.getSize()[0]; ++i ) {
        /* recursive convert inner vectors to multipliers */
        /* for division set upper diagonal, for multipl. set lower diagonal */
        mul.setDiagonal( toMultiplyer( v[i] ), n0-i );
    }
}

int main( void ) {
    /* the main idea is like this: 1+x = \begin{pmatrix} 0 \\ 1 \\ 2 \\ 0 \\ 0 \end{pmatrix} \cdot \begin{pmatrix} x^{-1} \\ x^0 \\ x^1 \\ x^2 \\ x^3 \end{pmatrix}\\
(1+x)\cdot \hat{=} \begin{pmatrix} 1 & 0 &0 & 0 &0 \\ 0 & 1 &0 & 0 &0 \\ 0 & 0 &1 & 0 &0 \\ 0 & 0 &0 & 1 &0 \\0 & 0 &0 & 0 &1 \\ \end{pmatrix} + \begin{pmatrix} 0 & 0 &0 & 0 &0 \\ 2 & 0 &0 & 0 &0 \\ 0 & 2 &0 & 0 &0 \\ 0 & 0 &2 & 0 &0 \\0 & 0 &0 & 2 &0 \\ \end{pmatrix}\\
\Rightarrow (1+x)\cdot(1+x) = \begin{pmatrix} x^{-1} \\ x^0 \\ x^1 \\ x^2 \\ x^3 \end{pmatrix} \begin{pmatrix} 1 & 0 &0 & 0 &0 \\ 2 & 1 &0 & 0 &0 \\ 0 & 2 &1 & 0 &0 \\ 0 & 0 &2 & 1 &0 \\0 & 0 &0 & 2 &1 \end{pmatrix}\begin{pmatrix} 0 \\ 1 \\ 2 \\ 0 \\ 0 \end{pmatrix} = \begin{pmatrix} x^{-1} \\ x^0 \\ x^1 \\ x^2 \\ x^3 \end{pmatrix} \begin{pmatrix} 0 \\ 1 \\ 4 \\ 4 \\ 0 \end{pmatrix} = 4x^2+4x^1+1 */

    /* ( d**2 - (x1-x2)**2 ) / ( x1*x2 ) =
         ( d**2*x1**(-1) - x1) * x2**(-1)
         +                   2 * x2**( 0)
         +         -1*x1**(-1) * x2**( 1) =
      ( ( d**2, 0, -1, 0, 0 )_x1,
        ( 0   , 2, 0 , 0, 0 )_x1,
        (   -1, 0, 0 , 0, 0 )_x1,
        (    0, 0, 0 , 0, 0 )_x1,
        (    0, 0, 0 , 0, 0 )_x1
      )_x2 = again polynomial space for d! */

    int np = 5;
    int n0 = 1;
    MathMatrix<double> spaced(np,1); spaced = 0;
/*
    MathMatrix<double> spacedconst(np,1); spacedconst = 0;
    MathMatrix<MathMatrix<double> > spacex1(np,1); spacex1 = spaced;
    MathMatrix<MathMatrix<MathMatrix<double> > > spacex2(np,1); spacex2 = spacex1;

    spacex2[n0-1][n0-1][n0+2] =  1;
    spacex2[n0-1][n0+1][n0  ] = -1;
    spacex2[n0  ][n0  ][n0  ] =  2;
    spacex2[n0+1][n0-1][n0  ] = -1;

    MathMatrix<double> MI(np,np);
    for ( int i = 0; i < np-1; ++i )
        MI[i+1,i] = 1.0 / double(i+1-n0); // it is wanted that 1/0 is INF
    MI[np-1,np-1] = INF;
    MathMatrix<double> Mx(np,np);
    for ( int i = 0; i < np-1; ++i )
        Mx[i+1,i] = 1;

    #define SCP(a,b) (((a).Transpose())*b)[0]
    MathMatrix<MathMatrix<double> > x2e;
    MathMatrix<MathMatrix<MathMatrix<double> > > x2evec;
    x2evec[n0][n0][n0] = 1; // !!! ToDo: alle n0 sind 1
    for ( int i = n0 + 1; i < x2evec.dim; ++i )
        x2evec[i] = x2evec[i-1] * x2e;

    //MI*Mx*Mx*SCP(MI*Mx*Mx*spacex2,x2evec);

    std::cout << "\frac{Q_1 Q_2}{4\pi\epsilon_0 d**2} (4 pi R**2)**2 R**2  rho0_n(R)^2\n";
*/

/****************************** Expand (1-2x)**2 ******************************/
    MathMatrix<int> x1(np,1);
    x1 = 0;
    std::cout << "x1: " << x1 << "\n";
    x1[1] =  1;
    x1[2] = -2;
    std::cout << "x1: " << x1 << "\n";

    MathMatrix<int> mx1(np,np);
    mx1 = 3;
    std::cout << "mx1: " << mx1 << "\n";
    mx1.setDiagonal(17,-2);
    std::cout << "mx1: " << mx1 << "\n";
    mx1 = 0;
    for ( int i = 0; i < x1.size[0]; ++i ) {
        mx1.setDiagonal( x1[i], n0-i );
        // for division set upper diagonal, for multiplication set lower diagonal
    }
    std::cout << "mx1: " << mx1 << "\n";
    std::cout << "mx1*x1: " << (mx1*x1) << "\n";
    std::cout << "  => [ ";

    {bool firstprinted = false;
    for ( int i = 0; i < x1.size[0]; ++i ) {
        if ( x1[i] == 0 )
            continue;
        if ( firstprinted )
            std::cout << " + ";
        std::cout << x1[i] << " x**(" << i-n0 << ")";
        firstprinted = true;
    }}

    std::cout << " ]**2 = ";

    MathMatrix<int> res = mx1*x1;
    //res = mx1*x1;

    {bool firstprinted = false;
    for ( int i = 0; i < res.size[0]; ++i ) {
        if ( res[i] == 0 )
            continue;
        if ( firstprinted )
            std::cout << " + ";
        std::cout << res[i] << " x**(" << i-n0 << ")";
        firstprinted = true;
    }}

    std::cout << "\n";

/****************************** Expand (d-3x)**2 ******************************/

    /* initialize polynomial: d-2x */
    MathMatrix<int> dspace(np,1);
        dspace = 0;
    MathMatrix<MathMatrix<int> > vxd(np,1);
        vxd = dspace;
    dspace[2] =  1;
    vxd[1]    =  dspace;
    dspace[1] = -3;
    dspace[2] =  0;
    vxd[2]    =  dspace;

    /* set up multiplier manually: (d-3x)* */
    MathMatrix<int> mzero(np,np);
    mzero = 0;
    MathMatrix<int> munit(np,np);
    munit = 0;
    munit.setDiagonal(1);
    MathMatrix<int> mmul (np,np);
    mmul  = 0;
    mmul.setDiagonal(1,-1);

    MathMatrix<MathMatrix<int> > mx(np,np);
    mx = mzero;
    mx.setDiagonal( -3 * munit,-1);
    MathMatrix<MathMatrix<int> > md(np,np);
    md = mzero;
    md.setDiagonal(mmul,0);

    /* convert vector to multiplication operator i.e. matrix */
    MathMatrix<MathMatrix<int> > mxd(np,np);
    for ( int i = 0; i< mxd.getSize().product(); ++i )
       mxd[i].setSize(np,np);
    mxd = 0;

    for ( int i = 0; i < vxd.getSize()[0]; ++i ) {
        /* recursive convert inner vectors to multipliers */
        MathMatrix<int> md(np,np);
        md = vxd[0][0] * 0;
        for ( int j = 0; j < vxd[i].getSize()[0]; ++j )
            md.setDiagonal( vxd[i][j], n0-j );

        mxd.setDiagonal( md, n0-i );
        // for division set upper diagonal, for multiplication set lower diagonal
    }
    /*********************************************************/

    auto sum = md+mx;
    std::cout << "md+mx: " << sum << "\n\n";
    std::cout << "mxd: " << mxd << "\n\n";
    std::cout << "md+mx == mxd: " << (mxd == sum) << "\n\n";

    std::cout << vxd << "\n";
    std::cout << "(md+mx)*vxd:\n";
    std::cout << "    [" << toPolynomial( vxd, n0, 1 );
    std::cout << "]**2\n  = " << toPolynomial( sum*vxd, n0, 1 );

    /************************* Integrate Coulomb Force ************************/
    /* Integrate[ x2**2
         Integrate[ x1**2 * ( d*d-(x1-x2)**2 )/( 4 x1 x2 ), {x1,0,1} ]
       , {x2,0,1} ]
       => we are working 3(!) different spaces: x1,x2,d
          After inserting the limits for x1, we leave that space!
          Then we leave the polynomial space for x2!
          At last we need to output the coefficients for d
       => the limits for x1 is a d-coeff-vector of x2-coeff-vectors of
          x1-coeff-vectors. If d and x2 don't appear in the limits, then
          the two outer vectors are only non-zero on their diagonals
     */
    typedef MathMatrix<MathMatrix<MathMatrix<int>>> SpaceX1;
    /* Todo: try to convert ShapesV6.nb into this C++ scheme, maybe even with sed or awk! */

    //std::string sr1( "(x1)**(2) * ( (((d)**(2))) - (x1-x2)**2 )" ); // "/( 4 x1 x2 )"
    std::string sr1( "x1**2 * ( d**2 - x1**2 - x2**2 + 2 x1 x2 )" );
    SpaceX1 result;

    std::vector<std::string> varnames;
    varnames.push_back("d");
    varnames.push_back("x1");
    varnames.push_back("x2");
    std::cout << "\n\n";
    fromString( "((x1))", varnames, result );
    std::cout << "\n\n";
    fromString( "((x1)+x2)", varnames, result );

    return 0;
}

template<typename T>
MathMatrix<T> & fromString
( std::string str, std::vector<std::string> varnames, MathMatrix<T> & result )
{
    std::cout << "Input String   : " << str << "\n";
    std::cout << "Variable Names : ";
    for ( int i = 0; i < varnames.size(); ++i ) {
        std::cout << varnames[i] << ",";
    }
    std::cout << "\n";

    /* strip whitespaces */
    {int i = 0;
    while ( i < str.length() ) {
        if ( str[i] == ' ' )
            str.erase(i,1);
        else
            ++i;
    }}
    std::cout << "Strip Whitespaces: " << str << "\n";

    /* strip unnecessary parentheses, e.g. '(a+b)' or '(a)+b' reduce to 'a+b' */
    {assert(    std::count( str.begin(), str.end(), '(' )
            == std::count( str.begin(), str.end(), ')' ) );
    std::stack<int> ppos;
    bool pmsignfound = false;
    for ( int i = 0; i < str.length(); ++i ) {
        switch ( str[i] ) {
            case '(':
                ppos.push(i);
                break;
            case ')':
                if ( not pmsignfound ) {
                    str.erase( i,1 );
                    str.erase( ppos.top(), 1 );
                    ppos.pop();
                    i -= 2;
                }
                pmsignfound = false;
                break;
            case '+':
            case '-':
                pmsignfound = true;
                break;
        }
    }
    while ( str[0] == '(' and str[ str.length()-1 ] == ')' ) {
        str.erase( str.length()-1, 1 );
        str.erase( 0, 1 );
    }}
    std::cout << "Strip unnecessary parentheses: " << str << "\n";

    /* analyze highest level only, meaning: 2a*(...)+b-(...) call recursively
     * for (...) */
    result = 0;
    MathMatrix<T> summand = result;
    std::stack<int> ppos;
    bool summandinitialized = true;
    for ( int i = 0; i < str.length(); ++i )
    {
        /* Recursively evaluate things in parentheses */
        if ( str[i] == '(' ) {
            ppos.push(i);
        } else if ( str[i] == ')' ) {
            if ( ppos.size() == 1 ) {
                std::string sexpr = str.substr( ppos.top()+1, i-1 );
                std::cout << "Recursively evaluate: " << sexpr << "\n";
            }
            ppos.pop();
        }

        MathMatrix<T> tmp = summand*0;

        /* Evaluate numbers appearing */
        if ( str[i] >= '0' and str[i] <= '9' ) {
            int j = i;
            int dotsfound = 0;
            while ( ( str[j] >= '0' and str[j] <= '9' ) or  str[j] == '.' ) {
                ++j;
                if ( str[j] == '.' )
                    ++dotsfound;
            }
            std::string number = str.substr( i, j-i );
            assert( dotsfound <= 0 );

            if ( dotsfound == 0 ) {
                tmp = std::stoi( number );
            } else if ( dotsfound == 1 ) {
                tmp = std::stod( number );
            }
        }

        /* variables may only begin with a letter */
        if ( ( str[i] >= 'A' and str[i] <= 'Z' ) or
             ( str[i] >= 'a' and str[i] <= 'z' ) )
        {

        }

        if ( summandinitialized ) {
            summand = tmp;
            summandinitialized = false;
        } else {
            summand =  * tmp;
            summandinitialized = false;
        }

        if ( str[i] == '+' ) {
            result += summand;
            summand = 0;
            summandinitialized = true;
        } else if ( str[i] == '-' ) {
            result += summand;
            //summand = toMultiplyer( -1 );
            summandinitialized = false;
        }
    }

    #if 1==0
    int foundpos = std::string::npos; /* stores last found position */
    /* what datatype do i even return -.- ?... can't specify auto return type and oh my god is this frustratingly difficult xD... -> need class which manages recursion in a non-recursive way! instead of using lowlevel matrices! */
    do {
        /* find powers */
        std::string substr = varnames[0] + std::string("**");
        foundpos = str.find( substr );
        std::cout << "Found " << varnames[0] << " to the power of "
                  << str[ foundpos + substr.length() ] << " with coefficients:";
        #if 1==0
        /* find coefficient string */
        int parentheses = 0;
        for ( int i = foundpos-1; i >= 0; --i ) {
            /* if not inside a pair of parentheses */
            if ( parentheses = 0 ) {
                if ( str[i] == '+' or str[i] == '-' )
            }
            if ( str[i] !=
        }
        #endif
    }
    while ( foundpos != std::string::npos );
    #endif

    return result;
}

