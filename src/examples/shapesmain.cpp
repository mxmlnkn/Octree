/*
 rm shapesmain.exe; g++ shapesmain.cpp -o shapesmain.exe -Wall -Wextra -Wchar-subscripts -Wcomment -Wconversion -Wformat -Wmissing-braces -Wmissing-noreturn -Wparentheses -Wreturn-type -Wshadow -Wsign-compare -Wstrict-aliasing -Wuninitialized -Wunknown-pragmas -Wunreachable-code -Wno-unused-parameter -Wno-unused-variable -Werror -std=c++0x -I ../libs/ -g 2>&1 | grep 'error'; ./shapesmain.exe
*/

#include <iostream>
#include <iomanip> //setw
#include <algorithm> // count
#include <stack>
#include <vector>
#include <string>   // std::stoi, std::stod
#include <sstream>
#include "CompileTime.h"
#include "math/Matrix.h"
#include "math/Polynomial.h"
#include "shapes.cpp"

int main( void ) {
    {
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
        MathMatrix<int> mdrec(np,np);
        mdrec = vxd[0][0] * 0;
        for ( int j = 0; j < vxd[i].getSize()[0]; ++j )
            mdrec.setDiagonal( vxd[i][j], n0-j );

        mxd.setDiagonal( mdrec, n0-i );
        // for division set upper diagonal, for multiplication set lower diagonal
    }
    /*********************************************************/

    auto sum = md+mx;
    std::cout << "md+mx: " << sum << "\n\n";
    std::cout << "mxd: " << mxd << "\n\n";
    std::cout << "md+mx == mxd: " << (mxd == sum) << "\n\n";

    /*std::cout << vxd << "\n";
    std::cout << "(md+mx)*vxd:\n";
    std::cout << "    [" << toPolynomial( vxd, n0, 1 );
    std::cout << "]**2\n  = " << toPolynomial( sum*vxd, n0, 1 );*/
    }

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

    std::vector<std::string> varnames;
    varnames.push_back("d");
    varnames.push_back("x1");
    varnames.push_back("x2");

    Pol3Vars result ( "d|x1|x2" );
    Pol2Vars result2( varnames  );
    Pol1Var  result1( varnames  );

    std::cout << "\n\n";

    std::cout << "3 -> " << result1.fromOneVarString("3") << "\n" << "raw: " << result1.data << "\n";
    std::cout << "3 -> " << result2.fromOneVarString("3") << "\n" << "raw: " << result2.data << "\n";
    std::cout << "((x1)+d) -> " << result2.fromString("((x1)+d)") << "\n" << "raw: " << result2.data << "\n";
    std::cout << "((x1)d)+1 -> " << result2.fromString("((x1)d)+1") << "\n" << "raw: " << result2.data << "\n";
    std::cout << "\n";

    std::cout << "Polynomials with 3 variables:\n";
    std::vector<std::string> teststrings = { "3","x1","x2**2" };
    for ( auto it = teststrings.begin(); it != teststrings.end(); ++it ) {
        std::cout << "\"" << *it << "\" -> " <<  result.fromOneVarString( *it ) << "\n";
         //std::cout << "raw: " << result.data << "\n";
         std::cout << "\n";
    }

    std::vector<std::string> teststrings2 = {"(3+(1+7))","((x1))","((x1)+x2)", "x1**2 * ( d**2 - x1**2 - x2**2 + 2 x1 x2 )" };
    for ( auto it = teststrings2.begin(); it != teststrings2.end(); ++it ) {
        std::cout << "\"" << *it << "\" -> " <<  result.fromString( *it ) << "\n";
        std::cout << "raw: " << result.data << "\n";
        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Calculate some Powers:\n";
    Pol2Vars a0( "x,r", "x-r" );
    Pol2Vars unit( "x,r", "1" );
    int newnp = a0.data.getSize()[0];
    MathMatrix<MathMatrix<MathMatrix<double>>> a( newnp, 1 );
    a[0] = unit.data;
    std::cout << "Calculate (x-r)**2 manually: \n";
    std::cout << a0*a0 << "\n";
    std::cout << "Done!\n";
    for ( int i = 1; i < newnp; ++i ) {
        std::cout << " Calculate power " << i << ": ";
        a[i] = ( a0 * Pol2Vars( a0.varnames , a[i-1] ) ).data;
        std::cout << Pol2Vars( a0.varnames , a[i] ) << "\n";
    }

    std::cout << "Shape-Functions:\n";
    for ( int order = 1; order < 3; ++order ) {
        std::cout << "w_" << order << " = \n";
        for ( int j = 0; j < order; ++j ) {
            std::cout << " " << std::setw(2) << j << " < r/R < " << j+1 << ": "
                      << Shapes::Clouds::Weighting(order,j) << "\n";
        }
    }

    return 0;
}

