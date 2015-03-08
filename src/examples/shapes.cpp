#pragma once

#define DEBUG_SHAPES 99

namespace Shapes {

namespace Clouds {

    /********************************************
     *     _:_     _:_     :  : /:\ :  :        *
     *    | : | * | : |  = :  /  :  \  :        *
     *    | : |   | : |    :/ :  :  : \:        *
     * r:-1 0 1     x   r:-2 -1  0  1  2        *
     *  j=-1:j=0             j=-1:j=0  :j=2     *
     *           j==section         :j=1        *
     * section <= x <= section+1                *
     ********************************************/
    Pol1Var Weighting( int order, int section ) {
        Pol1Var result("x"); /* need order+1 powers ! */
        if ( order == 1 )
        {
            if ( section == -1 or section == 0 )
                result = "1";
            else
                result = "0";
        }
        else
        {
            Pol1Var tmp("x");
            Pol1Var intoldvar("x");
            Pol2Vars integrand("x,r");

            integrand = Weighting( order-1, section-1 ).renameVariables("x","r");
            result    = integrand.integrate( "r", "x-1", toString(section) );

            integrand = Weighting( order-1, section ).renameVariables("x","r");
            result   += integrand.integrate( "r", toString(section), toString(section+1) );

            integrand = Weighting( order-1, section+1 ).renameVariables("x","r");
            result   += integrand.integrate( "r", toString(section+1), "x+1" );
        }
        return result;
    }

}

Pol3Vars Int( std::string integrand, std::string var, std::string pa,
std::string pb, std::string paceil, std::string pbfloor, int shapeorder )
{
    /* Module[{k, a, b, aceil, bfloor, x, g},                 *
     * a = lims[[2]]; b = lims [[3]]; aceil = intlims[[1]];   *
     * bfloor = intlims[[2]];                                 *
     * g[x_] := Unevaluated[ReplaceAll[f, {lims[[1]] -> x}]]; */
    Pol3Vars g      ( "d|x1|x2", integrand );
    Pol3Vars a      ( "d|x1|x2", pa        );
    Pol3Vars b      ( "d|x1|x2", pb        );
    Pol3Vars aceil  ( "d|x1|x2", paceil     );
    Pol3Vars bfloor ( "d|x1|x2", pbfloor    );
    /* Unevaluated[ If[aceil - 1 == bfloor,                                   *
     *  (MI[np].PadRight[CoefficientList[x^2 g[x] wx[n, bfloor, x], x], np]). *
     *  (Vx[np, lims[[3]]] - Vx[np, lims[[2]]])                               *
     * ,                                                                      */
    if ( aceil - 1 == bfloor ) {
        std::cerr << "aceil-1 == bfloor\n";
    }
    /*    (MI[np].PadRight[
          CoefficientList[x^2 g[x] wx[n, aceil - 1, x], x], np]).
              (Vx[np, aceil] -
                Vx[np, lims[[2]]]) + (MI[np].PadRight[
                 CoefficientList[x^2 g[x] wx[n, bfloor, x], x], np]).
              (Vx[np, lims[[3]]] - Vx[np, bfloor]) +
             \!\(
        \*UnderoverscriptBox[\(\[Sum]\), \(k = aceil\), \(bfloor -
                1\)]\(\((MI[np] . PadRight[CoefficientList[
        \*SuperscriptBox[\(x\), \(2\)] g[x]\ wx[n, k, x], x],
                  np])\) . \[IndentingNewLine]\((Vx[np, k + 1] -
                 Vx[np, k])\)\)\)
            ]]];*/
    return g;
}

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

struct {
    int n;
    int j;

    /* fDEHI = Unevaluated[Int[1,{x2, 0, d - x1},{0, floor["d-x1"]},n,np] *
     *    + Int[FSSInt, {x2, d - x1, n}, {floor["d-x1"] + 1, n}, n, np]]; */
    Pol3Vars fDEHI( std::string floordmx1 )
    {
        std::string FSSInt = "0.25*x1*x2*(d**2 - (x1 - x2)**2)";
        return Int( "x1x2", "x2", "0", "d-x1", "0", floordmx1    , n )
             + Int( FSSInt, "x2", "0", "d-x1", "0", floordmx1+"1", n );
    }

    Pol3Vars operator() ( int pn, int pj ) {
        n = pn; j = pj;
        /* FSSInt = (d^2 - (x1 - x2)^2)/(4 x1 x2); i = Ceiling[j/2];          *
         * floor["d-n"] = n - i + 1; floor["d"] = 2 n - i;                    *
         * x2G["d-x1"] = k + d - floor["d"];                                  */
        const int i = j / 2 + ( j % 2 != 0 ) ? 1 : 0 ;
        const std::string floord   = toString( 2*n-1 );
        const std::string floordmn = toString( n-i+1 );

        return
        /* (* AB *) Int[1, {x1, 0, d - n}, {0, floor["d-n"] - 1}, n, np] 1/(  *
         *  4 \[Pi] \[Rho]0[n, "3D", 1]) +                                    */
        Int( "1", "x1", "0", "d-"+toString(n), "0", floordmn+"-1", n ) // * ( 1.0 / rho3D4pi(n) )
        /* +(* C(D+E+H+I) *) Int[ (floor["d-x1"] = n - 1; fDEHI),             *
         *   {x1, d - n, floor["d-n"]}, {floor["d-n"], floor["d-n"]}, n, np]  */
        + Int( fDEHI(toString(n)+"-1"), "x1", "d-"+toString(n), floordmn, floordmn,
            floordmn, n )
        /* +(* F(D+E+H+I)1 *) SplitInt[(floor["d-x1"] = floor["d"] - k; fDEHI)*
         *   , {x1, floor["d-n"], n}, {k, k, x2G["d-x1"]}, n, np]             */
        /* +(* F(D+E+H+I)2 *) SplitInt[(floor["d-x1"] = floor["d"] - k - 1;   *
         *  fDEHI) , {x1, n - i + 1, n}, {k, x2G["d-x1"], k + 1}, n, np]      */
        ;
    }
} f1;

#if 1==0
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

} // namespace shapes
