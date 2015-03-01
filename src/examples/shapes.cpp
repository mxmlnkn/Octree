/*
 rm shapes.exe; g++ shapes.cpp -o shapes.exe -Wall -Wextra -Wchar-subscripts -Wcomment -Wconversion -Wformat -Wmissing-braces -Wmissing-noreturn -Wparentheses -Wreturn-type -Wshadow -Wsign-compare -Wstrict-aliasing -Wuninitialized -Wunknown-pragmas -Wunreachable-code -std=c++0x -I ../libs/ -g 2>&1 | grep 'error'; ./shapes.exe
*/

#include <iostream>
#include <cstdlib>  // malloc

#include "math/TVector.h"
#include "math/TBaseMatrix.h"
#include "CompileTime.h"

using namespace std;

//typedef BaseMatrix<int,2> Matrix;

template<typename T_DTYPE>
class MathMatrix : public BaseMatrix<T_DTYPE,2> {
public:
    typedef Vec<int,2> VecI;

    MathMatrix( void );
    MathMatrix( VecI psize );
    MathMatrix( int m, int n );
    MathMatrix( const MathMatrix & m );
    ~MathMatrix( void );
    MathMatrix & operator=(const MathMatrix m); // if I make this a call by reference, then it won't be recognized as a specialization of the next assignment and thereby breaking everything, but copy-arguments could be memory intensive
    template<typename T_ETYPE> MathMatrix & operator=(T_ETYPE a); // broadcast value to all elements

    /* Conversion Operators */
    //template<int T_DIM> operator Vec<double,T_DIM>() const;
    //template<int T_DIM> Matrix& operator=(const Vec<double,T_DIM> val);

    MathMatrix & operator+=(const MathMatrix &mat);    // Add up to matrices
    MathMatrix & operator*=(const MathMatrix &mat);    // Matrix Multiplication
    MathMatrix & operator-=(const MathMatrix &mat);    // Subtraction
    inline MathMatrix & operator/=(double a);          // scalar Division based on Multiplication
    MathMatrix & operator*=(double a);                 // scalar Multiplication

    MathMatrix operator+(const MathMatrix &mat) const; // Add up to matrices
    MathMatrix operator*(const MathMatrix &mat) const; // Matrix Multiplication
    MathMatrix operator-(const MathMatrix &mat) const; // Subtraction
    inline MathMatrix operator/(double a) const;       // scalar Division based on Multiplication
    MathMatrix operator*(double a) const;              // scalar Multiplication
    bool operator==(const MathMatrix &mat) const;

    double (minor)(int row, int col) const;       // Counting from 1. clash with minor-macro -.-
    double det(void) const;
    MathMatrix invert(void) const;
    MathMatrix adjugate(void) const;
    MathMatrix transpose(void) const;
    /* Returns Norm of Matrix if it is a Vector, else error (returns -1) */
    double norm(void) const;
    /* Returns matrix with only positive elements (aij -> |aij|) */
    MathMatrix abs(void) const;
    int rank(void) const;
    double trace(void) const;
    MathMatrix rowEchelon(void) const;
    inline bool isSquare(void) const;
    inline bool isVector(void) const;

    T_DTYPE & operator()( int i, int j );     // Get reference to element
    T_DTYPE   operator()( int i, int j ) const;     // Just read element
    VecI getSize(void) const;
    void setSize(int m, int n) const;

    int getVectorDim(void) const;
    int getSquareDim(void) const;

    void setDiagonal(T_DTYPE a=1, int k = 0);
    void setDiagonal(T_DTYPE a[], int k = 0);
    template<typename T_ETYPE> void setAll(T_DTYPE a);

    void delRow(int row, int amount=1);         // Deletes i-th Row counting from 1
    void delCol(int col, int amount=1);
    MathMatrix augment(const MathMatrix &mat) const;    //conjoines object with mat (if number of rows is equivalent)

    //explicit operator Vec();
};

/******************************** Constructors ********************************/

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( void )
 : BaseMatrix<T_DTYPE,2>()
{ };

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( int m, int n )
 : BaseMatrix<T_DTYPE,2>( VecI(m,n) )
{ }

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( VecI psize )
 : BaseMatrix<T_DTYPE,2>( VecI( psize[0], psize[1] ) )
{ }

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( const MathMatrix & m )
 : BaseMatrix<T_DTYPE,2>( m )
{ }

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::~MathMatrix( void )
/* no explicit call for parent class destructor needed */
{ }

/**************************** Assignment operators ****************************/
// broadcast value to all elements
template<typename T_DTYPE>
template<typename T_ETYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator= ( T_ETYPE a )
{
    BaseMatrix<T_DTYPE,2>::operator=(a);
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator= (const MathMatrix m)
{
    /* Can't call copy assignment from base-class, because it expects a
     * BaseMatrix argument instead of MathMatrix. We would have to convert it
     * first */
    BaseMatrix<T_DTYPE,2>::operator=( (BaseMatrix<T_DTYPE,2>) m );
    return *this;
}

template<typename T_DTYPE>
void MathMatrix<T_DTYPE>::setSize( int m, int n ) const {
    VecI psize(m,n);
    MathMatrix mnew(psize);

    /* copy values if possible into new matrix */
    VecI tocopy = psize.min( this->size );
    for ( int i = 0; i < mnew.size.product(); ++i ) {
        VecI vind = ConvertLinearToVectorIndex<2>(i, mnew.size);
        if ( vind < this->size )
            mnew[vind] = (*this)[vind];
    }

    (*this) = mnew;
}

/************************ Binary Arithmetic Operators *************************/

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator+=(const MathMatrix &mat)
{
    assert( this->size == mat.size );
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] += mat[i];
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator-=(const MathMatrix &mat)
{
    assert( this->size == mat.size );
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] -= mat[i];
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator*=(const MathMatrix &mat)
{
    *this = (*this) * mat;
    return *this;
}

template<typename T_DTYPE>
inline MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator/=(double a)
{
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] /= a;
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator*=(double a)
{
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] *= a;
    return *this;
}

/**************************** Arithmetic Operators ****************************/

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator+(const MathMatrix &mat) const
{
    MathMatrix<T_DTYPE> res( *this );
    res += mat;
    return res;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator-(const MathMatrix &mat) const
{
    MathMatrix<T_DTYPE> res( *this );
    res += mat;
    return res;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator*(const MathMatrix &mat) const
{
    assert( this->size[1] == mat.size[0] );
    MathMatrix res( VecI( this->size[0], mat.size[1] ) );
    for ( int j = 0; j < res.size[1]; j++ ) {       // go down every column
        for ( int i = 0; i < res.size[0]; i++ ) {   // go down every line
            T_DTYPE sum = (*this)(i,j); // copy (*this)(i,j), to to copy the dimension, this also works, if T_DTYPE is a double!
            sum = 0;
            if ( i==0 and j==0 ) {
                std::cout << "sum:" << sum << "\n";
            }
            for ( int k = 0; k < this->size[1]; k++ ) // dot product
	           sum += (*this)[VecI(i,k)] * mat[VecI(k,j)];
            res[VecI(i,j)] = sum;
        }
    }
    return res;
}

template<typename T_DTYPE>
inline MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator/(double divisor) const
{
    return (*this) * ( 1.0 / divisor );
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator*(double a) const
{
    MathMatrix<T_DTYPE> res( this->size );
    for ( int i = 0; i < this->size.product(); ++i )
        res[i] = (*this)[i] * a;
    return res;
}

template<typename T_DTYPE>
bool MathMatrix<T_DTYPE>::operator==(const MathMatrix &mat) const
{
    if ( this->size != mat.size ) {
        return false;
    } else {
        for ( int i = 0; i < this.size.product(); ++i )
            if ( (*this)[i] != mat[i] )
                return false;
    }
    return true;
}

/******************************** Test methods ********************************/

template<typename T_DTYPE>
inline bool MathMatrix<T_DTYPE>::isSquare( void ) const
{
    return this->size[0] == this->size[1];
}

template<typename T_DTYPE>
inline bool MathMatrix<T_DTYPE>::isVector( void ) const
{
    return ( this->size[0] == 1 || this->size[1] == 1 );
}

/****************************** Modifying methods *****************************/

/* k=0 is main diagonal, k>0 is above main diagonal, k<0 is below */
template<typename T_DTYPE>
void MathMatrix<T_DTYPE>::setDiagonal( T_DTYPE a, int k )
{
    if ( k >= 0 ) {
        assert( k < this->size[0] );
        for ( int i = 0; i < this->size[1] - k; ++i )
            (*this)( i,i+k ) = a;
    } else {
        assert( -k < this->size[1] );
        for ( int i = 0; i < this->size[0] + k; ++i )
            (*this)( i-k,i ) = a;
    }
    return;
}

template<typename T_DTYPE>
void MathMatrix<T_DTYPE>::setDiagonal( T_DTYPE a[], int k )
{
    if ( k >= 0 ) {
        assert( k < this->size[0] );
        for ( int i = 0; i < this->size[1] - k; ++i )
            (*this)( i,i+k ) = a[i];
    } else {
        assert( -k < this->size[1] );
        for ( int i = 0; i < this->size[0] + k; ++i )
            (*this)( i-k,i ) = a[i];
    }
    return;
}

/******************************* Reading methods*******************************/
template<typename T_DTYPE>
typename MathMatrix<T_DTYPE>::VecI MathMatrix<T_DTYPE>::getSize(void) const
{
    return this->size;
}

template<typename T_DTYPE>
T_DTYPE & MathMatrix<T_DTYPE>::operator()( int i, int j )
{
    return (*this)[ VecI(i,j) ];
}

template<typename T_DTYPE>
T_DTYPE MathMatrix<T_DTYPE>::operator()( int i, int j ) const
{
    return (*this)[ VecI(i,j) ];
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

int main( void ) {
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

/*    MathMatrix<double> spacedconst(np,1); spacedconst = 0;
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

/****************************** Expand (d-2x)**2 ******************************/
    MathMatrix<int> dspace(np,1);
    dspace = 0;
    MathMatrix<MathMatrix<int> > xdspace(np,1);
    xdspace = dspace;
     dspace[2] =  1;
    xdspace[1] =  dspace;
     dspace[2] =  0;
     dspace[3] = -2;
    xdspace[2] =  dspace;

    const char spacevarnames[3][10] = { "x1", "d" };
    const int nestedlevels = 2;
    Vec<int,nestedlevels> size;
    size[0] = xdspace.getSize()[0];
    size[1] = xdspace[0].getSize()[0];

    MathMatrix<int> mzero(np,np); mzero = 0;
    MathMatrix<int> munit(np,np); munit.setDiagonal(1);
    MathMatrix<int> mmul (np,np); mmul.setDiagonal(1,-1);
    std::cout << "mzero: " << mzero << "\n";
    std::cout << "munit: " << munit << "\n";
    std::cout << "mmul : " << mmul  << "\n";

    MathMatrix<MathMatrix<int> > mx(np,np);
    mx = mzero;
    mx.setDiagonal(munit,-1);
    MathMatrix<MathMatrix<int> > md(np,np);
    md = mzero;
    md.setDiagonal(mmul,0);
    std::cout << "mx: " << mx << " size: " << mx.getSize() << "\n";
    std::cout << "md: " << md << " size: " << md.getSize() << "\n";

    std::cout << "md+mx: " << (md+mx) << "\n";
    std::cout << "md+mx: " << (md+mx) << "\n";


    for ( int i = 0; i < nestedlevels; ++i )
        std::cout << " " << spacevarnames[i] << "_{";

    BaseMatrix<int,1> index(nestedlevels);
    index = 0;
    for ( int i = 0; i < size.product(); ++i )
    {
        /* recursively unpack xdspace which is Mathmatrix<Mathmatrix<...>...> *
         * using a local function i.e. struct lambda ...                      */

        std::cout << xdspace[index[0]][index[1]];

        /* increment vector index */
        int closedBrackets = 0;
        index[nestedlevels-1] += 1;
        for ( int j = nestedlevels-1; j >= 0; --j ) {
            if ( index[j] >= size[j] ) {
                if ( j > 0 ) {
                    index[j] = 0;
                    index[j-1] += 1;
                } else
                    assert( i == size.product()-1 );
                closedBrackets++;
                std::cout << "}_" << spacevarnames[nestedlevels - closedBrackets]
                          << " ";
            }
        }
        if ( i != size.product()-1 ) {
            std::cout << ",";
            for ( int j = 0; j < closedBrackets; ++j )
                std::cout << " " << spacevarnames[nestedlevels - closedBrackets]
                          << "_{";
        }
    }
    /**************************************************************************/

    return 0;
}
