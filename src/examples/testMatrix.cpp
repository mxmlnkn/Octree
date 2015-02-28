/*

rm testMatrix.exe; g++ testMatrix.cpp -o testMatrix.exe -Wall -std=c++0x; ./testMatrix.exe

*/

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc, srand, rand

#include "math/TVector.h"
#include "math/TBaseMatrix.h"
#include "yeesolver/YeeCell.h"

using namespace std;


int main( int argc, char **argv )
{
    cout << bool(0) << " " << bool(-1) << " " << bool(1) << "\n";

    int size3D[3] = {2,3,4};
    BaseMatrix<int,3> m3d( (Vec<int,3>)(size3D) );
    cout << m3d << "\n";

    cout << "\n To be interpolated:\n";
    int sizeraw[2] = {2,2};
    BaseMatrix<double,2> mraw( (Vec<int,2>)(sizeraw) );
    srand(2310748);
    for (int i=0; i<mraw.size.product(); i++)
        mraw[i] = double(rand()) / RAND_MAX;
    cout << mraw << "\n";
    int sizeint[2] = {30,30};
    BaseMatrix<double,2> mint( (Vec<int,2>)(sizeint) );
    mraw.NearestResizeTo( mint );
    cout << "\nAbove one interpolated:\n" << mint << "\n\n";

    int size[2] = {2,3};
    BaseMatrix<int,2> m( (Vec<int,2>)size /*Vec<int,2>(size) does not work! */ );
    Vec<int,2> ind(0);
    for (ind[1]=0; ind[1]<size[1]; ind[1]++ )
        for (ind[0]=0; ind[0]<size[0]; ind[0]++ ) {
            m[ind] = ind[1]*size[0]+ind[0];
        }
    cout << m << "\n";
    BaseMatrix<int,2> n(m);
    cout << n << "\n";

    cout << "Test indice range" << endl;

    /* Following will fail assertion. That all elements can be accessed is    *
     * proven by Print()     s                                                */
    /*{int ind[2] = {-1,-1};
    cout << m[ind];}
    {int ind[2] = {0,-1};
    cout << m[ind];}
    {int ind[2] = {3,0};
    cout << m[ind];}
    {int ind[2] = {0,4};
    cout << m[ind];}
    {int ind[2] = {3,4};
    cout << m[ind];}*/

    /* should especially print "->0" ... "->5" */
    for ( int i=0; i < m.getSize().product(); ++i ) {
       cout << m.getVectorIndex(i);
       cout << " -> " << m.getLinearIndex( m.getVectorIndex(i) ) << endl;
    }
    /* assertion Errors:
    m.getVectorIndex( m.getSize().product() );
    m.getVectorIndex( -1 );*/

    m.insertMatrix( Vec<int,2>(0), n );
    cout << m << "\n";

    size[1] = 1;
    BaseMatrix<int,2> p( (Vec<int,2>)size );
    p = 0;
    cout << p << "\n";
    int pos[2] = {0,1};
    m.insertMatrix( pos, p );
    cout << m << "\n";

    BaseMatrix<int,2> q;
    pos[0] = 1; pos[1] = 0;
    size[0]= 1; size[1]= 3;
    q = m.getPartialMatrix( pos, size );
    cout << q << "\n";

    size[0]=7;size[1]=7;
    BaseMatrix<int,2> r( (Vec<int,2>)size );
    for (int k=0; k<r.getSize().product(); ++k)
        r[k] = k;
    cout << r << "\n";

    Vec<int,2> newvec(size);
    cout << newvec << "\n";
    cout << (3+newvec) << "\n";
    cout << Vec<int,2>(1) << endl;
    cout << newvec << "\n";
}
