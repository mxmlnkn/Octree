/*

rm testMatrix.exe; g++ testMatrix.cpp -o testMatrix.exe -Wall -std=c++0x; ./testMatrix.exe
 
*/

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc

#include "Vector.h"
#include "BaseMatrix.h"


using namespace std;


int main( int argc, char **argv )
{
    int size[2] = {2,3};
    BaseMatrix<int,2> m( (Vec<int,2>)size /*Vec<int,2>(size) does not work! */ );
    Vec<int,2> ind(0);
    for (ind[1]=0; ind[1]<size[1]; ind[1]++ )
        for (ind[0]=0; ind[0]<size[0]; ind[0]++ ) {
            m[ind] = ind[1]*size[0]+ind[0];
        }
    m.Print();
    BaseMatrix<int,2> n(m);
    n.Print();
    
    cout << "Test indice range" << endl;
    
    /* Following will fail assertion. That all elements can be accessed is    *
     * provenby Print()                                                       */
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
       m.getVectorIndex(i).Print();
       cout << " -> " << m.getLinearIndex( m.getVectorIndex(i) ) << endl;
    }
    /* assertion Errors:
    m.getVectorIndex( m.getSize().product() );
    m.getVectorIndex( -1 );*/
    
    m.insertMatrix( Vec<int,2>(0), n );
    m.Print();
    
    size[1] = 1;
    BaseMatrix<int,2> p( (Vec<int,2>)size );
    p = 0;
    p.Print();
    int pos[2] = {0,1};
    m.insertMatrix( pos, p );
    m.Print();
    
    BaseMatrix<int,2> q;
    pos[0] = 1; pos[1] = 0;
    size[0]= 1; size[1]= 3;
    q = m.getPartialMatrix( pos, size );
    q.Print();
}
