/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe
make testOctree && ./testOctree.exe

*/

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc
#include <random>   // normal_distribution
#include "CompileTime.h"
#include "math/TVector.h"
#include "octree/Octree.h"
#include "octree/OctreeToSvg.h"
#include "teestream/TeeStream.h"

#define SIMDIM 2
typedef Vec<double,SIMDIM> VecD;


int main( int argc, char **argv )
{
    tout.Open("out");


    const int NValues = 100;
    int * data = (int*) malloc( sizeof(int)*NValues );
    for (int i=0; i<NValues; ++i)
        data[i] = i;

    VecD size(2), center(0);
    Octree::Octree<int,SIMDIM> tree( center, size );

    std::default_random_engine rng( 7542168 );
    std::normal_distribution<double> x_distribution(  0.2, 0.5 );
    std::normal_distribution<double> y_distribution( -0.3, 0.2 );

    for (int i=0; i<NValues; i++) {
        VecD position;
        position[0] = x_distribution(rng);
        if ( position[0] > +0.5*size[0] ) position[0] = +0.5*size[0] - 1e-6;
        if ( position[0] < -0.5*size[0] ) position[0] = -0.5*size[0];
        position[1] = y_distribution(rng);
        if ( position[1] > +0.5*size[1] ) position[1] = +0.5*size[1] - 1e-6;
        if ( position[1] < -0.5*size[1] ) position[1] = -0.5*size[1];
        tree.InsertData( position, &(data[i]) );
    }

    tout << tree;
    Octree::OctreeToSvg<int,SIMDIM> svgoutput( tree, std::string("octree") );
    svgoutput.PrintGrid();
    svgoutput.PrintPositions();
    /*//std::cout << "\n====== Compare Tree with itself =====" << std::endl;
    //svgoutput.AnimateUpdated( tree ); // works :)
    Octree::Octree<int,SIMDIM> tree2 = tree;
    std::cout << "======== Tree 2 ========" << std::endl;
    //std::cout << tree2;
    const int PID = 17;
    assert( PID < NValues );
    std::cout << "==== Find Data[17] ====" << std::endl;
    std::cout << tree2.FindData( &(data[PID]) );
    VecD diffvec(0);
    diffvec[1] = -0.1;
    VecD pos = tree2.FindData( &(data[PID]) );
    std::cout << "==== Moving Data[17] ====" << std::endl;
    tree2.MoveData( pos, &(data[PID]), pos+diffvec );
    std::cout << tree2;
    svgoutput.AnimateUpdated( tree2 ); */


    std::normal_distribution<double> dist( 0, 0.5 );
    for (int t=0; t<20; t++) {
        for (int i=0; i<NValues; i++) {
            VecD pos = tree.FindData( &(data[i]) );
            VecD diffvec(0);
            diffvec[0] = dist(rng);
            diffvec[1] = dist(rng);
            VecD newpos = pos+diffvec;
            if ( newpos[0] > +0.5*size[0] )
                newpos[0] = fmod( ( 0.5*size[0] + newpos[0] ) ,
                                  ( size[0]-1e-6 ) ) - 0.5*size[0];
            if ( newpos[0] < -0.5*size[0] )
                newpos[0] = fmod( ( 0.5*size[0] + -newpos[0] ) ,
                                  ( size[0]-1e-6 ) ) - 0.5*size[0];
            if ( newpos[1] > +0.5*size[1] ) newpos[1] =
                fmod( ( 0.5*size[1] + newpos[1] ) ,
                      ( size[1]-1e-6 ) ) - 0.5*size[1];
            if ( newpos[1] < -0.5*size[1] ) newpos[1] =
                fmod( ( 0.5*size[1] + -newpos[1] ) ,
                      ( size[1]-1e-6 ) ) - 0.5*size[1];

            tree.MoveData( pos, &(data[i]), newpos );
            if ( !tree.CheckIntegrity() ) {
                std::cout << "DIE HOELLE IST LOS!!!\n";
                return 1;
            }
        }
        //std::cout << "==t:" << t << "==" << std::endl;
        //tout << tree;
        svgoutput.AnimateUpdated( tree );
    }

    free(data);
}