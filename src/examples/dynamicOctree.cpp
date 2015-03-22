/*
rm dynamicOctree.exe; g++ dynamicOctree.cpp -o dynamicOctree.exe -I../libs -I.. -Wall -std=c++0x; ./dynamicOctree.exe
*/

#include <iostream>
#include <cmath>    // sin
#include "octree/Octree.h"
#include "octree/OctreeToSvg.h"

#define SIMDIM 2
typedef Vec<double,SIMDIM> VecD;


int main( int argc, char **argv )
{
    VecD size(1), center(0);
    Octree::Octree<SIMDIM> tree( center, size );
    Octree::OctreeToSvg<SIMDIM> svgoutput( tree, std::string("dynamicOctree.svg") );

    const double T_MAX = 10.0;
    const int N_FRAMES = 200;
    for ( double t = 0; t < T_MAX; t += T_MAX / double(N_FRAMES) )
    {
        std::cout << "t = " << t << " tree integrity:" << tree.CheckIntegrity() << "\n";
        /* refine all cells to a homogenous level, rejuvenating those with finer level */
        const int HOM_LEVEL = 4;
        bool didchange = true;

        while ( didchange ) {
            didchange = false;
            for ( auto it = tree.begin(); it != it.end(); /* do nothing */ )
                if ( it->IsLeaf() )
                {
                    if ( it->getLevel() < HOM_LEVEL ) {
                        (it++)->GrowUp();
                        didchange = true;
                    } else
                        ++it;
                } else if ( it->getLevel() >= HOM_LEVEL ) {
                    it->Rejuvenate();
                    didchange = true;
                    break;
                    /* Rejuvenating thrashes iterator, so break iterator loop */
                } else
                    ++it;
        }

        const VecD  center = VecD( fmod(0.5 + t/(0.5*T_MAX), 1.0) - 0.5, 0.0 );
        const float radius = 0.15 * tree.size[0] * ( 1.8 - cos( 2.0*3.141592653*t/(0.5*T_MAX) ) );
        std::cout << "center:" << center << " radius:" << radius << "\n";

        /* refine a subsection ( 2D-Ball ) of the cells one more time */
        for ( auto it = tree.begin(); it != it.end(); /* do nothing */ )
            if ( it->IsLeaf() )
            {
                if ( (it->center-center).norm() < radius or
                     (it->center-center - VecD(1.0,0.0)).norm() < radius or
                     (it->center-center + VecD(1.0,0.0)).norm() < radius
                )
                    (it++)->GrowUp();
                else
                    ++it;
            } else
                ++it;

        if ( t == 0 ) {
            svgoutput.PrintGrid();
        } else
            svgoutput.AnimateUpdated( tree, t );
    }
    svgoutput.close();

    return 0;
}