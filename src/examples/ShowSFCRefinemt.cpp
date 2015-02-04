/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe
make testOctree && ./testOctree.exe

*/

#include <iostream>
#include <cstdlib>  // malloc, rand
#include <ctime>    // time

#include "octree/Octree.h"
#include "octree/OctreeToSvg.h"
#include "Colors.h"
#include "Directions.h"

const int SIMDIM = 2;

int main( int argc, char **argv )
{
    /* Create timestamp and basefolder for filenames */
    time_t tmp_time = time(0);
    struct tm * now = localtime( &tmp_time );
    std::stringstream timestamp;
    timestamp  << 1900 + now->tm_year << "-" << std::setfill('0')
               << std::setw(2) << 1 + now->tm_mon << "-"
               << std::setw(2) << now->tm_mday << "_"
               << std::setw(2) << now->tm_hour << "-"
               << std::setw(2) << now->tm_min  << "-"
               << std::setw(2) << now->tm_sec;

    /* Run this Programm for several Ordering Methods ! */
    for ( int ORDERING = 0; ORDERING <= 2; ++ORDERING ) {
        std::ofstream resultsFile;
        std::stringstream sOrdering;
        switch (ORDERING) {
            case 0: sOrdering << "Morton"  ; break;
            case 1: sOrdering << "GrayCode"; break;
            case 2: sOrdering << "Hilbert" ; break;
            case 3: sOrdering << "Rows" ; break;
            case 4: sOrdering << "Four-Color-Theorem" ; break;
        }

        Octree::Octree<SIMDIM> tree( 0, 1 );
        Octree::OctreeToSvg<SIMDIM> svgoutput( tree, timestamp.str() + 
                            sOrdering.str() + std::string(".svg"), false );

        int MAX_OCTREE_REFINEMENT = 5;
        for ( int i=0; i < MAX_OCTREE_REFINEMENT; i++)
            for ( Octree::Octree<SIMDIM>::iterator it = tree.begin();
            it != it.end(); )
                if ( it->IsLeaf() and it->getLevel() == i )
                    (it++)->GrowUp();
                else
                    ++it;
        //svgoutput.PrintGrid();
        tree.root->DeleteChildren();
    
        for( int curlvl = 1; curlvl <= MAX_OCTREE_REFINEMENT; curlvl++ ) {
            /* refine all cells to some value */
            for ( int i=0; i < curlvl; i++)
                for ( Octree::Octree<SIMDIM>::iterator it = tree.begin();
                it != it.end(); )
                    if ( it->IsLeaf() and it->getLevel() == i )
                        (it++)->GrowUp();
                    else
                        ++it;

            struct ColorFunc {
                int color;
                ColorFunc(int pcolor) : color(pcolor) {}
                Vec<int,3> operator() ( double x ) {
                    //int index = int( floor( x * (double) Colors::Own1Length ) );
                    int index = color % int(Colors::Own1Length);
                    return Colors::getColorVector( Colors::Own1[index] );
                }
            } colorfunc(curlvl-1);
            svgoutput.PrintTraversal(ORDERING,colorfunc,1./pow(2,curlvl+1));
        }
        svgoutput.close();

    } /* Run this Programm for several Ordering Methods ! */
}
