/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe
make testOctree && ./testOctree.exe

*/

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc, rand
#include <random>   // normal_distribution
#include <ctime>    // time

namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime

#include "paramset/Parameters_2015-01-16.cpp"
#include "math/TVector.h"
#include "octree/Octree.h"
#include "octree/OctreeToSvg.h"
#include "teestream/TeeStream.h"

#define SIMDIM 2
typedef Vec<double,SIMDIM> VecD;
typedef Vec<double,SIMDIM> VecI;

#define RANDOM_SEED 24756139



int main( int argc, char **argv )
{
    tout.Open("out");
    srand( RANDOM_SEED );

/* Run this Programm for several Ordering Methods ! */
for ( int ORDERING = 0; ORDERING < 3; ++ORDERING ) {
/* Run this Programm for several world sizes ! */
    std::ofstream resultsFile;
    std::stringstream sOrdering;
    switch (ORDERING) {
        case 0: sOrdering << "Morton_"  ; break;
        case 1: sOrdering << "GrayCode_"; break;
        case 2: sOrdering << "Hilbert_" ; break;
    } sOrdering << "Ordering_";
    resultsFile.open( std::string("Octree_Benchmark_point_minRecursion4_") + sOrdering.str() + std::string(".dat") );
    resultsFile << "# worldsize totalTraffic interTraffic messageCount\n" << std::flush;

for ( int worldsize = 1; worldsize < 2; ++worldsize ) {
    tout << "World Size      : " << worldsize << "\n";

    VecD size(100), center(50), globSize(size), globCenter(center);
    typedef Octree::Octree<SIMDIM> OctreeType;
    Octree::Octree<SIMDIM> tree( center, size );
    Octree::Octree<SIMDIM>::iterator it;

    #define INITSETUP 6
    #if INITSETUP == 1
    /* refine all cells to some value */
        tout << "Initial homogenous Refinement\n";
        int maxRecursion = 2;
        for ( int i=0; i<maxRecursion; i++) {
            tout << i << "th Refinement\n";
            it = tree.begin();
            while ( it!=tree.end() ) {
                if ( it->IsLeaf() )
                    (it++)->GrowUp();
                else
                    ++it;
            }
        }
    #endif
    #if INITSETUP == 2
        /* refine a subsection of the cells one more time */
        tout << "Refinement of a subsection\n";
        it = tree.begin();
        while ( it!=tree.end() ) {
            if ( it->IsLeaf() and it->center.norm() < 0.4 )
                (it++)->GrowUp();
            else
                ++it;
        }
        it = tree.begin();
        while ( it!=tree.end() ) {
            if ( it->IsLeaf() and it->center.norm() < 0.25 )
                (it++)->GrowUp();
            else
                ++it;
        }
    #endif
    #if INITSETUP == 3
        /* refine a subsection of the cells one more time */
        tout << "Refinement of a subsection\n";
        it = tree.begin();
        while ( it!=tree.end() ) {
            if ( it->IsLeaf() and it->center[0] > 0 )
                (it++)->GrowUp();
            else
                ++it;
        }
        VecD pos(0); pos[0]=0.375; pos[1]=-0.125;
        tree.FindLeafContainingPos( pos*tree.size )->GrowUp(); */
    #endif
    #if INITSETUP == 4
        /* refine cells including fixed point many times */
        int maxRecursion = 10;
        VecD pos(0); pos[0]=0.1; pos[1]=0.1;
        for ( int i = 0; i < maxRecursion; ++i )
            tree.FindLeafContainingPos( pos*tree.size )->GrowUp();
    #endif
    #if INITSETUP == 5
        int targetCells  = 200;
        int currentCells = 1;
        while ( currentCells < targetCells ) {
            VecD pos(0);
            pos[0] = double(rand()) / double(RAND_MAX) - 0.5;
            pos[1] = double(rand()) / double(RAND_MAX) - 0.5;
            tree.FindLeafContainingPos( pos*tree.size )->GrowUp();
            currentCells += 3;
        }
    #endif
    #if INITSETUP == 6
    /**************************************************************************/
    /* (1) Setup Octree Refinement ********************************************/
    /**************************************************************************/

    /********* refine all cells to initial homogenous min-Refinement **********/
    for ( int lvl=0; lvl<INITIAL_OCTREE_REFINEMENT; lvl++) {
        for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it )
            if ( it->IsLeaf() and it->getLevel()==lvl ) it->GrowUp();
    }
    /*********************** Refine certain boundaries ************************/
    assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
    VecD M(0.5*globSize);          // center of circle
    double R = 0.2*globSize.min(); // radius of circle
    for ( int lvl=INITIAL_OCTREE_REFINEMENT; lvl<MAX_OCTREE_REFINEMENT; lvl++) {
        /* Get all circle angles, where it intersects with a cell border */
        std::list<double> lphi;
        std::list<double>::iterator it;
        VecD cellsize = globSize / pow(2,lvl);
        double linexmin = ceil ( (M[0]-R)/cellsize[0] ) * cellsize[0];
        double linexmax = floor( (M[0]+R)/cellsize[0] ) * cellsize[0];
        for ( double linex=linexmin; linex<=linexmax; linex += cellsize[0] ) {
            /* acos in [0,2*pi] */
            double phi = acos( (linex-M[0])/R );
            lphi.push_back( phi );
            /* also add value mirrored at y-axis to stack */
            lphi.push_back( 2*M_PI-phi );
        }
        double lineymin = ceil ( (M[1]-R)/cellsize[1] ) * cellsize[1];
        double lineymax = floor( (M[1]+R)/cellsize[1] ) * cellsize[1];
        for ( double liney=lineymin; liney<=lineymax; liney += cellsize[1] ) {
            /* asin in [-pi,pi] */
            double phi = asin( (liney-M[1])/R );
            lphi.push_back( phi < 0 ? 2*M_PI+phi : phi );
            /* also add value mirrored at x-axis to stack */
            lphi.push_back( M_PI-phi );
        }
        lphi.sort();

        /* Echo all found angles */
        #if DEBUG_MAIN_YEE >= 100
            tout << "Angle list contains:";
            for (it=lphi.begin(); it!=lphi.end(); ++it)
                tout << ' ' << *it;
            tout << '\n';
        #endif

        /* Grow up all cells, with which the circle intersects. Find them by  *
         * using an angle between to successive circle intersection angles    */
        for (it=lphi.begin(); it!=lphi.end(); ++it) {
            VecD pos(0);
            std::list<double>::iterator itnext = it;
            double phi;
            if ( ++itnext == lphi.end() ) {
                itnext = lphi.begin();
                phi = 0.5 * (2*M_PI + *itnext + *it);
            } else
                phi = 0.5 * (*itnext + *it);
            pos[0] = M[0] + R*cos(phi);
            pos[1] = M[1] + R*sin(phi);
            OctreeType::Node * node = tree.FindLeafContainingPos(pos);
            if ( node->getLevel() == lvl )
                node->GrowUp();
        }
    }
    #endif
    std::cout << "Tree-Integrity: " << tree.CheckIntegrity() << "\n";

    for ( OctreeType::iterator it=tree.begin(1); it!=tree.end(); ++it)
    if ( it->IsLeaf() )
        tout << it->center << "\n";
    return 1;

    /* Count Cells */
    tout << "Count Cells...";
    int NValues = 0;
    it = tree.begin();
    while ( it!=tree.end() ) {
        if ( it->IsLeaf() )
            NValues++;
        ++it;
    }
    tout << NValues << std::endl;
    tout << "Count Cells internally..." << tree.root->countLeaves() << std::endl;

    if ( worldsize > NValues )
        break;

    /* allocate data (which stores assigned ranks) to which pointers will     *
     * given to octree. And default it to the last rank                       */
    int * data = (int*) malloc( sizeof(int)*NValues );
    for (int i=0; i<NValues; ++i)
        data[i] = worldsize-1;

    /* Insert testDate (later YeeCell-Data or Absorbercelldata, or Guard) at  *
     * the center of every leaf node. By default all cells will be assigned   *
     * to last rank                                                           */
    int dataInserted = 0;
    for ( it = tree.begin( ORDERING ); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) {
            assert( dataInserted < NValues );
            assert( it->data.empty() );
            it->data.push_back( &(data[dataInserted]) );
            ++dataInserted;
        }

    /* Calculate total costs of all cells. Could be done when counting cells  */
    double totalCosts = 0;
    it = tree.begin();
    while ( it!=tree.end() ) {
        if ( it->IsLeaf() )
            totalCosts += 1. / it->size.min();
        ++it;
    }
    double optimalCosts = totalCosts / double(worldsize);
    tout << "Total Costs: " << totalCosts << " => Optimal Costs: " << optimalCosts << std::endl;

    /* Print tree to SVG */
    std::stringstream sWorldsize; sWorldsize << worldsize;
    std::stringstream sOrdering;
    switch (ORDERING) {
        case 0: sOrdering << "Morton_"  ; break;
        case 1: sOrdering << "GrayCode_"; break;
        case 2: sOrdering << "Hilbert_" ; break;
    } sOrdering << "Ordering_";
    Octree::OctreeToSvg<SIMDIM> svgoutput( tree, std::string("Octree_") + sOrdering.str() + std::string("worldsize_") + sWorldsize.str() );
    svgoutput.PrintGrid();

    /* Assign cells to all the processes */
    double cumulativeCosts = 0;
    int curRank = 0;
    int currentTime  = 0;
    double delay     = 1./64.; // 8 frames per second at max achievable with SVG
    double rankDelay = 1./64.;

    Octree::Octree<SIMDIM>::iterator it0 = tree.begin();
    Octree::Octree<SIMDIM>::iterator it1 = tree.begin();
    for ( it = tree.begin( ORDERING ); it!=tree.end(); ++it ) {
        if ( it->IsLeaf() ) {
            cumulativeCosts += 1. / it->getSize().min();
            if ( cumulativeCosts >= optimalCosts and curRank != worldsize-1) {
                cumulativeCosts = 1. / it->getSize().min();
                curRank++;
            }
            *((int*)it->data[0]) = curRank;

            /* Graphical output of traversal line */
            it0 = it1;
            it1 = it;
            if ( it0->IsLeaf() and it1->IsLeaf() ) {
                size_t id = reinterpret_cast<size_t>( it0->data[0] );
                VecD r0 = svgoutput.convertToImageCoordinates( it0->center );
                VecD r1 = svgoutput.convertToImageCoordinates( it1->center );
                /* Spawn invisible line element */
                svgoutput.out
                  << "<line"                                               "\n"
                  << " id=\"path" << id << "\""                            "\n"
                  << " x1=\"" << r0[0] << "px\" y1=\"" << r0[1] << "px\""  "\n"
                  << " x2=\"" << r1[0] << "px\" y2=\"" << r1[1] << "px\""  "\n"
                  << " style=\"stroke:none;stroke-width:3px\""             "\n"
                  << "/>"                                                  "\n";
                /* Animate line element to become visible after some time */
                svgoutput.out
                  << "<set"                                            "\n"
                  << " xlink:href=\"#path" << id << "\""               "\n"
                  << " attributeName=\"stroke\""                       "\n"
                  << " begin=\"" << delay*double(currentTime) << "s\"" "\n"
                  << " to   =\"#008000\""                              "\n"
                  << "/>"                                              "\n";
                currentTime += rankDelay / delay;
            }

            /* Graphical output of and cell-rank-mapping */
            svgoutput.boxesDrawnIt = svgoutput.boxesDrawn.find( it->center );
            int id = svgoutput.boxesDrawnIt->second.id;
            int r  = int( 128 * double(curRank+1) / double(worldsize) );
            int g  = r;
            int b  = r;
            /*int r  = 255;;
            int g  = 0;
            int b  = 255 * double(curRank) / double(worldsize); */
            svgoutput.out
                << "<set"                                                "\n"
                << " xlink:href=\"#" << id << "\""                       "\n"
                << " attributeName=\"fill\""                             "\n"
                << " begin=\"" << delay*double(currentTime) << "s\""     "\n"
                << " to   =\"rgb(" << r << "," << g << "," << b << ")\"" "\n"
                << "/>"                                                  "\n";
        }
    }

    /* Count Neighbors intra- and interprocessdata to transmit */
    double * costs = new double[worldsize];
    for (int i=0; i<worldsize; ++i)
        costs[i] = 0;

    double interTraffic = 0;
    double totalTraffic = 0;
    const int bytesPerCell = (6+4)*8;

    it = tree.begin( ORDERING );
    while ( it!=tree.end() ) {
        if ( it->IsLeaf() ) {
            costs[ *((int*)it->data[0]) ] += 1. / it->getSize().min();
            int nNeighbors = 0;
            int nLeavesOnOtherNodes = 0;
            int thisRank = *((int*)it->data[0]);

            VecI dir[4];
            dir[0][0]=+1; dir[0][1]= 0;
            dir[1][0]=-1; dir[1][1]= 0;
            dir[2][0]= 0; dir[2][1]=+1;
            dir[3][0]= 0; dir[3][1]=-1;

            for ( int lindir = 0; lindir < 4; lindir++ ) {
                Octree::Node<SIMDIM> * neighbor = it->getNeighbor( dir[lindir] );
                if ( neighbor == NULL )
                    continue;

                nNeighbors += neighbor->countLeaves();

                Octree::Octree<SIMDIM>::iterator itn = neighbor->begin();
                while ( itn != neighbor->end() ) {
                    if ( itn->IsLeaf() ) {
                        int neighborRank = *((int*)itn->data[0]);
                        if ( thisRank != neighborRank )
                            nLeavesOnOtherNodes++;
                    }
                    ++itn;
                }
            }

            totalTraffic += bytesPerCell * nNeighbors;
            interTraffic += bytesPerCell * nLeavesOnOtherNodes;

            //tout << it->center << " needs data from "
            //     << nNeighbors << " neighbors. " << nLeavesOnOtherNodes << " of those are not on this process\n";
        }
        ++it;
    }

    tout << "Number of Cells : " << tree.root->countLeaves() << "\n";
    for (int i=0; i<worldsize; i++) {
        tout << "Cost assigned to rank " << i << " is " << costs[i] << "\n";
    }
    tout << "Total data to be read from neighbors : " << totalTraffic << "\n";
    tout << "Data which needs to be communicated  : " << interTraffic << "\n";

    resultsFile << worldsize << " " << totalTraffic << " " << interTraffic << "\n" << std::flush;

    delete[] costs;
    free(data);

} /* Run this Programm for several world sizes ! */
    resultsFile.close();
} /* Run this Programm for several Ordering Methods ! */
}
