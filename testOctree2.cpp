/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe
make testOctree && ./testOctree.exe

*/

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc
#include <random>   // normal_distribution

namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime

#include "Vector.h"
#include "Octree.h"
#include "OctreeToSvg.h"
#include "TeeStream.h"

#include "OctreeToSvg.tpp"

#define SIMDIM 2
typedef Vec<double,SIMDIM> VecD;
typedef Vec<double,SIMDIM> VecI;

int main( int argc, char **argv )
{
    tout.Open("out");
    
    //int ORDERING = Octree::Ordering::Morton;
    
    
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
    resultsFile.open( std::string("Octree_Benchmark_circle_minRecursion4_") + sOrdering.str() + std::string(".dat") );
    resultsFile << "# worldsize totalTraffic interTraffic\n";
    
for ( int worldsize = 1; worldsize < 128; ++worldsize ) {
    
    VecD size(100), center(0);
    Octree::Octree<int,SIMDIM> tree( center, size );

    /* refine all cells to some value */
    tout << "Initial homogenous Refinement\n";
	Octree::Octree<int,SIMDIM>::iterator it;
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

    /* refine a subsection of the cells one more time */
    /*tout << "Refinement of a subsection\n";
	it = tree.begin();
	while ( it!=tree.end() ) {
		if ( it->IsLeaf() and it->center[0] > 0 )
			(it++)->GrowUp();
		else
			++it;
	}
    VecD pos(0); pos[0]=0.375; pos[1]=-0.125;
    tree.FindLeafContainingPos( pos*tree.size )->GrowUp(); */

    /* Count Cells */
    tout << "Count Cells\n";
    int NValues = 0;
	it = tree.begin();
	while ( it!=tree.end() ) {
		if ( it->IsLeaf() )
			NValues++;
        ++it;
	}

    /* allocate data (which stores assigned ranks) to which pointers will     *
     * given to octree. And default it to the last rank                       */
    int * data = (int*) malloc( sizeof(int)*NValues );
    for (int i=0; i<NValues; ++i)
        data[i] = worldsize-1;

    /* Insert testDate (later YeeCell-Data or Absorbercelldata, or Guard) at  *
     * the center of every lead node. By default all cells will be assigned   *
     * to last rank                                                           */
    int dataInserted = 0;
	it = tree.begin();
	while ( it!=tree.end() ) {
		if ( it->IsLeaf() ) {
            ++dataInserted;
			it->InsertData( it->center, &(data[dataInserted]) );
        }
        ++it;
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
    Octree::OctreeToSvg<int,SIMDIM> svgoutput( tree, std::string("Octree_") + sOrdering.str() + std::string("worldsize_") + sWorldsize.str() );
    svgoutput.PrintGrid();

    /* Assign cells to all the processes */
    double cumulativeCosts = 0;
    int curRank = 0;
    int currentTime  = 0;
    double delay     = 1./64.; // 8 frames per second at max achievable
    double rankDelay = 1./64.;
    
	it = tree.begin( ORDERING );
    Octree::Octree<int,SIMDIM>::iterator it0 = tree.begin();
    Octree::Octree<int,SIMDIM>::iterator it1 = tree.begin();
	while ( it!=tree.end() ) {
		if ( it->IsLeaf() ) {
			cumulativeCosts += 1. / it->getSize().min();
            if ( cumulativeCosts >= optimalCosts and curRank != worldsize-1) {
                cumulativeCosts = 1. / it->getSize().min();
                curRank++;
            }
            *(it->getDataPtr(0).object) = curRank;

            /* Graphical output of traversal line */
            it0 = it1;
            it1 = it;
            if ( it0->IsLeaf() and it1->IsLeaf() ) {
                size_t id = reinterpret_cast<size_t>( it0->getDataPtr(0).object );
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
        ++it;
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
            costs[ *(it->getDataPtr(0).object) ] += 1. / it->getSize().min();
            int nNeighbors = 0;
            int nLeavesOnOtherNodes = 0;
            int thisRank = *(it->getDataPtr(0).object);

            VecI dir[4];
            dir[0][0]=+1; dir[0][1]= 0;
            dir[1][0]=-1; dir[1][1]= 0;
            dir[2][0]= 0; dir[2][1]=+1;
            dir[3][0]= 0; dir[3][1]=-1;

            for ( int lindir = 0; lindir < 4; lindir++ ) {
                Octree::Node<int,SIMDIM> * neighbor = it->getNeighbor( dir[lindir] );
                if ( neighbor == NULL )
                    continue;

                nNeighbors += neighbor->countLeaves();

                Octree::Octree<int,SIMDIM>::iterator itn = neighbor->begin();
                while ( itn != neighbor->end() ) {
                    if ( itn->IsLeaf() ) {
                        int neighborRank = *(itn->getDataPtr(0).object);
                        if ( thisRank != neighborRank )
                            nLeavesOnOtherNodes++;
                    }
                    ++itn;
                }
            }

            totalTraffic += bytesPerCell * nNeighbors;
            interTraffic += bytesPerCell * nLeavesOnOtherNodes;

            /*tout << it->center << " needs data from "
                 << nNeighbors << " neighbors. " << nLeavesOnOtherNodes << " of those are not on this process\n"; */
		}
        ++it;
	}

    tout << "World Size      : " << worldsize << "\n";
    tout << "Number of Cells : " << tree.root->countLeaves() << "\n";
    for (int i=0; i<worldsize; i++) {
        tout << "Cost assigned to rank " << i << " is " << costs[i] << "\n";
    }
    tout << "Total data to be read from neighbors : " << totalTraffic << "\n";
    tout << "Data which needs to be communicated  : " << interTraffic << "\n";
    
    resultsFile << worldsize << " " << totalTraffic << " " << interTraffic << "\n";

    delete[] costs;
    free(data);

} /* Run this Programm for several world sizes ! */
    resultsFile.close();
} /* Run this Programm for several Ordering Methods ! */

}