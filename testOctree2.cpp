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

    int worldsize = 2;

    VecD size(100), center(0);
    Octree::Octree<int,SIMDIM> tree( center, size );

    /* refine all cells to some value */
    tout << "Initial homogenous Refinement\n";
	Octree::Octree<int,SIMDIM>::iterator it;
	int maxRecursion = 1;
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
    /*tout << "Refinement of a subsection\n";
	it = tree.begin();
	while ( it!=tree.end() ) {
		if ( it->IsLeaf() and it->center.norm() < 0.25 )
			(it++)->GrowUp();
		else
			++it;
	}*/
    
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
    tree.FindLeafContainingPos( pos*tree.size )->GrowUp();

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
    Octree::OctreeToSvg<int,SIMDIM> svgoutput( tree, std::string("octree") );
    svgoutput.PrintGrid();

    /* Assign cells to all the processes */
    double cumulativeCosts = 0;
    int curRank = 0;
	it = tree.begin( Octree::Ordering::GrayCode );
    Octree::Octree<int,SIMDIM>::iterator it0 = tree.begin();
    Octree::Octree<int,SIMDIM>::iterator it1 = tree.begin();
    int currentTime  = 0;
    double delay     = 0.125; // 8 frames per second at max achievable
    double rankDelay = 0.25;
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
                std::cout << "Leaf at " << it0->center << "and at " << it1->center << "\n";
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
                  << " to   =\"#d000ff\""                              "\n"
                  << "/>"                                              "\n";
                currentTime += rankDelay / delay;
            }
            
            /* Graphical output of and cell-rank-mapping */
            
		}
        ++it;
	}

    /* Count Neighbors intra- and interprocessdata to transmit */
    //const int bytesPerCell = (6+4)*8;
    double * costs = new double[worldsize];
	it = tree.begin();
	while ( it!=tree.end() ) {
		if ( it->IsLeaf() ) {
            costs[ *(it->getDataPtr(0).object) ] += 1. / it->getSize().min();
            /* this should only work if the neighbors are smaller or equal size */
            int nNeighbors = 0;
            VecI dir(0);
            dir[0]=+1; dir[1]= 0;
            nNeighbors += ( it->getNeighbor(dir) == NULL ) ? 0 : it->getNeighbor(dir)->countLeaves();
            dir[0]=-1; dir[1]= 0;
            nNeighbors += ( it->getNeighbor(dir) == NULL ) ? 0 : it->getNeighbor(dir)->countLeaves();
            dir[0]= 0; dir[1]=+1;
            nNeighbors += ( it->getNeighbor(dir) == NULL ) ? 0 : it->getNeighbor(dir)->countLeaves();
            dir[0]= 0; dir[1]=-1;
            nNeighbors += ( it->getNeighbor(dir) == NULL ) ? 0 : it->getNeighbor(dir)->countLeaves();

            std::cout << "Leaf at " << it->center << " needs data from "
                      << nNeighbors << " neighbors:\n";
            dir[0]=+1; dir[1]= 0;
            /*std::cout << "   " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->center ) << " -> " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->countLeaves() ) << "\n";
            dir[0]=-1; dir[1]= 0;
            std::cout << "   " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->center ) << " -> " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->countLeaves() ) << "\n";
            dir[0]= 0; dir[1]=+1;
            std::cout << "   " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->center ) << " -> " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->countLeaves() ) << "\n";
            dir[0]= 0; dir[1]=-1;
            std::cout << "   " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->center ) << " -> " << ( ( it->getNeighbor(dir) == NULL ) ? 333 : it->getNeighbor(dir)->countLeaves() ) << "\n";
            */
            /* only correct up to a certain cellsize !!! */
            /*const double smallestCellsize = 1e-8; // >> DBL_EPSILON
            VecD neighborPos( it->center );
            neighborPos[0] += 0.5*smallestCellsize;
            tree->FindLeafContaining( neighborPos );
            neighborPos = it->center;
            neighborPos[0] -= 0.5*smallestCellsize;
            tree->FindLeafContaining( neighborPos );
            neighborPos = it->center;
            neighborPos[1] += 0.5*smallestCellsize;
            tree->FindLeafContaining( neighborPos );
            neighborPos = it->center;
            neighborPos[1] -= 0.5*smallestCellsize; */
		}
        ++it;
	}

    for (int i=0; i<worldsize; i++) {
        std::cout << "Cost assigned to rank " << i << " is " << costs[i] << "\n";
    }

    delete[] costs;
    free(data);
}