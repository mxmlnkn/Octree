/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe
make testOctree && ./testOctree.exe

*/

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc, rand
#include <random>   // normal_distribution
#include <ctime>    // time
#include "getopt.h"
#include <ctime>
#include <boost/filesystem.hpp>

namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime

#include "paramset/Parameters_2015-01-26_testOctree.cpp"
#include "math/TVector.h"
#include "octree/Octree.h"
#include "octree/OctreeToSvg.h"
#include "teestream/TeeStream.h"
#include "Colors.h"
#include "Directions.h"

typedef Vec<double,SIMDIM> VecD;
typedef Vec<double,SIMDIM> VecI;

#define RANDOM_SEED 24756139



int main( int argc, char **argv )
{
    time_t t = time(0);
    struct tm * now = localtime( &t );
    std::stringstream basefolder;
    basefolder << "output/" << 1900 + now->tm_year << "-" << std::setfill('0')
               << std::setw(2) << 1 + now->tm_mon << "-"
               << std::setw(2) << now->tm_mday << "_"
               << std::setw(2) << now->tm_hour << "-"
               << std::setw(2) << now->tm_min;
    boost::filesystem::create_directory(
        boost::filesystem::absolute(basefolder.str()) );
    basefolder << "/";

    tout.Open( std::string("out"), -1, basefolder.str() );

    bool PRINT_SVG = true;
    while ( true ) {
        static struct option long_options[] = {
            {"init-refinement" , required_argument, 0, 'i'},
            {"max-refinement"  , required_argument, 0, 'm'},
            {"number-of-cells" , required_argument, 0, 'n'},
            {"octree-setup"    , required_argument, 0, 'o'},
            {"svg"             , required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "i:m:n:o:s:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
            case 'i':
                INITIAL_OCTREE_REFINEMENT = atoi(optarg);
                break;
            case 'm':
                MAX_OCTREE_REFINEMENT = atoi(optarg);
                break;
            case 'o':
                OCTREE_SETUP = atoi(optarg);
                break;
            case 's':
                if ( strstr( argv[optind-2], "--svg\0" ) != NULL ) {
                    if ( argv[optind-1][0] == '0' )
                        PRINT_SVG = false;
                    else if ( argv[optind-1][0] == '1' )
                        PRINT_SVG = true;
                }
                break;
            default:
                tout << "Wrong Parameters!\n";
                return 1;
        }
    }

tout << "OCTREE_SETUP : " << OCTREE_SETUP << "\n";
tout << "INITIAL_OCTREE_REFINEMENT: " << INITIAL_OCTREE_REFINEMENT << "\n";
tout << "MAX_OCTREE_REFINEMENT: " << MAX_OCTREE_REFINEMENT << "\n";
tout << "Print SVG: " << PRINT_SVG << "\n";

/* Run this Programm for several Ordering Methods ! */
for ( int ORDERING = 3; ORDERING <= 3; ++ORDERING ) {
/* Run this Programm for several world sizes ! */
    std::ofstream resultsFile;
    std::stringstream sOrdering, filename;
    switch (ORDERING) {
        case 0: sOrdering << "Morton"  ; break;
        case 1: sOrdering << "GrayCode"; break;
        case 2: sOrdering << "Hilbert" ; break;
        case 3: sOrdering << "Rows" ; break;
        case 4: sOrdering << "Four-Color-Theorem" ; break;
    }
    filename << basefolder.str() << (SIMDIM == 2 ? "Quadtree" : "Octree")
             << "-Setup-" << OCTREE_SETUP
             << "_Initial-" << INITIAL_OCTREE_REFINEMENT << "_Max-Refinement-"
             << MAX_OCTREE_REFINEMENT << "_" << sOrdering.str()
             << "_Ordering";
    resultsFile.open( filename.str() + std::string(".dat") );
    resultsFile << "# worldsize totalTraffic interTraffic messageCount\n" << std::flush;

const int NWSIZES = 400;
int NVALUESGUESSED = 93184;
int LINUPTO = 16;
int wsizes[NWSIZES];
for (int i=0; i<LINUPTO; i++)
    wsizes[i] = i+1;
/* a*exp[b*xe] = NVALUESGUESSED; xe=NWSIZES-1-LINUPTO
 * a*exp[b*xa] = LINUPTO; xa=0 => a=LINUPTO */
double bcoeff = log(double(NVALUESGUESSED)/double(LINUPTO))/double(NWSIZES-LINUPTO-1);
for (int i=0; i < NWSIZES - LINUPTO; i++)
    wsizes[LINUPTO+i] = int(LINUPTO * exp(bcoeff*i));
/*tout << "wsizes:\n";
for (int i=0; i<NWSIZES; i++)
    tout << wsizes[i] << "\n";*/


int worldsize = 1;
for ( int iw = 0; iw < NWSIZES; worldsize = wsizes[iw++] ) {
    tout << "World Size      : " << worldsize << "\n";

    VecD globSize   = VecD(NUMBER_OF_CELLS);
    VecD globCenter = 0;
    typedef Octree::Octree<SIMDIM> OctreeType;
    Octree::Octree<SIMDIM> tree( globCenter, globSize );

    /* refine all cells to some value */
    if ( OCTREE_SETUP == 1 ) {
        for ( int i=0; i < INITIAL_OCTREE_REFINEMENT; i++)
            for ( Octree::Octree<SIMDIM>::iterator it = tree.begin();
            it != it.end(); )
                if ( it->IsLeaf() and it->getLevel() == i )
                    (it++)->GrowUp();
                else
                    ++it;
    }
    /* refine a subsection ( 2D-Ball ) of the cells one more time */
    if ( OCTREE_SETUP == 2 ) {
        for ( Octree::Octree<SIMDIM>::iterator it = tree.begin();
        it != it.end(); ) if ( it->IsLeaf() and it->center.norm() < 0.4 )
            (it++)->GrowUp();
        else
            ++it;
        for ( Octree::Octree<SIMDIM>::iterator it = tree.begin();
        it != it.end(); ) if ( it->IsLeaf() and it->center.norm() < 0.25 )
            (it++)->GrowUp();
        else
            ++it;
    }
    /* refine upper right corner one more time */
    if ( OCTREE_SETUP == 3 ) {
        for ( Octree::Octree<SIMDIM>::iterator it = tree.begin();
        it != it.end(); ) if ( it->IsLeaf() and it->center[0] > 0 )
            (it++)->GrowUp();
        else
            ++it;
        VecD pos(0); pos[0]=0.375; pos[1]=-0.125;
        tree.FindLeafContainingPos( pos*tree.size )->GrowUp();
    }
    /* refine cells at fixed point many times */
    if ( OCTREE_SETUP == 4 ) {
        VecD pos(0); pos[0]=0.1; pos[1]=0.1;
        for ( int i = 0; i < MAX_OCTREE_REFINEMENT; ++i )
            tree.FindLeafContainingPos( pos*tree.size )->GrowUp();
    }
    /* refine some randomly chosen cells until target cell count reached */
    if ( OCTREE_SETUP == 5 ) {
        srand( RANDOM_SEED );
        int targetCells  = int(pow(pow(2,SIMDIM),INITIAL_OCTREE_REFINEMENT));
        int currentCells = 1;
        while ( currentCells < targetCells ) {
            VecD pos(0);
            pos[0] = double(rand()) / double(RAND_MAX);
            if ( SIMDIM >= 2 )
                pos[1] = double(rand()) / double(RAND_MAX);
            if ( SIMDIM >= 3 )
                pos[2] = double(rand()) / double(RAND_MAX);
            tree.FindLeafContainingPos( tree.center - 0.5*tree.size + pos*tree.size )->GrowUp();
            currentCells += int(pow(2,SIMDIM))-1;
        }
    }
    if ( OCTREE_SETUP == 6 ) {
    /**************************************************************************/
    /* (1) Setup Octree Refinement ********************************************/
    /**************************************************************************/

    /********* refine all cells to initial homogenous min-Refinement **********/
    for ( int lvl=0; lvl < INITIAL_OCTREE_REFINEMENT; lvl++) {
        for ( OctreeType::iterator it=tree.begin(); it != tree.end(); )
            if ( it->IsLeaf() and it->getLevel()==lvl )
                (it++)->GrowUp();
            else
                ++it;
    }
    /*********************** Refine certain boundaries ************************/
    assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
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
    }
    tout << "Tree-Integrity  : " << tree.CheckIntegrity() << "\n";
    tout << "Number of Leaves: " << tree.root->countLeaves() << "\n";

    int NValues = tree.root->countLeaves();

    /* allocate data (which stores assigned ranks) to which pointers will be  *
     * given to octree. And default it to the last rank                       */
    int * data = (int*) malloc( sizeof(int)*NValues );
    for (int i=0; i<NValues; ++i)
        data[i] = worldsize-1;

    /* Insert testData (later YeeCell-Data or Absorbercelldata, or Guard) at  *
     * the center of every leaf node. By default all cells will be assigned   *
     * to last rank                                                           */
    int dataInserted = 0;
    for ( Octree::Octree<SIMDIM>::iterator it = tree.begin(); it!=tree.end();
    ++it ) if ( it->IsLeaf() )
    {
        assert( dataInserted < NValues );
        assert( it->data.empty() );
        it->data.push_back( &(data[dataInserted]) );
        ++dataInserted;
    }

    /* cost function */
    struct T_WEIGHT_FUNC {
        double operator() (const Octree::Octree<SIMDIM>::Node & curnode) {
            return 1;
            //return 1. / node.getSize().min();
        }
    } weightfunc;

    /* Calculate total costs of all cells. Could be done when counting cells  */
    double totalCosts = 0;
    for ( Octree::Octree<SIMDIM>::iterator it=tree.begin(); it!=it.end(); ++it )
        if ( it->IsLeaf() )
            totalCosts += weightfunc(*it);
    double optimalCosts = totalCosts / double(worldsize);
    tout << "Total Costs: " << totalCosts << " => Optimal Costs: " << optimalCosts << std::endl;

    /* Print tree to SVG */
    std::stringstream sWorldsize;
    sWorldsize << filename.str();
    if (PRINT_SVG)
        sWorldsize << "_worldsize-" << worldsize << ".svg";
    Octree::OctreeToSvg<SIMDIM> svgoutput( tree, sWorldsize.str(), false );
    tout << "Open: " << sWorldsize.str() << "\n";
    if (PRINT_SVG)
        svgoutput.PrintGrid();

    /* Assign cells to all the processes */
    double cumulativeCosts = 0;
    int curRank = 0;
    double currentTime = 0;
    double rankDelay = 1./64.;// 8 frames per second at max achievable with SVG

    //tout << "Tree:\n" << tree << "\n\n";

    /*
        case 3: sOrdering << "Column-wise_" ; break;
        case 4: sOrdering << "Four-Color-Theorem_" ; break;
    */
    Octree::Octree<SIMDIM>::iterator it0 = tree.begin();
    Octree::Octree<SIMDIM>::iterator it1 = tree.begin();
    for ( Octree::Octree<SIMDIM>::iterator it = tree.begin( ORDERING );
          it != tree.end(); ++it ) if ( it->IsLeaf() )
    {
        cumulativeCosts += weightfunc(*it);
        if ( cumulativeCosts >= optimalCosts and curRank != worldsize-1) {
            cumulativeCosts = weightfunc(*it);
            curRank++;
        }
        *((int*)it->data[0]) = curRank;

        if (PRINT_SVG) {
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
                  << " begin=\"" << currentTime << "s\"" "\n"
                  << " to   =\"#008000\""                              "\n"
                  << "/>"                                              "\n";
                currentTime += rankDelay;
            }

            /* Graphical output of and cell-rank-mapping */
            int id = svgoutput.boxesDrawn.find( it->center )->second.id;
            int r  = Colors::getRed  ( Colors::BuPu[8 - curRank % 9] );
            int g  = Colors::getGreen( Colors::BuPu[8 - curRank % 9] );
            int b  = Colors::getBlue ( Colors::BuPu[8 - curRank % 9] );
            svgoutput.out
                << "<set"                                                "\n"
                << " xlink:href=\"#" << id << "\""                       "\n"
                << " attributeName=\"fill\""                             "\n"
                << " begin=\"" << currentTime << "s\""     "\n"
                << " to   =\"rgb(" << r << "," << g << "," << b << ")\"" "\n"
                << "/>"                                                  "\n";
        }
    }
    svgoutput.close();

    /* Count Neighbors intra- and interprocessdata to transmit */
    double * costs = new double[worldsize];
    for (int i=0; i<worldsize; ++i)
        costs[i] = 0;

    double interTraffic = 0;
    double totalTraffic = 0;
    const int bytesPerCell = 1;//(6+4)*8;

    for ( Octree::Octree<SIMDIM>::iterator it = tree.begin( ORDERING );
          it!=tree.end(); ++it ) if ( it->IsLeaf() )
    {
        costs[ *((int*)it->data[0]) ] += weightfunc(*it);
        int nNeighbors = 0;
        int nLeavesOnOtherNodes = 0;
        int thisRank = *((int*)it->data[0]);

        int lindirs[4] = {RIGHT,LEFT,TOP,BOTTOM};
        VecI dir[4];
        dir[0][0]=+1; dir[0][1]= 0;
        dir[1][0]=-1; dir[1][1]= 0;
        dir[2][0]= 0; dir[2][1]=+1;
        dir[3][0]= 0; dir[3][1]=-1;

        for ( int lindir = 0; lindir < 4; lindir++ ) {
            assert( getDirectionVector<SIMDIM>(lindirs[lindir]) == dir[lindir]);
            Octree::Node<SIMDIM> * neighbor = it->getNeighbor( dir[lindir], VecI(0) );
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

    /*tout << "Number of Cells : " << tree.root->countLeaves() << "\n";
    for (int i=0; i<worldsize; i++) {
        tout << "Cost assigned to rank " << i << " is " << costs[i] << "\n";
    }
    tout << "Total data to be read from neighbors : " << totalTraffic << "\n";
    tout << "Data which needs to be communicated  : " << interTraffic << "\n";*/

    resultsFile << worldsize << " " << totalTraffic << " " << interTraffic << "\n" << std::flush;

    delete[] costs;
    free(data);

} /* Run this Programm for several world sizes ! */
    resultsFile.close();
} /* Run this Programm for several Ordering Methods ! */
}
