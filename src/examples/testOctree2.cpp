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
//#include <omp.h> // not needed for pragmas
#include "CompileTime.h"
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
    /* Create timestamp and basefolder for filenames */
    time_t tmp_time = time(0);
    struct tm * now = localtime( &tmp_time );
    std::stringstream basefolder, timestamp, basefilename;
    timestamp  << 1900 + now->tm_year << "-" << std::setfill('0')
               << std::setw(2) << 1 + now->tm_mon << "-"
               << std::setw(2) << now->tm_mday << "_"
               << std::setw(2) << now->tm_hour << "-"
               << std::setw(2) << now->tm_min  << "-"
               << std::setw(2) << now->tm_sec;
    basefolder << "output/" << timestamp.str();
    boost::filesystem::create_directory(
        boost::filesystem::absolute(basefolder.str()) );
    basefilename << basefolder.str() << "/" << timestamp.str() << "_";

    tout.Open( basefilename.str() + std::string("out"), -1 );
    terr.Open( basefilename.str() + std::string("err"), -1 );
    tout << "Write logs to: " << basefolder.str() << "\n";

    bool PRINT_SVG = true;
    int ORDERING = 0;
    int NUMBER_OF_WORLDSIZES = 400;
    while ( true ) {
        static struct option long_options[] = {
            {"init-refinement"     , required_argument, 0, 'i'},
            {"max-refinement"      , required_argument, 0, 'm'},
            {"number-of-cells"     , required_argument, 0, 'n'},
            {"octree-setup"        , required_argument, 0, 'o'},
            {"svg"                 , required_argument, 0, 's'},
            {"ordering"            , required_argument, 0, 'r'},
            {"number-of-worldsizes", required_argument, 0, 'n'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "i:m:n:o:r:s:", long_options, &option_index);

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
            case 'r':
                ORDERING = atoi(optarg);
                break;
            case 'n':
                //if ( strstr( argv[optind-2], "--number-of-worldsizes\0" ) != NULL )
                    NUMBER_OF_WORLDSIZES = atoi(optarg);
                break;
            case 's':
                if ( strstr( argv[optind-2], "--svg\0" ) != NULL ) {
                    if ( optarg[0] == '0' )
                        PRINT_SVG = false;
                    else if ( optarg[0] == '1' )
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
    tout << "PRINT_SVG: " << PRINT_SVG << "\n";
    tout << "NUMBER_OF_WORLDSIZES: " << NUMBER_OF_WORLDSIZES << "\n";

/* Run this Programm for several Ordering Methods ! */

//for ( int ORDERING = 0; ORDERING <= 2; ++ORDERING ) {
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
    filename << basefolder.str() << "/" << (SIMDIM == 2 ? "Quadtree" : "Octree")
             << "-Setup-" << OCTREE_SETUP
             << "_Initial-" << INITIAL_OCTREE_REFINEMENT << "_Max-Refinement-"
             << MAX_OCTREE_REFINEMENT << "_" << sOrdering.str()
             << "_Ordering";
    resultsFile.open( filename.str() + std::string(".dat") );
    resultsFile << "# worldsize totalTraffic interTraffic messageCount\n" << std::flush;

    VecD globSize   = VecD(NUMBER_OF_CELLS);
    VecD globCenter = 0;
    typedef Octree::Octree<SIMDIM> OctreeType;
    Octree::Octree<SIMDIM> tree( globCenter, globSize );

    /* refine all cells to some value */
    if ( OCTREE_SETUP == 1 or OCTREE_SETUP == 6 or OCTREE_SETUP == 7 ) {
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
    if ( OCTREE_SETUP == 6 or OCTREE_SETUP == 7 ) {
    /**************************************************************************/
    /* (1) Setup Octree Refinement ********************************************/
    /**************************************************************************/

    /********************** Refine spherical boundaries ***********************/
    assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
    for ( int lvl=INITIAL_OCTREE_REFINEMENT; lvl<MAX_OCTREE_REFINEMENT; lvl++) {
        std::stack<OctreeType::Node*> torefine;
        for ( OctreeType::iterator it = tree.begin(); it != tree.end(); ++it )
            if ( it->IsLeaf() and ( M - it->center*tree.size ).norm() <= R ) {
                bool oneneighboroutside = false;
                for ( int i = 0; i < pow(3,SIMDIM); ++i ) {
                    OctreeType::Node * neighbor =
                        it->getNeighbor( getDirectionVector<SIMDIM>(i), VecI(0) );
                    if ( neighbor == NULL ) {
                        oneneighboroutside = true;
                        continue;
                    }
                    if ( ( M - neighbor->center*tree.size ).norm() > R ) {
                        oneneighboroutside = true;
                        break;
                    }
                }
                if ( oneneighboroutside )
                    torefine.push( &(*it) );
            }
        while ( not torefine.empty() ) {
            torefine.top()->GrowUp();
            torefine.pop();
        }
    }
    } // OCTREE_SETUP == 6

    int NValues = tree.root->countLeaves();
    tout << "Tree-Integrity  : " << tree.CheckIntegrity() << "\n";
    tout << "Number of Leaves: " << NValues               << "\n";
    clock_t bufferstart = clock();
    Octree::Octree<SIMDIM>::iterator itordering = tree.begin( ORDERING );
    tout << "Buffering traversal took " << (double) (clock() - bufferstart) / CLOCKS_PER_SEC << "s\n";

    /* allocate data (which stores assigned ranks) to which pointers will be  *
     * given to octree. And default it to the last rank                       */
    int * data = (int*) malloc( sizeof(int)*NValues );
    for (int i=0; i<NValues; ++i)
        data[i] = 0;

    /* Insert testData (later YeeCell-Data or Absorbercelldata, or Guard) at  *
     * the center of every leaf node. By default all cells will be assigned   *
     * to last rank                                                           */
    int dataInserted = 0;
    for ( Octree::Octree<SIMDIM>::iterator it = tree.begin(); it!=tree.end();
    ++it ) if ( it->IsLeaf() ) {
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

    /* Buffer traversal Order, to save time! */
    /*typedef Octree::Octree<SIMDIM>::Node Node;
    Node ** traversalOrder = (Node**) malloc( sizeof(Node**) * NValues );
    int curLeaf = 0;
    for ( Octree::Octree<SIMDIM>::iterator it = tree.begin( ORDERING );
          it != tree.end(); ++it ) if ( it->IsLeaf() )
    {
        assert( curLeaf < NValues );
        traversalOrder[ curLeaf++ ] = &(*it);
        tout << "curLeaf: " << curLeaf << "\n";
    }*/

int LINUPTO = 16;
if ( LINUPTO > NUMBER_OF_WORLDSIZES )
    LINUPTO = NUMBER_OF_WORLDSIZES;

int * wsizes = (int*) malloc( NUMBER_OF_WORLDSIZES*sizeof(int) );

for (int i=0; i<LINUPTO; i++)
    wsizes[i] = i+1;
/* a*exp[b*xe] = NValues; xe=NUMBER_OF_WORLDSIZES-1-LINUPTO
 * a*exp[b*xa] = LINUPTO; xa=0 => a=LINUPTO */
double bcoeff = log(double(NValues)/double(LINUPTO))/double(NUMBER_OF_WORLDSIZES-LINUPTO-1);
int lastindex = LINUPTO-1;
for (int i=0; i < NUMBER_OF_WORLDSIZES - LINUPTO; i++) {
    int curindex = LINUPTO+i;
    int curvalue = int(LINUPTO * exp(bcoeff*i));
    if ( wsizes[lastindex] >= curvalue )
        curvalue = wsizes[lastindex]+1;
    wsizes[curindex] = curvalue;
    lastindex = curindex;
}
tout << "wsizes:\n";
for (int i=0; i<NUMBER_OF_WORLDSIZES; i+=20)
    tout << "  " << wsizes[i] << "\n";


int worldsize = wsizes[NUMBER_OF_WORLDSIZES-1];
for ( int iw = NUMBER_OF_WORLDSIZES-1; iw >= 0; worldsize = wsizes[--iw] ) {
    clock_t twstart = clock();

    double optimalCosts = totalCosts / double(worldsize);

    /* Print tree to SVG */
    std::stringstream sWorldsize;
    sWorldsize << filename.str();
    if (PRINT_SVG)
        sWorldsize << "_worldsize-" << worldsize << ".svg";
    Octree::OctreeToSvg<SIMDIM> svgoutput( tree, sWorldsize.str(), false );
    if (PRINT_SVG) {
        tout << "Open: " << sWorldsize.str() << "\n";
        if ( SIMDIM == 2 ) {
            svgoutput.PrintGrid();
        }
    }

    /* Assign cells to all the processes */
    double cumulativeCosts = 0;
    int curRank = 0;
    double currentTime = 0;
    double rankDelay = 8./64.;// 8 frames per second at max achievable with SVG
    int curCell = 0;

    //tout << "Tree:\n" << tree << "\n\n";

    /*
        case 3: sOrdering << "Column-wise_" ; break;
        case 4: sOrdering << "Four-Color-Theorem_" ; break;
    */
    /* Using optimalCosts, assign cells to processes */
    Octree::Octree<SIMDIM>::iterator it0 = tree.begin();
    Octree::Octree<SIMDIM>::iterator it1 = tree.begin();
    for ( Octree::Octree<SIMDIM>::iterator it = itordering.begin();
          it != tree.end(); ++it ) if ( it->IsLeaf() )
    {
        curRank = int( cumulativeCosts / optimalCosts );
        cumulativeCosts += weightfunc(*it);
        /* only should happen if last cell(s) in traversal has 0 weighting */
        assert( curRank < worldsize );
        *((int*)it->data[0]) = curRank;

        if (PRINT_SVG) {
            /* Graphical output of traversal line */
            it0 = it1;
            it1 = it;
            if ( it0->IsLeaf() and it1->IsLeaf() ) {
                size_t id = reinterpret_cast<size_t>( it0->data[0] );
                Vec<double,2> r0 = svgoutput.convertToImageCoordinates( it0->center );
                Vec<double,2> r1 = svgoutput.convertToImageCoordinates( it1->center );
                /* Spawn invisible line element */
                svgoutput.out
                  << "<line"                                               "\n"
                  << " id=\"path" << id << "\""                            "\n"
                  << " x1=\"" << r0[0] << "px\" y1=\"" << r0[1] << "px\""  "\n"
                  << " x2=\"" << r1[0] << "px\" y2=\"" << r1[1] << "px\""  "\n"
                  << " style=\"stroke:none;stroke-width:3px\""             "\n"
                  << "/>"                                                  "\n";

                int ivalue = 4 * curCell;
                int ncolor = int( floor( ivalue / double(NValues) ) );
                assert( size_t(ncolor) < Colors::Own1Length );
                /*ncolor = int( floor( ivalue / double(NValues) ) );
                uint32_t rgb0 = Colors::BuPu[ncolor];
                uint32_t rgb1 = Colors::BuPu[ncolor+1];
                uint32_t rgb  = rgb0 + ( (rgb1 - rgb0) * (ivalue % NValues)) / NValues; //Colors::BuPu[ncolor];
                int r  = Colors::getRed  ( rgb );
                int g  = Colors::getGreen( rgb );
                int b  = Colors::getBlue ( rgb );*/
                Vec<double,3> rgb0 = Colors::getColorVector( Colors::Own1[ncolor+0] );
                Vec<double,3> rgb1 = Colors::getColorVector( Colors::Own1[ncolor+1] );
                Vec<double,3> rgb  = rgb0 + (rgb1 - rgb0) * (ivalue % NValues) / NValues;
                int r  = (int) rgb[0];
                int g  = (int) rgb[1];
                int b  = (int) rgb[2];

                /* Animate line element to become visible after some time */
                svgoutput.out
                  << "<set"                                            "\n"
                  << " xlink:href=\"#path" << id << "\""               "\n"
                  << " attributeName=\"stroke\""                       "\n"
                  << " begin=\"" << currentTime << "s\"" "\n"
                  << " to   =\"#008000\""                              "\n"
                  //<< " to   =\"rgb(" << r << "," << g << "," << b << ")\"" "\n"
                  << "/>"                                              "\n";
                currentTime += rankDelay;
            }
            if ( SIMDIM == 2 ) {
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
        curCell++;
    }
    svgoutput.close();

    tout << "Printing SVG took " << double(clock()-twstart)/CLOCKS_PER_SEC << "s\n";

    /* Count Neighbors intra- and interprocessdata to transmit */
    double * costs = new double[worldsize];
    for (int i=0; i<worldsize; ++i)
        costs[i] = 0;

    int interTraffic    = 0;
    int totalTraffic    = 0;
    const int bytesPerCell = 256;//(6+4)*8;

    /* Count foreign and overall neighbors of cells, also count/distribute costs */
    for ( Octree::Octree<SIMDIM>::iterator it = itordering.begin();
          it != it.end(); ++it ) if ( it->IsLeaf() )
    {
        costs[ *((int*)it->data[0]) ] += weightfunc(*it);
        int nNeighbors          = 0;
        int nLeavesOnOtherNodes = 0;
        int thisRank            = *((int*)it->data[0]);
        int tmpinterTraffic     = interTraffic;
        int tmptotalTraffic     = totalTraffic;

        /* Iterate over all directions and get neighbors there */
        for ( int lindir = 0; lindir < pow(3,SIMDIM); lindir++ )  {
            VecI dir = getDirectionVector<SIMDIM>(lindir);
            /* Only count direct neighbors, no diagonal ones. Also exclude    *
             * lindir == 0, because that is the Node itself !                 */
            if ( dir.abs().sum() != 1 )
                continue;
            std::list<OctreeType::Node*> neighbors = it->getNeighbors( dir, VecI(0) );

            /* Iterate over all neighbors in a given direction */
            for ( std::list<OctreeType::Node*>::iterator itn = neighbors.begin();
                  itn != neighbors.end(); ++itn ) if ( (*itn)->IsLeaf() )
            {
                assert( (*itn)->data[0] != NULL );
                nNeighbors += 1;
                const int lvldiff = abs( (*itn)->getLevel() - it->getLevel() );
                totalTraffic += bytesPerCell / int(pow( 2, lvldiff ));

                const int neighborRank = *( (int*) (*itn)->data[0] );
                if ( thisRank != neighborRank ) {
                    nLeavesOnOtherNodes++;
                    interTraffic += bytesPerCell / int(pow( 2, lvldiff ));
                }
            }
        }
    }
    resultsFile << worldsize << " " << totalTraffic << " " << interTraffic << "\n" << std::flush;

    delete[] costs;

    clock_t twend = clock();
    tout << "[" << sOrdering.str() << "] World Size : " << worldsize
         << ", Total Costs: " << totalCosts
         << " => Optimal Costs: " << optimalCosts
         << " took " << double(twend-twstart)/CLOCKS_PER_SEC << "s\n"
         << std::flush;

} /* Run this Programm for several world sizes ! */
    free(data);
    free(wsizes);
    //free(traversalOrder);
    resultsFile.close();

//} /* Run this Programm for several Ordering Methods ! */
}
