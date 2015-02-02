/*

ToDo:
  - Check for energy conservation automatically (problematic for initial wave)
  - draw double slit into png in gray and shift wave from red-green to only red
  - calculate and draw into the png the intensity after the obstacle
  - for CELL_SIZE_X = 5.0 the wave doesn't move to the right :S?
  - why no tprev needed like in 1D wave equation ...
Done:
  - check courant criterium !!! (automatically)
  - Introduce SI and internal Units!
  - above is for CELL_SIZE_X = 1.0

*/

namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const int exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime

#include <iostream>
#include <cmath>    // sin
#include <cfloat>   // FLT_EPSILON
#include <cstdlib>  // malloc, srand, rand, RAND_MAX
#include "getopt.h"
#include <pngwriter.h>
#include <list>
#include "mpi.h"
#include "math/TVector.h"
#include "math/TVector.tpp"
#include "math/TBaseMatrix.h"
#include "teestream/TeeStream.h"
#include "paramset/Parameters_2015-01-16.cpp"
#include "YeeCell.h"
#include "YeeCellColors.h"
#include "YeeSolver.h"
#include "octree/Octree.h"
#include "octree/OctreeToSvg.h"
#define DEBUG_COMMUNICATOR 0
#include "Communicator.h"
#include "Colors.h"
#include "Sources.h"

#define DEBUG_MAIN_YEE 99

#define N_CELLS_X NUMBER_OF_CELLS_X // should be power of two, because of octree !
#define N_CELLS_Y NUMBER_OF_CELLS_Y // should be same as above, because octree isn't able to have a different amount of equal sized cells in two directions !!!
#define N_CELLS_Z NUMBER_OF_CELLS_Z


int main( int argc, char **argv )
{
    /* NUMBER_OF_PARTICLES_PER_CELL, BOUNDARY_CONDITION, SPECIES, PNG_INTERVAL  -> Watch out for dependent Variables in Parameters.cpp!!! :  NUMBER_OF_EONS_PER_CELL, NUMBER_OF_IONS_PER_CELL       */

    while ( true ) {
        static struct option long_options[] = {
            {"timesteps"       , required_argument, 0, 't'},
            {"init-refinement" , required_argument, 0, 'i'},
            {"max-refinement"  , required_argument, 0, 'm'},
            {"number-of-cells" , required_argument, 0, 'n'},
            {"octree-setup"    , required_argument, 0, 'o'},
            {"simulation-setup", required_argument, 0, 's'},
            {"png-interval"    , required_argument, 0, 'p'},
            {"wave-source"     , required_argument, 0, 'w'},
            {"absorber"        , required_argument, 0, 'a'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "t:i:m:n:o:s:p:w:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
            case 't':
                NUMBER_OF_STEPS = atoi(optarg);
                break;
            case 'i':
                INITIAL_OCTREE_REFINEMENT = atoi(optarg);
                break;
            case 'm':
                MAX_OCTREE_REFINEMENT = atoi(optarg);
                break;
            case 'n':
                for (int i=0; i<SIMDIM; i++) {
                    assert( argv[optind-1+i][0] != '-' );
                    NUMBER_OF_CELLS[i]  = atoi(argv[optind-1+i]);
                }
                optind += SIMDIM-1; // extra arguments taken
                NUMBER_OF_CELLS_X   = NUMBER_OF_CELLS[0];
                if (SIMDIM > 1)
                    NUMBER_OF_CELLS_Y = NUMBER_OF_CELLS[1];
                else
                    NUMBER_OF_CELLS_Y = 1;
                if (SIMDIM > 2)
                    NUMBER_OF_CELLS_Z = NUMBER_OF_CELLS[2];
                else
                    NUMBER_OF_CELLS_Z = 1;
                NUMBER_OF_PARTICLES = NUMBER_OF_PARTICLES_PER_CELL *
                    NUMBER_OF_CELLS_X * NUMBER_OF_CELLS_Y * NUMBER_OF_CELLS_Z;
                SIM_SIZE = Vec<double,SIMDIM>( NUMBER_OF_CELLS ) * CELL_SIZE;
                break;
            case 'o':
                OCTREE_SETUP = atoi(optarg);
                break;
            case 's':
                SIMULATION_SETUP = atoi(optarg);
                break;
            case 'p':
                PNG_INTERVAL = atoi(optarg);
                break;
            case 'w':
                WAVE_SPAWN_SETUP = atoi(optarg);
                if ( WAVE_SPAWN_SETUP == 3 ) {
                if ( SIMDIM == 2 and argv[optind+0][0] != '-' and argv[optind+1][0] != '-' )
                {
                    SPAWN_POS = Vec<double,SIMDIM>( atoi(argv[optind+0]), atoi(argv[optind+1]) );
                    optind += 2;
                }
                if ( SIMDIM == 2 and argv[optind+0][0] != '-' and argv[optind+1][0] != '-' )
                {
                    SPAWN_AREA_SIZE = Vec<double,SIMDIM>( atoi(argv[optind+0]), atoi(argv[optind+1]) );
                    optind += 2;
                }
                }
                break;
            case 'a':
                ABSORBER_STRENGTH = atoi(optarg);
                if ( argv[optind-2][0] != '-' ) {
                    ABSORBING_BORDER_THICKNESS = atoi(argv[optind-2]);
                    optind += 1;
                }
                break;
            default:
                abort();
        }
    }

    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    typedef Vec<double,SIMDIM> VecD;
    typedef Vec<int   ,SIMDIM> VecI;

    double tProgramStart = MPI_Wtime();

    /* Call (indirectly) Basic Communicator-, Octree- and File-Constructors */
    VecD cellSize(CELL_SIZE);
    VecD NCells( NUMBER_OF_CELLS );
    VecD globSize(cellSize*NCells), globCenter(0.5*globSize);
    typedef typename Octree::Octree<SIMDIM> OctreeType;
    OctreeType tree( globCenter, globSize );
    //typedef BaseMatrix<YeeCell,SIMDIM> CellMatrix;
    typedef OctreeCommunicator<SIMDIM,OctreeType,YeeCell> OctreeCommType;
    OctreeCommType combox(tree,VecI(true));
    typedef SimulationBox::SimulationBox<SIMDIM,YeeCell> OctCell;
    if ( combox.rank != 0) {
        std::cerr << "Disable tout verbosity on rank " << combox.rank << "\n";
        tout.verbosity = 0;
    }
    tout.Open("out",combox.rank);
    terr.Open("err",combox.rank);
    /* Initialize SVG output file */
    std::stringstream svgfn;
    svgfn << "Octree_worldsize-" << combox.worldsize << "_rank-" << combox.rank;
    Octree::OctreeToSvg<SIMDIM> svgoutput( tree, svgfn.str() );

    /**************************** Print Parameters ****************************/
    srand(RANDOM_SEED);
	const double S = SPEED_OF_LIGHT * DELTA_T * (1./( Vec<double, 2>(CELL_SIZE) ) ).norm();
    tout << "MUE0                 : " << MUE0                       << "\n";
    tout << "EPS0                 : " << EPS0                       << "\n";
    tout << "SPEED_OF_LIGHT       : " << SPEED_OF_LIGHT             << "\n";
    tout << "CELL_SIZE_SI         : " << CELL_SIZE_SI               << "\n";
    tout << "CELL_SIZE            : " << CELL_SIZE_SI / UNIT_LENGTH << "\n";
    tout << "PARTICLE_RADIUS      : " << PARTICLE_RADIUS            << "\n";
    tout << "NUMBER_OF_CELLS_X    : " << NUMBER_OF_CELLS_X          << "\n";
    tout << "NUMBER_OF_CELLS_Y    : " << NUMBER_OF_CELLS_Y          << "\n";
    tout << "NUMBER_OF_CELLS_Z    : " << NUMBER_OF_CELLS_Z          << "\n";
    tout << "DELTA_T_SI           : " << DELTA_T_SI                 << "\n";
    tout << "UNIT_ENERGY          : " << UNIT_ENERGY                << "\n";
    tout << "UNIT_MOMENTUM        : " << UNIT_MOMENTUM              << "\n";
    tout << "UNIT_ANGULAR_MOMENTUM: " << UNIT_ANGULAR_MOMENTUM      << "\n";
    tout << "ELECTRON_MASS        : " << ELECTRON_MASS              << "\n";
    tout << "S                    : " << S                          << "\n";
    const double xoverc = 1. / (SPEED_OF_LIGHT * (1./( Vec<double, 2>(CELL_SIZE) ) ).norm() );
	tout << "For stability in vacuum Delta_T=" << DELTA_T << " =< " << xoverc <<  "=DELTA_X/sqrt(2)/c_M\n";
    if ( xoverc < DELTA_T )
        tout << " NOT FULFILLED!!!\n";
    const double xovercm = 1. / (SPEED_OF_LIGHT/1.33 * (1./( Vec<double, 2>(CELL_SIZE) ) ).norm() );
	tout << "For stability in glass  Delta_T=" << DELTA_T << " =< " << xovercm << "=DELTA_X/sqrt(2)/c_M\n";
    if ( xovercm < DELTA_T )
        tout << " NOT FULFILLED!!!\n";
    tout << "NUMBER_OF_STEPS          : " << NUMBER_OF_STEPS            << "\n";
    tout << "INITIAL_OCTREE_REFINEMENT: " << INITIAL_OCTREE_REFINEMENT  << "\n";
    tout << "MAX_OCTREE_REFINEMENT    : " << MAX_OCTREE_REFINEMENT      << "\n";
    tout << "OCTREE_SETUP             : " << OCTREE_SETUP               << "\n";
    tout << "SIMULATION_SETUP         : " << SIMULATION_SETUP           << "\n";
    tout << "WAVE_SPAWN_SETUP         : " << WAVE_SPAWN_SETUP           << "\n";
    tout << "SPAWN_POS                : " << SPAWN_POS                  << "\n";
    tout << "SPAWN_AREA_SIZE          : " << SPAWN_AREA_SIZE            << "\n";
    tout << "\n";


    /**************************************************************************/
    /* (1) Setup Octree Refinement ********************************************/
    /**************************************************************************/

    /********* refine all cells to initial homogenous min-Refinement **********/
    if ( OCTREE_SETUP == 7 or OCTREE_SETUP == 6 or OCTREE_SETUP == 9 ) {
        for ( int lvl=0; lvl<INITIAL_OCTREE_REFINEMENT; lvl++) {
            for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it )
                if ( it->IsLeaf() and it->getLevel()==lvl ) it->GrowUp();
        }
    }
    /*********************** Refine certain boundaries ************************/
    if ( OCTREE_SETUP == 6 ) {
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
    /******* Just grow up one point in upper left area (minimal setup) ********/
    if ( OCTREE_SETUP == 7 ) {
        tree.root->GrowUp();
        tree.FindLeafContainingPos( 0.75*globSize )->GrowUp();
    }
    /******* Grow up Spawning Area ********/
    if ( (OCTREE_SETUP == 8 or OCTREE_SETUP == 6) and false ) {
        double minCellSizeY = tree.size[1] / pow(2.,MAX_OCTREE_REFINEMENT);
        for ( VecD pos = SPAWN_POS; pos < SPAWN_POS + VecD(SPAWN_AREA_SIZE); pos[1] += minCellSizeY ) {
            OctreeType::Node * curNode = tree.FindLeafContainingPos(pos);
            if ( curNode->getLevel() < MAX_OCTREE_REFINEMENT ) {
                curNode->GrowUp();
                pos[1] -= minCellSizeY;
            }
        }
    }
    /**************************************************************************/
    if ( OCTREE_SETUP == 7 ) {
        assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
        for ( int lvl=INITIAL_OCTREE_REFINEMENT; lvl<MAX_OCTREE_REFINEMENT; lvl++)
            for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it ) {
                bool insideLense = false;
                if ( it->IsLeaf() and it->getLevel()==lvl and insideLense )
                    it->GrowUp();
            }
    }
    
    tout << "Tree-Integrity: " << tree.CheckIntegrity() << "\n\n";
    svgoutput.PrintGrid();

    /**************************************************************************/
    /* (2) Distribute weighting and Octree cells to processes *****************/
    /**************************************************************************/;
    VecI cellsPerOctreeCell = VecD(NCells) / pow( 2, tree.getMinLevel() );
    /* the number of cells may have to be adjusted, but don't touch cellsize  *
     * and time step, as those are stability critical parameters !!           */
    for (int i=0; i < cellsPerOctreeCell.dim; ++i) {
        bool toAssert =  fmod( (VecD(NCells) / pow( 2, tree.getMinLevel()))[i], 1.0 ) == 0.0;
        if ( not toAssert )
            terr << "Please choose NUMBER_OF_CELLS in multiples of "
                 << pow( 2, tree.getMinLevel()) <<"!\n";
        assert( toAssert );
    }
    tout << "Initial Refinement            : " << INITIAL_OCTREE_REFINEMENT << "\n"
         << "Minimum refinement Level found: " << tree.getMinLevel() << "\n"
         << "Maximum refinement Level found: " << tree.getMaxLevel() << "\n"
         << "Cells in Octree               : " << tree.countLeaves() << "\n"
         << "Cells per Octree Cell         : " << cellsPerOctreeCell << "\n\n";
    combox.initCommData( cellsPerOctreeCell, GUARDSIZE, 1 /*timestepbuffer*/ );

    const int ORDERING = Octree::Ordering::Hilbert; // -> Parameters -> need to include Octree there :S ?
    combox.distributeCells( ORDERING );

    /* Graphical output of cell-rank-mapping */
    for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
    if ( it->IsLeaf() ) {
        int id = svgoutput.boxesDrawn.find( it->center )->second.id;
        int curRank = ((typename OctreeCommType::CommData *)it->data[OctreeCommType::COMM_HEADER_INDEX])->rank;
        int r  = Colors::getRed  ( Colors::BuPu[8 - curRank % 9] );
        int g  = Colors::getGreen( Colors::BuPu[8 - curRank % 9] );
        int b  = Colors::getBlue ( Colors::BuPu[8 - curRank % 9] );
        svgoutput.out
            << "<set"                                                "\n"
            << " xlink:href=\"#" << id << "\""                       "\n"
            << " attributeName=\"fill\""                             "\n"
            << " to   =\"rgb(" << r << "," << g << "," << b << ")\"" "\n"
            << "/>"                                                  "\n";
    }

    /* Graphical output of traversal line */
    double currentTime = 0;
    double delay     = 1./64.; // 8 frames per second at max achievable with SVG
    double rankDelay = 1./64.;
    typename OctreeType::iterator it0 = tree.begin();
    typename OctreeType::iterator it1 = tree.begin();
    for ( typename OctreeType::iterator it = tree.begin( ORDERING );
          it!=tree.end(); ++it ) if ( it->IsLeaf() )
    {
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
              << "<set"                                    "\n"
              << " xlink:href=\"#path" << id << "\""       "\n"
              << " attributeName=\"stroke\""               "\n"
              << " begin=\"" << delay*currentTime << "s\"" "\n"
              << " to   =\"#008000\""                      "\n"
              << "/>"                                      "\n";
            currentTime += rankDelay / delay;
        }
    }
    svgoutput.close();

    /* Set Core, Border and Guard to different test values */
    for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
    if ( it->IsLeaf() ) { if ( combox.rank == ((OctreeCommType::CommData*)it->
         data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
    {
        OctCell & data = *((OctCell*)it->data[combox.CELL_DATA_INDEX]);
        typename OctCell::IteratorType itm = data.getIterator( 0, SimulationBox::CORE );
        itm = data.getIterator( 0, SimulationBox::BORDER );
        for ( itm = itm.begin(); itm != itm.end(); ++itm ) {
            itm->E[X] = -1.0;
        }
    } else if ( it->data.size() > combox.CELL_DATA_INDEX ) {
        assert(false);
    }}

    combox.PrintPNG( 0, "Border_Grid", returnEandH );

    /**************************************************************************/
    /* (5) Setup Simulation Start Data like Material **************************/
    /**************************************************************************/
    for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
    if ( it->IsLeaf() ) if ( combox.rank == ((OctreeCommType::CommData*)it->
         data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
    {
        OctCell & simbox = *((OctCell*)it->data[combox.CELL_DATA_INDEX]);

        for ( typename OctCell::IteratorType itm = simbox.getIterator( 0,
              SimulationBox::GUARD ).begin(); itm != itm.end(); ++itm )
        {
            /* it traverses octree nodes, while itm traverses matrix cells! */
            itm->E       = 0;
            itm->H       = 0;
            itm->epsilon = EPS0;
            itm->mu      = MUE0;
            itm->sigmaE  = 0;
            itm->sigmaM  = 0;
        }

        for ( typename OctCell::IteratorType itm = simbox.getIterator( 0,
              SimulationBox::CORE + SimulationBox::BORDER ).begin();
              itm != itm.end(); ++itm )
        {
            /* it traverses octree nodes, while itm traverses matrix cells! */
            itm->E       = 0;
            itm->H       = 0;
            itm->epsilon = EPS0;
            itm->mu      = MUE0;
            itm->sigmaE  = 0;
            itm->sigmaM  = 0;

            /* cur Pos is in internal units, meaning curPos in [0,1.0], but   *
             * it still points to the center of the current cell, not the     *
             * lower left corner. Also it's in unitless physical units        */
            for ( int i=0; i < SIMDIM; ++i ) {
                assert( itm.icell[i] > 0 );
                assert( itm.icell[i] - itm.guardsize < cellsPerOctreeCell[i] );
                bool toAssert = (tree.root->size/2 + it->center - it->size/2)[i] >= 0;
                if ( ! toAssert )
                    terr << "Lowerleft of node at " << it->center << " sized " << it->center << " is out of bounds!!\n";
                assert( toAssert );
            }
            VecD curPos = it->center - it->size/2 + (VecD(itm.icell) - itm.guardsize + 0.5) /
                          VecD(cellsPerOctreeCell) * it->size;
            curPos = tree.toGlobalCoords( curPos );
            curPos = simbox.getGlobalPosition(itm);
            assert( itm.ncells - 2*itm.guardsize == cellsPerOctreeCell );

            /********** Result 009: absorbing Material on right side **********/
            if ( SIMULATION_SETUP == 1 ) {
                if ( curPos[0] > 0.2 ) {
                    /* For INF instead of 2e8 it diverges! -> characteristic *
                     * wave length for exponential decay -> 0 for INF ?       */
                    itm->sigmaE  = ABSORBER_STRENGTH;
                    itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
                }
            }
            /*************** Result 014: broken total reflexion ***************/
            /* two glass plates with small vacuum/air slit inbetween          */
            /******************************************************************/
            if ( SIMULATION_SETUP == 2 ) {
                const double n  = 1.33; // = sqrt( eps_r * mue_r )
                if ( curPos[X] < 0.333*tree.size[X]
                or   curPos[X] > 0.350*tree.size[X]  )
                    itm->epsilon = EPS0 * n*n;
            }
            /********** Perfectly reflecting barrier with one slit ************/
            if ( SIMULATION_SETUP == 3 ) {
                const double wy = LAMBDA;
                if ( curPos[X] > LAMBDA and curPos[X] < 2*LAMBDA )
                if ( curPos[Y] > 0.5 - wy/2 and curPos[Y] < 0.5 + wy/2 ) {
                        itm->epsilon  = INF;//2*EPS0;
                        //itm->mu       = INF;
                }
            }
            /************************* Circular Lense *************************/
            if ( SIMULATION_SETUP == 4 ) {
                const double n  = 1.33; // = sqrt( eps_r * mue_r )
                if ( (curPos - M).norm() < R )
                    itm->epsilon = EPS0 * n*n;
                if ( (curPos[X] - 0 < ABSORBING_BORDER_THICKNESS) or
                     (tree.center[X] + tree.size[X]/2 - curPos[X] < ABSORBING_BORDER_THICKNESS) or
                     (tree.center[Y] + tree.size[Y]/2 - curPos[Y] < ABSORBING_BORDER_THICKNESS) or
                     (curPos[Y] - 0 < ABSORBING_BORDER_THICKNESS)
                   )
                {
                    /* For INF instead of 2e8 it diverges! -> characteristic *
                     * wave length for exponential decay -> 0 for INF ?       */
                    itm->sigmaE  = ABSORBER_STRENGTH;
                    itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
                }
                /* waveguide absorbers for WAVE_SPAWN_SETUP 3 */
                if ( ( (curPos >= SPAWN_POS + VecD( 0.0, SPAWN_AREA_SIZE[1] ) ) and
                       (curPos <  SPAWN_POS + VecD( 0.0, SPAWN_AREA_SIZE[1] ) + WAVE_GUIDE_SIZE ) ) or
                     ( (curPos >= SPAWN_POS - VecD( 0.0, WAVE_GUIDE_SIZE[1] ) ) and
                       (curPos < SPAWN_POS + VecD( WAVE_GUIDE_SIZE[0], 0.0 ) ) )
                   )
               {
                    itm->sigmaE  = ABSORBER_STRENGTH;
                    itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
               }
            }
            /************************ Parabolic Mirror ************************/
            if ( SIMULATION_SETUP == 5 ) {
                const double xcalc = pow( curPos[1] - MIRROR_CENTER[1], 2.0 ) / ( 4.0*FOCAL_LENGTH );
                if ( curPos[0] < xcalc + MIRROR_CENTER[0] ) {
                    itm->sigmaE  = INF;
                    itm->sigmaM  = INF / SPEED_OF_LIGHT * MUE0/EPS0;
               }
            }
            /********************** Gaussian Wave Pulse ***********************/
            if ( WAVE_SPAWN_SETUP == 5 ) {
                /**************************************************************
                 * Sine plane Wave going to Direction alpha and beginning     *
                 * line going through pos0  y                                 *
                 *               ^                                            *
                 *               | \     e.g. p0 ( line includes p0! )        *
                 *               |  --  /                                     *
                 *               |    \     alpha                             *
                 *               |     --  /                                  *
                 *               |_______\__________ x                        *
                 **************************************************************/
                double alpha    = 0 / 360. * 2*M_PI; //0.9*asin(1./n); // radian
                double NUMERICAL_SLOWING = 0.996797;
                double n                 = 1.0;
                double cM                = SPEED_OF_LIGHT*NUMERICAL_SLOWING/n;
                double T        = LAMBDA / cM;
                double kx       = 2*M_PI/LAMBDA * cos(alpha);
                double ky       = 2*M_PI/LAMBDA * sin(alpha);
                SPAWN_POS       = VecD( 0.1*SIM_SIZE[0], 0.75*SIM_SIZE[1] );
                SPAWN_AREA_SIZE = VecD( SIM_SIZE[0]/32, SIM_SIZE[1]/64 );
                assert( SPAWN_POS[1]+SPAWN_AREA_SIZE[1] < tree.center[1] + tree.size[1] );
                static bool firstIteration = true;
                if (firstIteration == true) {
                    firstIteration = false;
                    tout << "Spawn Position: " << SPAWN_POS << ", Size of Pulse: " << SPAWN_AREA_SIZE << "\n";
                }
                /* B = µH = E/c */
                /*itm->H[Y] = -1.0 / ( cM * MUE0 ) * 
                    TIME_SPAWN_FUNCTIONS::sinewave2d( 
                    T, 0.0, kx, curPos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
                    simbox.cellsize[0]), ky, curPos[1]-SPAWN_POS[1] ) *
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );*/
                itm->H[Y] = - TIME_SPAWN_FUNCTIONS::GaussianBeam( 
                    curPos[0] - SPAWN_POS[0], curPos[1] - SPAWN_POS[1], 0, 
                    SPAWN_AREA_SIZE[1], LAMBDA, T, 1 ) / ( cM * MUE0 ) * 
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] )* 
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
                /* E is half a time step later ! */
                /*itm->E[Z] = TIME_SPAWN_FUNCTIONS::sinewave2d( 
                    T, 0.5*DELTA_T, kx, curPos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
                    simbox.cellsize[0]), ky, curPos[1]-SPAWN_POS[1] ) *
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );*/
                itm->E[Z] = TIME_SPAWN_FUNCTIONS::GaussianBeam( curPos[0] - SPAWN_POS[0],
                    curPos[1] - SPAWN_POS[1], 0.5, SPAWN_AREA_SIZE[1], LAMBDA, T, 1 ) *
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] )*
                    TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
                if ( std::abs(itm->E[Z]) < 0.001) {
                    itm->E[Z] = 0;
                    itm->H[Y] = 0;
                }
            }
            /******* Two Gaussian Wave Pulses interfering destructively *******/
            if ( WAVE_SPAWN_SETUP == 6 ) {
                /**************************************************************
                 * Sine plane Wave going to Direction alpha and beginning     *
                 * line going through pos0  y                                 *
                 *               ^                                            *
                 *               | \     e.g. p0 ( line includes p0! )        *
                 *               |  --  /                                     *
                 *               |    \     alpha                             *
                 *               |     --  /                                  *
                 *               |_______\__________ x                        *
                 **************************************************************/
                const int NUM_WAVES      = 2;
                double alpha[NUM_WAVES]  = { 0, M_PI };
                double NUMERICAL_SLOWING = 0.996797;
                double n                 = 1.0;
                double cM                = SPEED_OF_LIGHT*NUMERICAL_SLOWING/n;
                VecD SPAWNPOS[NUM_WAVES] = { VecD( 0.1*SIM_SIZE[0], 0.5*SIM_SIZE[1] ),
                                             VecD( 0.9*SIM_SIZE[0], 0.5*SIM_SIZE[1] ) };
                VecD SPAWNAREASIZE[NUM_WAVES] = { VecD( SIM_SIZE[0]/32, SIM_SIZE[1]/64 ),
                                                  VecD( SIM_SIZE[0]/32, SIM_SIZE[1]/64 ) };
                assert( SPAWN_POS[1]+SPAWN_AREA_SIZE[1] < tree.center[1] + tree.size[1] );
                static bool firstIteration = true;
                if (firstIteration == true) {
                    firstIteration = false;
                    tout << "Spawn Position: " << SPAWN_POS << ", Size of Pulse: " << SPAWN_AREA_SIZE << "\n";
                }
                for (int iWave = 0; iWave < NUM_WAVES; iWave++) {
                    double kx = 2*M_PI/LAMBDA * cos(alpha[iWave]);
                    double ky = 2*M_PI/LAMBDA * sin(alpha[iWave]);
                    double T  = LAMBDA / cM;
                    /* B = µH = E/c */
                    itm->H[Y] = -1.0 / ( cM * MUE0 ) * 
                        TIME_SPAWN_FUNCTIONS::sinewave2d( 
                        T, 0.0, kx, curPos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
                        simbox.cellsize[0]), ky, curPos[1]-SPAWN_POS[1] ) *
                        TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
                        TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
                    /* E is half a time step later ! */
                    itm->E[Z] = TIME_SPAWN_FUNCTIONS::sinewave2d( 
                        T, 0.5*DELTA_T, kx, curPos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
                        simbox.cellsize[0]), ky, curPos[1]-SPAWN_POS[1] ) *
                        TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
                        TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
                    if ( std::abs(itm->E[Z]) < 0.001) {
                        itm->E[Z] = 0;
                        itm->H[Y] = 0;
                    }
                }
            }
            /********* Using Lens to make point source to plane wave **********/
            if ( WAVE_SPAWN_SETUP == 7 ) {
                /**************************************************************
                 *           y                                                *
                 * Absorber  ^    Air/Vacuum                                  *        
                 *         __|______         o...Source = M + R = (x0,y0)     *     
                 *       -   |  -   --                                        *
                 *      - Air|   -    \                                       *
                 *      |    o    | nG |                                      * 
                 *      -    |   -    /                                       *   
                 *       - __|__-___--                                        *
                 *           |                                                *
                 *           |---------------> x                              *
                 **************************************************************/
                 /* We want the center to be exactly the center there of the  *
                  * wave which will be spawned! That's why we search for a    *
                  * wanted pos and get back the center of the cell containing *
                  * that position                                             */
                /* static bool firstCall = true;
                 if ( firstCall ) {
                     firstCall = false;
                     R = 2*LAMBDA;
                     VecD guessedPos = VecD(ABSORBING_BORDER_THICKNESS + R, SIM_SIZE[1]/2 );
                     OctCell & t_simbox = *((OctCell*)(combox.tree.FindLeafContainingPos(guessedPos)->data[OctreeCommType::CELL_DATA_INDEX]));
                     M = t_simbox.getGlobalPosition( t_simbox.findCellContaining( guessedPos ) );
                     SPAWN_POS = M; // doesn't work with M, why !!! It's a border case, but anyway ...
                     tout << "R: " << R << ", M: " << M << ", M_again: " << t_simbox.getGlobalPosition( t_simbox.findCellContaining( M ) ) << "\n";
                 }*/
                 const double nVacuum = 1.0;
                 const double nLense  = 1.33;
                 const double e    = nVacuum / nLense;
                 const double & b  = 4*R;
                 const double a    = b / sqrt( 1 - e*e );
                 const double & x  = curPos[0];
                 const double & y  = curPos[1];
                 const double & x0 = M[0];
                 const double & y0 = M[1];
                 const double r    = (curPos-M).norm();
                 const double phi  = atan( (y-y0)/(x-x0) );
                 double xLeftCirc  = fabs(y-y0) < R ? x0 - R*sqrt( 1 - pow( (y-y0)/R, 2 ) ) : x0;
                 double xRightCirc = fabs(y-y0) < R ? x0 + R*sqrt( 1 - pow( (y-y0)/R, 2 ) ) : x0;
                 double xEllipseRight = fabs(y-y0) < b ? x0 + e*a + a*sqrt( 1 - pow( (y-y0)/b ,2 ) ) : x0;
                 double xEllipseLeft  = fabs(y-y0) < b ? x0 + e*a - a*sqrt( 1 - pow( (y-y0)/b ,2 ) ) : x0;
                 bool isAbsorber   = x < x0 and x < xLeftCirc;
                 bool isLense      = x >= xRightCirc and x <= xEllipseRight and x >= xEllipseLeft;
                 bool isSpawnGuard = r >= R and r <= R + 2*ABSORBING_BORDER_THICKNESS and fabs(phi) > M_PI/6;
                 /* comparison with NaN above will be always false! */
                 if ( isAbsorber or isSpawnGuard ) {
                    itm->sigmaE  = ABSORBER_STRENGTH;
                    itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
                 }
                 /*if ( isLense ) {
                    itm->epsilon = EPS0 * nLense*nLense;
                 }*/
                 /* Set AbsorberBorders */
                 if ( (curPos[X] - 0 < ABSORBING_BORDER_THICKNESS) or
                      (tree.center[X] + tree.size[X]/2 - curPos[X] < ABSORBING_BORDER_THICKNESS) or
                      (tree.center[Y] + tree.size[Y]/2 - curPos[Y] < ABSORBING_BORDER_THICKNESS) or
                      (curPos[Y] - 0 < ABSORBING_BORDER_THICKNESS)
                    )
                 {
                     itm->sigmaE = ABSORBER_STRENGTH;
                     itm->sigmaM = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
                 }
            }
            /******************************************************************/
        }
    }
    combox.PrintPNG( 0, "n_init", returnn );
    combox.PrintPNG( 0, "Ez_init", returnEz );
    combox.PrintPNG( 0, "Hy_init", returnHy );

    /**************************************************************************/
    /* (6) Actual Timestepping ************************************************/
    /**************************************************************************/
    MPI_Barrier(MPI_COMM_WORLD);
    double tStart = MPI_Wtime();
    double tLast  = tStart;
	for ( int timestep=0; timestep < NUMBER_OF_STEPS; ++timestep ) {
    for ( double internaltimestep = 0; internaltimestep < 1; internaltimestep +=
        1./pow( 2, combox.maxLevel - combox.minLevel ) )
    {
        #if DEBUG_MAIN_YEE >= 100
            if (timestep == 0) for (int lvl = combox.minLevel; lvl<=combox.maxLevel; lvl++)
                tout << "internaltimestep: " << internaltimestep << " mod " << 1/pow(2,lvl) << " = " << 1 / pow(2, combox.maxLevel - lvl ) << " == 0 ? " << (fmod( internaltimestep, 1 / pow(2, combox.maxLevel - lvl ) ) == 0) << "\n";
        #endif
        double t = timestep + internaltimestep;

        if ( WAVE_SPAWN_SETUP == 1 or WAVE_SPAWN_SETUP == 2 or WAVE_SPAWN_SETUP == 7 ) {
            /* Function Generator on Cell in the Center */
            OctreeType::Node * node = tree.FindLeafContainingPos( SPAWN_POS );
            if ( ((OctreeCommType::CommData*)node->data[OctreeCommType::COMM_HEADER_INDEX])->rank == combox.rank ) {
                OctCell & simbox = *((OctCell*)node->data[OctreeCommType::CELL_DATA_INDEX]);
                #if DEBUG_MAIN_YEE >= 100
                    tout << "Find position " << SPAWN_POS << " in tree sized " << tree.size << " positioned at center " << tree.center << " returned OctreeCell at " << node->center << " => Searching for position in simbox of that node with abspos " << simbox.abspos << ", localcells " << simbox.localcells << " and cellsize " << simbox.cellsize << "\n";
                #endif
                VecI targetIndex = simbox.findCellContaining( SPAWN_POS );
                simbox.t[0]->cells[targetIndex].E[Z] = 20*t_spawn_func( t * DELTA_T_SI );
                #if DEBUG_MAIN_YEE >=100
                    tout << "Write to source in cell " << targetIndex << " in node at " << node->center << "\n";
                #endif
            }
        }
        if ( WAVE_SPAWN_SETUP == 2 ) {
            /* shield function generator in one direction */
            OctreeType::Node * node = tree.FindLeafContainingPos( SPAWN_POS );
            if ( ((OctreeCommType::CommData*)node->data[OctreeCommType::COMM_HEADER_INDEX])->rank == combox.rank ) {
                OctCell & simbox = *((OctCell*)node->data[OctreeCommType::CELL_DATA_INDEX]);
                VecI targetIndex = simbox.findCellContaining( SPAWN_POS );
                /*simbox.t[0]->cells[targetIndex + VecI(8,-4)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,-3)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,-2)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,-1)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,0)   ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,1)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,2)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,3)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(8,4)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,-4)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,-3)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,-2)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,-1)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,0)   ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,1)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,2)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,3)  ].E[Z] = 0;
                simbox.t[0]->cells[targetIndex + VecI(9,4)  ].E[Z] = 0;*/
            }
        }
        if ( WAVE_SPAWN_SETUP == 3 ) {
            /**********************************************************************
             * Sine plane Wave going to Direction alpha and beginning line going  *
             * through pos0  y                                                    *
             *               ^                                                    *
             *               | \     e.g. p0 ( line includes p0! )                *
             *               |  --  /                                             *
             *               |    \     alpha                                     *
             *               |     --  /                                          *
             *               |_______\__________ x                                *
             **********************************************************************/
            //double const n = 1.0;
            double alpha   = 0 / 360. * 2*M_PI; //0.9*asin(1./n); // radian
            //VecI pos0(0); pos0[X] = 6*lambda;
            double T0x     = LAMBDA / SPEED_OF_LIGHT;
            double T0y     = LAMBDA / SPEED_OF_LIGHT;
            double kx      = 2*M_PI/LAMBDA * cos(alpha);
            double ky      = 2*M_PI/LAMBDA * sin(alpha);
            double width   = SPAWN_AREA_SIZE[1];
            VecD pos       = SPAWN_POS;// + SPAWN_AREA_SIZE - 1;
            assert( SPAWN_POS[1]+SPAWN_AREA_SIZE[1] < tree.center[1] + tree.size[1] );

            #if DEBUG_MAIN_YEE >= 100
            if ( timestep == 0 ) 
                tout << "\nInternal Timestep: " << internaltimestep;
            #endif
            /* while inside spawning pos and inside simulation area */
            while ( (pos[1] >= SPAWN_POS[1]) and (pos[1] < SPAWN_POS[1] + SPAWN_AREA_SIZE[1]) and
                     pos >= tree.center - tree.size/2 and pos <= tree.center + tree.size/2 )
            {
                OctreeType::Node * curnode = tree.FindLeafContainingPos( pos );
                if ( ((OctreeCommType::CommData*)curnode->data[OctreeCommType::COMM_HEADER_INDEX])->rank == combox.rank ) {
                    OctCell & simbox = *((OctCell*)curnode->data[OctreeCommType::CELL_DATA_INDEX]);
                    VecI targetIndex = simbox.findCellContaining( pos );
                    VecI index = targetIndex;
                    if ( timestep == 0 )
                        tout << "\n-------------------------------------------";
                    /*for ( index[0] = targetIndex[0]; index[0] - targetIndex[0] < SPAWN_AREA_SIZE[0]*pow(2,curnode->getLevel() - combox.minLevel) and simbox.inArea(index, SimulationBox::CORE + SimulationBox::BORDER ); index[0]++ )
                        simbox.t[0]->cells[index].E[Z] = 0;*/
                    for ( pos[1] = pos[1]; simbox.inArea(index, SimulationBox::CORE + SimulationBox::BORDER ) and
                          (pos[1] < SPAWN_POS[1] + SPAWN_AREA_SIZE[1]); pos[1] += simbox.cellsize[1] ) {
                        index = simbox.findCellContaining( pos );
                        /*for ( index[0] = targetIndex[0]; index[0] - targetIndex[0] < pow(2,curnode->getLevel() - combox.minLevel) and simbox.inArea(index, SimulationBox::CORE + SimulationBox::BORDER ); index[0]++ )
                           simbox.t[0]->cells[index].E[Z] = 0;*/
                        pos[0] = SPAWN_POS[0];
                        for ( index[0] = targetIndex[0]; index[0] - targetIndex[0] < SPAWN_AREA_SIZE[0] and simbox.inArea(index, SimulationBox::CORE + SimulationBox::BORDER ); index[0]++ ) {
                           simbox.t[0]->cells[index].E[Z] = TIME_SPAWN_FUNCTIONS::sinewave2d( T0x, t * DELTA_T, kx, pos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] - simbox.cellsize[0]),
                           ky, pos[1]-SPAWN_POS[1] ) /* * TIME_SPAWN_FUNCTIONS::PSQ_STEP( T0y, t * DELTA_T ) *//* * 50* TIME_SPAWN_FUNCTIONS::gauss( pos[1], SPAWN_POS[1] + SPAWN_AREA_SIZE[1]/2, SPAWN_AREA_SIZE[1]/4 )*/;
                            #if DEBUG_MAIN_YEE >= 100
                            if ( timestep == 0 ) {
                                if ( pos[0] == SPAWN_POS[0] )
                                    tout << "\n";
                                tout << index << " -> " << simbox.t[0]->cells[index].E[Z] << " | ";
                            }
                            #endif
                            pos[0] += simbox.cellsize[0];
                        }
                    }
                    /*if ( pos[1] > SPAWN_POS[1] + SPAWN_AREA_SIZE[1] ) {
                        index[1]--;
                        for ( index[0] = targetIndex[0]; index[0] - targetIndex[0] < 3*pow(2,curnode->getLevel() - combox.minLevel) and simbox.inArea(index, SimulationBox::CORE + SimulationBox::BORDER ); index[0]++ )
                            simbox.t[0]->cells[index].E[Z] = 0;
                        index[1]++;
                    }*/
                }
                pos[0] = SPAWN_POS[0];
                pos[1] = tree.toGlobalCoords( curnode->center + curnode->size/2 + pow( 0.5, combox.maxLevel ) / VecD(cellsPerOctreeCell) / 2. )[1];
                if ( pos > tree.center - tree.size/2 and pos < tree.center + tree.size/2 )
                    assert( curnode != tree.FindLeafContainingPos( pos ) );
            }
            #if 1==0
            //std::cout << "Spawning slanted sine wave: \n";
            for (int j=0; j<2; j++) {
                VecI pos( pos0+GUARDSIZE ); pos[Y]+=j;
                for (int i=0; i<10*LAMBDA; i++) {
                    pos[X]++;
                    data[pos].t[0].E[Z] = TIME_SPAWN_FUNCTIONS::sinewave2d( T0x, t * DELTA_T, kx, i*CELL_SIZE_X, ky, j*CELL_SIZE_Y ) * TIME_SPAWN_FUNCTIONS::PSQ_STEP( T0y, t * DELTA_T );
                        //* TIME_SPAWN_FUNCTIONS::sinewave( T0y, t * DELTA_T );
                    //std::cout << data[pos].t[0].E[Z] << " ";
                }
                //std::cout << std::endl;
            }
            #endif
        }

        /* 0 is timestep to be calculated, 1 is the previous one */
        combox.StartGuardUpdate( 0 );

        /* Traverse all Octree Nodes and apply Yee-Solver there ( Do this     *
         * with a combox iterator ???                                         */
        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( combox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[combox.CELL_DATA_INDEX]);
            double dt = DELTA_T / pow( 2, it->getLevel() - combox.minLevel );
            if ( fmod( internaltimestep, dt ) == 0.0 )
                YeeSolver::CalcH( simBox, 0, 0, SimulationBox::CORE, dt );
        }

        combox.FinishGuardUpdate( 0 );

        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( combox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[combox.CELL_DATA_INDEX]);
            double dt = DELTA_T / pow( 2, it->getLevel() - combox.minLevel );
            if ( fmod( internaltimestep, dt ) == 0.0 )
                YeeSolver::CalcH( simBox, 0, 0, SimulationBox::BORDER, dt );
        }

        /* In the next halfstep, do the same for E-Field, but stay in timestep*/
        combox.StartGuardUpdate( 0 );

        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( combox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[combox.CELL_DATA_INDEX]);
            double dt = DELTA_T / pow( 2, it->getLevel() - combox.minLevel );
            if ( fmod( internaltimestep, dt ) == 0.0 )
                YeeSolver::CalcE( simBox, 0, 0, SimulationBox::CORE, dt );
        }

        combox.FinishGuardUpdate( 0 );

        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( combox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[combox.CELL_DATA_INDEX]);
            double dt = DELTA_T / pow( 2, it->getLevel() - combox.minLevel );
            if ( fmod( internaltimestep,dt ) == 0.0 )
                YeeSolver::CalcE( simBox, 0, 0, SimulationBox::BORDER, dt );
        }
        MPI_Barrier(MPI_COMM_WORLD);
    } // internaltimestep

		if (timestep+1 % PNG_INTERVAL == 0) {
			static int framecounter = 0;
			framecounter++;
			char filename[100];
			/*sprintf( filename, "Ex_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnEx );
            sprintf( filename, "Ey_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnEy );*/
            sprintf( filename, "Ez_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnEz );
			/*sprintf( filename, "Hx_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnHx );*/
			sprintf( filename, "Hy_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnHy );
			/*sprintf( filename, "Hz_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnHz );*/
			/*sprintf( filename, "All_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnEandH );*/
            /*sprintf( filename, "n_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnn );*/
		}

        MPI_Barrier(MPI_COMM_WORLD);
        #if DEBUG_MAIN_YEE >= 90
            tout << "Timestep " << timestep << " took " << MPI_Wtime() - tLast << " seconds\n" << std::flush;
        #endif
        tLast = MPI_Wtime();
	}

    tout << "All " << NUMBER_OF_STEPS << " together took " << tLast - tStart << " seconds\n";
    tout << "The whole program took " << tLast - tProgramStart << " seconds\n";

    if ( combox.maxLevel == 0 ) {
        typename OctreeType::iterator it=tree.begin();
        OctCell & simbox = *((OctCell*)it->data[combox.CELL_DATA_INDEX]);
        tout << "\nE.z: {";
        for ( typename OctCell::IteratorType itm = simbox.getIterator( 0,
              SimulationBox::CORE + SimulationBox::BORDER ).begin();
              itm != itm.end(); ++itm )
        if ( itm.icell[1] == 10 )
            tout << itm->E[Z] << ", ";
        tout << "}\nH.y: {";
        for ( typename OctCell::IteratorType itm = simbox.getIterator( 0,
              SimulationBox::CORE + SimulationBox::BORDER ).begin();
              itm != itm.end(); ++itm )
        if ( itm.icell[1] == 10 )
            tout << itm->H[Y] << ", ";
        tout << "}\n";
    }
    
    /* doesn't work in destructor, because it would be called too late */
    MPI_Finalize();
    return 0;
}
