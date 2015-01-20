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
#include <random>   // normal_distribution
#include <pngwriter.h>
#include <list>
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


int main( void /*int argc, char **argv*/ )
{
    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    typedef Vec<double,SIMDIM> VecD;
    typedef Vec<int   ,SIMDIM> VecI;
    
    double tProgramStart = MPI_Wtime();

    /* Call (indirectly) Basic Communicator-, Octree- and File-Constructors */
    VecD cellSize(CELL_SIZE);
    VecD NCells(NUMBER_OF_CELLS);
    VecD globSize(cellSize*NCells), globCenter(0.5*globSize);
    typedef typename Octree::Octree<SIMDIM> OctreeType;
    OctreeType tree( globCenter, globSize );
    //typedef BaseMatrix<YeeCell,SIMDIM> CellMatrix;
    typedef OctreeCommunicator<SIMDIM,OctreeType,YeeCell> OctreeCommType;
    OctreeCommType comBox(tree,VecI(true));
    typedef SimulationBox::SimulationBox<SIMDIM,YeeCell> OctCell;
    if ( comBox.rank != 0) {
        std::cerr << "Disable tout verbosity on rank " << comBox.rank << "\n";
        tout.verbosity = 0;
    }
    tout.Open("out",comBox.rank);
    terr.Open("err",comBox.rank);
    /* Initialize SVG output file */
    std::stringstream svgfn;
    svgfn << "Octree_worldsize-" << comBox.worldsize << "_rank-" << comBox.rank;
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
    tout << "\n";

    /**************************************************************************/
    /* (1) Setup Octree Refinement ********************************************/
    /**************************************************************************/
#define INITSETUP 6
#if INITSETUP == 5 or INITSETUP == 6
    /********* refine all cells to initial homogenous min-Refinement **********/
    for ( int lvl=0; lvl<INITIAL_OCTREE_REFINEMENT; lvl++) {
        for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it )
            if ( it->IsLeaf() and it->getLevel()==lvl ) it->GrowUp();
    }
#endif
#if INITSETUP == 6
    /*********************** Refine certain boundaries ************************/
    assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
    VecD M(0.5*tree.size);          // center of circle
    double R = 0.4*tree.size.min(); // radius of circle
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
#if INITSETUP == 7
    tree.root->GrowUp();
    tree.FindLeafContainingPos( 0.75*globSize )->GrowUp();
#endif

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
    tout << "Initial Refinement: " << INITIAL_OCTREE_REFINEMENT << "\n"
         << "Minimum refinement Level found: " << tree.getMinLevel() << "\n"
         << "Maximum refinement Level found: " << tree.getMaxLevel() << "\n"
         << "Cells per Octree: " << cellsPerOctreeCell << "\n\n";
    comBox.initCommData( cellsPerOctreeCell, GUARDSIZE, 3 /*timestepbuffer*/ );

    const int ORDERING = Octree::Ordering::Hilbert; // -> Parameters -> need to include Octree there :S ?
    comBox.distributeCells( ORDERING );

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
    if ( it->IsLeaf() ) { if ( comBox.rank == ((OctreeCommType::CommData*)it->
         data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
    {
        OctCell & data = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
        BaseMatrix<YeeCell,SIMDIM> & matrix = data.t[0]->cells;
        typename OctCell::IteratorType itm = data.getIterator( 0, SimulationBox::CORE );
        for ( itm = itm.begin(); itm != itm.end(); ++itm ) {
            matrix[itm.icell].E[X] = 1.0;
        }
        itm = data.getIterator( 0, SimulationBox::BORDER );
        for ( itm = itm.begin(); itm != itm.end(); ++itm ) {
            matrix[itm.icell].E[X] = double(2*rand()-1.0)/double(RAND_MAX);
            matrix[itm.icell].H[X] = double(2*rand()-1.0)/double(RAND_MAX);
        }
    } else if ( it->data.size() > comBox.CELL_DATA_INDEX ) {
        assert(false);
        OctCell & data = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
        BaseMatrix<YeeCell,SIMDIM> & matrix = data.t[0]->cells;
        typename OctCell::IteratorType itm = data.getIterator( 0, SimulationBox::BORDER );
        for ( itm = itm.begin(); itm != itm.end(); ++itm ) {
            matrix[itm.icell].E[X] = +1.0 * ( 1+ comBox.rank ) / (comBox.worldsize+1);
        }
    }}

    comBox.PrintPNG( 0, "TestGuardCommunication_a", returnEandH );

    tout << "Beginning asynchronous communication operations\n";
    comBox.StartGuardUpdate(0);
    tout << "Finishing up asynchronous communication operations\n";
    comBox.FinishGuardUpdate(1);

    comBox.PrintPNG( 1, "TestGuardCommunication_b", returnEandH );

    /* Copy Guard to Border by adding up all neighbors (only B-field) */
    for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
    if ( it->IsLeaf() ) if ( it->data.size() > 1 )
    {
        OctCell & data = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
        typename OctCell::IteratorType itm = data.getIterator( 1 /* timestep */, SimulationBox::BORDER );
        for ( itm = itm.begin(); itm != itm.end(); ++itm )
        {
            itm->H[X] = 0.0;
            itm->E[X] = 0.0;
            int dirs[ compileTime::pow(2,SIMDIM) ] = {RIGHT,LEFT,TOP,BOTTOM};
            int nNeighboringGuards = 0;
            for (int i=0; i < compileTime::pow(2,SIMDIM); ++i)
            {
                VecI neighbor = itm.icell + getDirectionVector<SIMDIM>(dirs[i]);
                if ( data.inArea( neighbor, SimulationBox::GUARD ) )
                {
                    nNeighboringGuards++;
                    itm->H[X] += data.t[1]->cells[neighbor].H[X];
                    itm->E[X] += data.t[1]->cells[neighbor].E[X];
                }
            }
            itm->H[X] /= nNeighboringGuards;
            itm->E[X] /= nNeighboringGuards;
        }
    }
    comBox.PrintPNG( 1, "TestGuardCommunication_c", returnEandH );

    /**************************************************************************/
    /* (5) Setup Simulation Start Data like Material **************************/
    /**************************************************************************/
    for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
    if ( it->IsLeaf() ) if ( comBox.rank == ((OctreeCommType::CommData*)it->
         data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
    {
        OctCell & data = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);

        
        for ( typename OctCell::IteratorType itm = data.getIterator( 0,
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

        for ( typename OctCell::IteratorType itm = data.getIterator( 0,
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
            VecD curPos = it->center + (VecD(itm.icell) - itm.guardsize + 0.5) /
                          VecD(cellsPerOctreeCell) * it->size;
            curPos = tree.toGlobalCoords( curPos );
            assert( itm.ncells - 2*itm.guardsize == cellsPerOctreeCell );

            #define YEE_INIT_SETUP 4
            #if YEE_INIT_SETUP == 1
                /* Result 009: absorbing Material on right side */
                if ( curPos[0] > 0.2 ) {
                    /* For INF instead of 2e8 it diverges! -> characteristic *
                     * wave length for exponential decay -> 0 for INF ?       */
                    itm->sigmaE  = 2e8;
                    itm->sigmaM  = 2e8 * MUE0/EPS0;
                }
            #endif
            #if YEE_INIT_SETUP == 2
                /* Result: 014 - broken total reflexion (two glass plates     *
                 *               with small vacuum/air slit inbetween         */
                const double n  = 1.33; // = sqrt( eps_r * mue_r )
                if ( curPos[X] < 0.333*tree.size[X]
                or   curPos[X] > 0.350*tree.size[X]  )
                    itm->epsilon = EPS0 * n*n;
            #endif
            #if YEE_INIT_SETUP == 3
                /* Spawn Barrier with one slit and perfectly reflecting material else */
                const double wy = LAMBDA;
                if ( curPos[X] > LAMBDA and curPos[X] < 2*LAMBDA )
                if ( curPos[Y] > 0.5 - wy/2 and curPos[Y] < 0.5 + wy/2 ) {
                        itm->epsilon  = INF;//2*EPS0;
                        //itm->mu       = INF;
                }
            #endif
            #if YEE_INIT_SETUP == 4 and INIT_SETUP == 3
                const double n  = 1.33; // = sqrt( eps_r * mue_r )
                if ( (curPos - M).norm() < R )
                    itm->epsilon = EPS0 * n*n;
            #endif
        }
    }
    comBox.PrintPNG( 0, "output/n_init", returnn, false );

    /**************************************************************************/
    /* (6) Actual Timestepping ************************************************/
    /**************************************************************************/
    double tStart = MPI_Wtime();
    double tLast  = tStart;
	for ( int timestep=0; timestep < NUMBER_OF_TIMESTEPS; ++timestep )
	{
		#define TIME_SPAWN_SETUP 1
        #if TIME_SPAWN_SETUP == 1
            /* Function Generator on Cell in the Center */
            OctreeType::Node * node = tree.FindLeafContainingPos( tree.center + 0.11*tree.size );
            if ( ((OctreeCommType::CommData*)node->data[OctreeCommType::COMM_HEADER_INDEX])->rank == comBox.rank ) {
                OctCell & cellMatrix = *((OctCell*)node->data[OctreeCommType::CELL_DATA_INDEX]);
                VecI targetIndex = cellMatrix.findCellContaining( tree.center + 0.11*tree.size );
                cellMatrix.t[0]->cells[targetIndex].E[Z] = t_spawn_func( timestep * DELTA_T_SI );
                #if DEBUG_MAIN_YEE >= 100
                    tout << "Write to source in cell " << targetIndex << " in node at " << node->center << "\n";
                #endif
            }
        #endif
        #if TIME_SPAWN_SETUP == 2
            /* Function Generator on left side creates sine wave */
            /*VecI pos(GUARDSIZE);
            for ( int y = 0; y < N_CELLS_Y - 2*GUARDSIZE; y++ ) {
                pos[Y] = GUARDSIZE + y;
                //if( (pos[Y]>20 and pos[Y]<N_CELLS_Y/2-20) or (pos[Y]>N_CELLS_Y/2+20 and pos[Y]<N_CELLS_Y-20)
                    data[pos].E[tcur][Z] = t_spawn_func( timestep * DELTA_T_SI );
                //if( (pos[Y]>20 and pos[Y]<N_CELLS_Y/2-20) or (pos[Y]>N_CELLS_Y/2+20 and pos[Y]<N_CELLS_Y-20)
                //	data[pos].E[tcur][Z] = t_spawn_func( timestep * DELTA_T_SI );
            }*/
        #endif
        #if TIME_SPAWN_SETUP == 3
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
            double const n = 1.33;
            double alpha   = 45. / 360. * 2*M_PI; //0.9*asin(1./n); // radian
            double lambda  = 10e-9 / UNIT_LENGTH;
            VecI pos0(0); pos0[X] = 6*lambda;
            double T0x     = lambda / (SPEED_OF_LIGHT/n);
            double T0y     = lambda / (SPEED_OF_LIGHT/n);
            double kx      = 2*M_PI/lambda * cos(alpha);
            double ky      = 2*M_PI/lambda * sin(alpha);
            //std::cout << "Spawning slanted sine wave: \n";
            for (int j=0; j<2; j++) {
                VecI pos( pos0+GUARDSIZE ); pos[Y]+=j;
                for (int i=0; i<10*lambda; i++) {
                    pos[X]++;
                    data[pos].t[0].E[Z] = TIME_SPAWN_FUNCTIONS::sinewave2d( T0x, timestep * DELTA_T, kx, i*CELL_SIZE_X, ky, j*CELL_SIZE_Y ) * TIME_SPAWN_FUNCTIONS::PSQ_STEP( T0y, timestep * DELTA_T );
                        //* TIME_SPAWN_FUNCTIONS::sinewave( T0y, timestep * DELTA_T );
                    //std::cout << data[pos].t[0].E[Z] << " ";
                }
                //std::cout << std::endl;
            }
        #endif

        /* Swap timestep buffers for alle SimulationBoxes before Calculating */
        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( comBox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
            simBox.copyCurrentToPriorTimestep();
        }

        /* 0 is timestep to be calculated, 1 is the previous one */
        comBox.StartGuardUpdate( 1 );

        /* Traverse all Octree Nodes and apply Yee-Solver there ( Do this     *
         * with a comBox iterator ???                                         */
        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( comBox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
            YeeSolver::CalcH( simBox, 0, 1, SimulationBox::CORE );
        }

        comBox.FinishGuardUpdate( 1 );

        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( comBox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
             YeeSolver::CalcH( simBox, 0, 1, SimulationBox::BORDER );
        }

        /* In the next halfstep, do the same for E-Field, but stay in timestep*/
        comBox.StartGuardUpdate( 0 );

        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( comBox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
             YeeSolver::CalcE( simBox, 0, 1, SimulationBox::CORE );
        }

        comBox.FinishGuardUpdate( 0 );

        for ( typename OctreeType::iterator it=tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) if ( comBox.rank == ((OctreeCommType::CommData*)it->
             data[OctreeCommType::COMM_HEADER_INDEX])->rank ) /* owned cells */
        {
            OctCell & simBox = *((OctCell*)it->data[comBox.CELL_DATA_INDEX]);
            YeeSolver::CalcE( simBox, 0, 1, SimulationBox::BORDER );
        }

        #ifndef NPNGOUTPUT
		if (timestep % 1 == 0) {
			static int framecounter = 0;
			framecounter++;
			char filename[100];
			sprintf( filename, "output/Ex_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnEx, false );
            sprintf( filename, "output/Ey_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnEy, false );
            sprintf( filename, "output/Ez_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnEz, false );
			sprintf( filename, "output/Hx_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnHx, false );
			sprintf( filename, "output/Hy_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnHy, false );
			sprintf( filename, "output/Hz_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnHz, false );
			sprintf( filename, "output/All_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnEandH, false );
            sprintf( filename, "output/n_%05i", framecounter );
            comBox.PrintPNG( 0, filename, returnn, false );
		}
        #endif

        MPI_Barrier(MPI_COMM_WORLD);
        #if DEBUG_MAIN_YEE >= 90
            tout << "Timestep " << timestep << " took " << MPI_Wtime() - tLast << " seconds\n";
        #endif
        tLast = MPI_Wtime();
	}

    tout << "All " << NUMBER_OF_TIMESTEPS << " together took " << tLast - tStart << " seconds\n";
    tout << "The whole program took " << tLast - tProgramStart << " seconds\n";

    /* doesn't work in destructor, because it would be called too late */
    MPI_Finalize();
    return 0;
}
