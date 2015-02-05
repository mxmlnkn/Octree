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

int argc;
char ** argv;

#include <iostream>
#include <cmath>    // sin
#include <cfloat>   // FLT_EPSILON
#include <cstdlib>  // malloc, srand, rand, RAND_MAX
#include "getopt.h"
#include <pngwriter.h>
#include <list>
#include <algorithm>
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


int main( int pargc, char **pargv )
{
    argc = pargc;
    argv = pargv;
    double tProgramStart = MPI_Wtime();

    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    typedef Vec<double,SIMDIM> VecD;
    typedef Vec<int   ,SIMDIM> VecI;

    #include "InitArguments.cpp"

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

    /* Call (indirectly) Basic Communicator-, Octree- and File-Constructors */
    VecD cellSize(CELL_SIZE);
    VecD NCells( NUMBER_OF_CELLS );
    VecD globSize(cellSize*NCells), globCenter(0.5*globSize);
    typedef typename Octree::Octree<SIMDIM> OctreeType;
    OctreeType tree( globCenter, globSize );
    //typedef BaseMatrix<YeeCell,SIMDIM> CellMatrix;
    typedef OctreeCommunicator<SIMDIM,OctreeType,YeeCell> OctreeCommType;
    OctreeCommType combox(tree,VecI(true));
    typedef typename OctreeCommType::CellIterator CellIterator;
    typedef SimulationBox::SimulationBox<SIMDIM,YeeCell> OctCell;
    if ( combox.rank != 0) {
        std::cerr << "Disable tout verbosity on rank " << combox.rank << "\n";
        tout.verbosity = 0;
    }
    tout.Open( basefilename.str() + std::string("out.txt"), combox.rank );
    terr.Open( basefilename.str() + std::string("err.txt"), combox.rank );

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
    tout << "PNG_INTERVAL             : " << PNG_INTERVAL               << "\n";
    tout << "combox.worldsize         : " << combox.worldsize           << "\n";
    tout << "combox.rank              : " << combox.rank                << "\n";
    tout << "\n";

    #include "InitOctreeRefinement.cpp"

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
    combox.distributeCells(ORDERING);

    /*********************** Initialize SVG output file ***********************/
    std::stringstream svgfilename;
    svgfilename << basefilename.str() << "Octree_worldsize-"
                << combox.worldsize << "_rank-" << combox.rank << ".svg";
    Octree::OctreeToSvg<SIMDIM> svgoutput( tree, svgfilename.str() );

    svgoutput.PrintGrid();

    /* Graphical output of cell-rank-mapping */
    for ( typename OctreeType::iterator it = tree.begin(); it != it.end(); ++it)
    if ( it->IsLeaf() ) {
        assert( it->data.size() > 0 );
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

    svgoutput.PrintTraversal(ORDERING);
    svgoutput.close();

    /* Set Core, Border and Guard to different test values */
    tout << "Set Borders..." << std::flush;
    for ( CellIterator it = combox.getCellIterator( SimulationBox::BORDER, 0 ); it != it.end(); ++it )
            it->E[0] = -1.0;
    tout << "OK\n" << std::flush;

    tout << "Write marked Border to png..." << std::flush;
    combox.PrintPNG( 0, basefilename.str() + std::string("Border_Grid"), returnEandH );
    tout << "OK\n" << std::flush;

    tout << "Setup Initial Simulation Data..." << std::flush;
    /**************************************************************************/
    /* (5) Setup Simulation Start Data like Material **************************/
    /**************************************************************************/
    for ( CellIterator it = combox.getCellIterator( SimulationBox::GUARD, 0 ); it != it.end(); ++it )
    {
        /* it traverses octree nodes, while itm traverses matrix cells! */
        it->E       = 0;
        it->H       = 0;
        it->epsilon = EPS0;
        it->mu      = MUE0;
        it->sigmaE  = 0;
        it->sigmaM  = 0;
    }

    for ( CellIterator itm = combox.getCellIterator( SimulationBox::CORE +
    SimulationBox::BORDER, 0 ); itm != itm.end(); ++itm )
    {
        itm->E       = 0;
        itm->H       = 0;
        itm->epsilon = EPS0;
        itm->mu      = MUE0;
        itm->sigmaE  = 0;
        itm->sigmaM  = 0;

        VecD curpos = itm.getGlobalPosition();

        /* needs ( write as function ? ):                                 *
         *  - itm: must be deferencable to YeeCell                        *
         *  - curpos: holds global position of center of cell (no         *
         *            internal units)                                     *
         *  - tree: Octree, especially for size and other attributes      *
         *  - Parameters.cpp                                              */
        #include "InitSimulationSetup.cpp"
        #include "InitWaveSpawn.cpp"
    }
    tout << "OK\n" << std::flush;

    tout << "Write initial n,Ez,Hy to png..." << std::flush;
    combox.PrintPNG( 0, basefilename.str() + std::string("n_init" ), returnn  );
    combox.PrintPNG( 0, basefilename.str() + std::string("Ez_init"), returnEz );
    combox.PrintPNG( 0, basefilename.str() + std::string("Hy_init"), returnHy );
    tout << "OK\n" << std::flush;

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
        double t = timestep + internaltimestep;
        #if DEBUG_MAIN_YEE >= 90
            tout << "t = " << t << ", " << std::flush;
        #endif

        #include "UpdateWaveSpawn.cpp"

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

		if ( (timestep+1) % PNG_INTERVAL == 0) {
			static int framecounter = 0;
			framecounter++;
			char filename[100];
			/*sprintf( filename, "Ex_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnEx );
            sprintf( filename, "Ey_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnEy );*/
            sprintf( filename, "Ez_%05i", framecounter );
            combox.PrintPNG( 0, basefilename.str() + std::string(filename), returnEz );
			/*sprintf( filename, "Hx_%05i", framecounter );
            combox.PrintPNG( 0, filename, returnHx );*/
			sprintf( filename, "Hy_%05i", framecounter );
            combox.PrintPNG( 0, basefilename.str() + std::string(filename), returnHy );
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
