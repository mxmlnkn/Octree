#pragma once

#include <cstdlib>   // malloc
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <assert.h>
#include <list>
#include <algorithm> // find
#include <mpi.h>
#include <pngwriter.h>
#include "math/TVector.h"
#include "octree/Octree.h"
#include "simbox/SimulationBox.h"
#include "Directions.h"
#include "teestream/TeeStream.h"
#include "math/TVector.h"
#include "math/TBaseMatrix.h"

#ifndef DEBUG_COMMUNICATOR
    #define DEBUG_COMMUNICATOR 1
#endif

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
class OctreeCommunicator{
public:
    /* Make Copy- and Assignment-Constructors unavailable, to prevent       *
     * copies. If there is an error at                                      *
     *   Singleton var = Singleton::getInstance                             *
     * use a reference!                                                     *
     *   Singleton & var = Singleton::getInstance                           */
    OctreeCommunicator(const OctreeCommunicator&) { assert(false); };
    OctreeCommunicator & operator=(const OctreeCommunicator & ) { assert(false); };

    typedef Vec<int,T_DIM> VecI;
    typedef Vec<double,T_DIM> VecD;
    typedef BaseMatrix<T_CELLTYPE,T_DIM> CellMatrix;
    typedef SimulationBox::SimulationBox<T_DIM,T_CELLTYPE> OctCell;

    struct CommData {
        int rank;
        double weighting;
    };

    int  rank, worldsize;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int  processor_name_length;

    T_OCTREE & tree;
    VecI periodic;
    int minLevel, maxLevel;

    static const int COMM_HEADER_INDEX;
    static const int CELL_DATA_INDEX;
    CommData* comDataPtr;
    int distOrdering; // ordering used when distributing cells. Can be used to traverse effectively all own cells

    int NLeaves, NOwnLeaves;

    VecI cellsPerOctreeCell;
    int guardsize, timestepbuffers;

    /* stores pointer to node + which direction / side to send or recv */
    struct BorderData {
        typename T_OCTREE::Node * node;
        int direction;
        bool operator==( const struct BorderData & rhs ) const {
            return ( this->node == rhs.node and this->direction == rhs.direction );
        }
        /* Strict weak Ordering to sort send and receive buffers */
        bool operator<(const struct BorderData & b) const {
            /* return true, if a comes before b */
            for (int i=0; i<T_DIM; i++) {
                if (this->node->center[i] < b.node->center[i]) return true;
                if (this->node->center[i] > b.node->center[i]) return false;
            }
            if ( this->direction < b.direction ) return true;
            if ( this->direction > b.direction ) return false;
            return false; // if they are equal
        }
    };
    typedef std::list<struct BorderData> ToCommList;

    struct NeighborData {
        /* this has nsendmatrices * cellsPerOctreeCell*guardsize elements! */
        T_CELLTYPE * sendData;
        int cellsToSend;
        std::list<struct BorderData> sSendData; // s->source

        /* this has nrecvmatrices * cellsPerOctreeCell*guardsize elements!    *
         * Meaning it stores all the borders to send in an order both threads *
         * agree on. recvData is basically packed data                        */
        T_CELLTYPE * recvData;
        int cellsToRecv;
        std::list<struct BorderData> sRecvData; // s->source

        /* simple constructor / initializer */
        NeighborData(void) : sendData(NULL), cellsToSend(0), sSendData(), 
                             recvData(NULL), cellsToRecv(0), sRecvData() {};
        NeighborData & operator=( const NeighborData & src ) { assert(false); }
        NeighborData( const NeighborData & src ) { assert(false); }
    };

    /* worldsize elements, meaning rank includes itself. sRecvData entries in *
     * neighbors[rank] are meant to be received from rank by us.              *
     * sSendData entries in neighbors[rank] are meant to be sent to rank      */
    NeighborData * neighbors;
    MPI_Request * sendrequests;
    MPI_Request * recvrequests;
    int timestepSent;

    /******************************** Methods *********************************/
    
    /* squash the cube (CORE+BORDER) we want to send to guardsize in each     *
     * direction we want to sent it to, e.g.                                  *
     *   guard = 1, core+border = 3 => begin with 3x3x3 cube.                 *
     *   send to LEFT or RIGHT => make it a 1x3x3 surface                     *
     *   send to TOP or BOTTOM => make it a 3x1x3 surface                     *
     *   send to LEFT+TOP      => make it a 1x1x3 line of cells (edge)        */
    VecI getBorderSizeInDirection( const int & direction ) const;
    VecI getBorderPositionInDirection( const int & direction ) const;
    VecI getGuardSizeInDirection( const int & direction ) const;
    VecI getGuardPositionInDirection( const int & direction ) const;

    /* Start asynchronous sends of all Border Data belonging to other threads *
     * If neighbor-cell owned by this thread itself, interpolate and copy     *
     * data to guards.                                                        */
    void StartGuardUpdate( int timestep );

    /* Waits for MPI_Recv_Requests to finish and then copies the received     *
     * partial matrix from the buffer to the simulation matrix. Because of    *
     * complicated stride in all directions it can't be communicated directly */
    void FinishGuardUpdate( int timestep );

#if DEBUG_COMMUNICATOR >= 100
    void Print( void ) {
        terr << "Printing from Communicator on rank " << this->rank
             << " with cartesian coordinates: " << VecI( this->coords )
             << " My neighbors are: " << endl;

        for (int direction = 1; direction <= this->nneighbors; direction++) {
            terr << (VecI(this->coords) + getDirectionVector<T_DIM>( direction ))
                 << "(lin=" << direction << ") -> Rank: "
                 << this->neighbors[direction] << endl << flush;
        }

        /* Wait a bit till everything is flushed out, to not get scrambled output */
        double t0 = MPI_Wtime();
        while( MPI_Wtime() - t0 < 0.1 ) {}
    }
#endif

    OctreeCommunicator(void);
    OctreeCommunicator(T_OCTREE& tree, const VecI periodic);
    ~OctreeCommunicator();

    /* this allocates memory for commdata, which will be hooked into every    *
     * octree cell. After this function therefore GrowUp shouldn't and can't  *
     * be called anymore, because GrowUp assumes empty Octree cells. This     *
     * therefore 'finalizes' the octree. This also sets the weighting         *
     * todo: take argument for weighting function/operator                    */
    void initCommData( VecI cellsPerOctreeCell, int guardsize, int timesteps );
    /* initCommData needs to be called before this ! */
    void distributeCells( int ordering = Octree::Ordering::Hilbert );

    /* The function must take YeeCell as an argument and should return a      *
     * 3-dimensional double Vector VecD (R,G,B)                               */
    template<typename T_FUNC>
    void PrintPNG( int timestep, const char * name, T_FUNC function, bool timestamp = true );
};

#include "Communicator.tpp"
