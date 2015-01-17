#pragma once

#include <cstdlib>   // malloc
#include <iostream>
#include <list>
#include <algorithm> // find
#include <mpi.h>
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
private:
      /* Make Copy- and Assignment-Constructors unavailable, to prevent       *
       * copies. If there is an error at                                      *
       *   Singleton var = Singleton::getInstance                             *
       * use a reference!                                                     *
       *   Singleton & var = Singleton::getInstance                           */
      OctreeCommunicator(const OctreeCommunicator&);                // Don't Implement
      OctreeCommunicator& operator=(const OctreeCommunicator&);     // Don't implement

public:
    typedef Vec<int,T_DIM> VecI;
    typedef BaseMatrix<T_CELLTYPE,T_DIM> CellMatrix;
    typedef SimulationBox::SimulationBox<T_DIM,T_CELLTYPE> OctCell;

    struct CommData {
        int rank;
        int weighting;
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
    double totalCosts, optimalCosts;

    VecI cellsPerOctreeCell;
    int guardsize, timestepbuffers;
    
    /* stores pointer to node + which direction / side to send or recv */
    struct BorderData {
        typename T_OCTREE::Node * node;
        VecI direction;
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
            if (getLinearDirection(this->direction) < getLinearDirection(b.direction)) return true;
            if (getLinearDirection(this->direction) > getLinearDirection(b.direction)) return false;
            return false; // if they are equal
        }
    };
    typedef std::list<struct BorderData> ToCommList;
    
    struct NeighborData {
        /* this has nrecvmatrices * cellsPerOctreeCell*guardsize elements!    *
         * Meaning it stores all the borders to send in an order both threads *
         * agree on. recvmats is basically packed data                        */
        double recvmats[];
        std::list<struct BorderData> srecvmats; // s->source
        
        /* this has nsendmatrices * cellsPerOctreeCell*guardsize elements! */
        double sendmats[];
        std::list<struct BorderData> ssendmats; // s->source
    };
    
    /* worldsize elements, meaning rank includes itself. srecvmats entries in *
     * neighbors[rank] are meant to be received from rank by us.              *
     * ssendmats entries in neighbors[rank] are meant to be sent to rank      */
    NeighborData * neighbors;
    MPI_Request * recvrequests;
    MPI_Request * sendrequests;


#if 1==0
    VecI getBorderSizeInDirection( const int & direction ) const {
        VecI size  = simbox.localcells;
        /* squash the cube (CORE+BORDER) we want to send to guardsize in each *
         * direction we want to sent it to, e.g.                              *
         *   guard = 1, core+border = 3 => begin with 3x3x3 cube.             *
         *   send to LEFT or RIGHT => make it a 1x3x3 surface                 *
         *   send to TOP or BOTTOM => make it a 3x1x3 surface                 *
         *   send to LEFT+TOP      => make it a 1x1x3 line of cells (edge)    */
        VecI v( getDirectionVector<T_DIM>( direction ) );
        for (int axis=0; axis<T_DIM; ++axis)
            if ( v[axis] != 0 )
                size[axis] = simbox.guardsize;
        return size;
    }

    VecI getBorderPositionInDirection( const int & direction ) const {
        VecI pos ( simbox.guardsize );
        VecI v( getDirectionVector<T_DIM>( direction ) );
        for ( int axis=0; axis<T_DIM; ++axis )
            if ( v[axis] == 1 )
                pos[axis] += simbox.localcells[axis] - simbox.guardsize;
        assert( simbox.inArea( pos, SimulationBox::BORDER ) );
        return pos;
    }

    VecI getGuardSizeInDirection( const int & direction ) const {
        return this->getBorderSizeInDirection( direction );
    }

    VecI getGuardPositionInDirection( const int & direction ) const {
        VecI pos = this->getBorderPositionInDirection( direction );
        VecI v( getDirectionVector<T_DIM>( direction ) );
        for ( int axis=0; axis<T_DIM; ++axis ) {
            assert( abs(v[axis]) <= 1 );
            pos[axis] += v[axis]*simbox.guardsize;
        }
        return pos;
    }
#endif


    /* Start asynchronous sends of all Border Data belonging to other threads *
     * If neighbor-cell owned by this thread itself, interpolate and copy     *
     * data to guards.                                                        */
    void StartGuardUpdate( int timestep = 1 );

#if 1 == 0
    /* Waits for MPI_Recv_Requests to finish and then copies the received     *
     * partial matrix from the buffer to the simulation matrix. Because of    *
     * complicated stride in all directions it can't be communicated directly */
    void FinishGuardUpdate( int timestep = 1 ) {
        while(true) {
            int index;
            MPI_Waitany(nneighbors, &(recvrequests[1]), &index, MPI_STATUSES_IGNORE);
            /* Because we don't let MPI_Waitany watch the 0th request, which  *
             * would be a recv from the process itself, we get the index      *
             * shifted by one back !                                          */
            int direction = index + 1;

            /* this index is returned if all recvrequests are MPI_REQUEST_NULL*/
            if ( index == MPI_UNDEFINED )
                break;
            assert( direction >= 1 and direction <= this->nneighbors );
            /* delete request (neccessary?) for MPI_UNDEFINED */
            this->recvrequests[direction] = MPI_REQUEST_NULL;

            VecI size = getGuardSizeInDirection    (direction);
            VecI pos  = getGuardPositionInDirection(direction);

            #if DEBUG_COMMUNICATOR >= 2
                terr << "[Rank " << this->rank << " in ComBox] Recv from Direction ";
                     << getDirectionVector<T_DIM>( direction )
                     << "(=lin:" << direction << "=Rank:" << neighbors[direction]
                     << ") => pos: " << pos << " size: " << size << endl
                     << "Received Matrix: " << endl << recvmatrices[direction]
                     << endl << flush;
            #endif

            simbox.t[timestep]->cells.insertMatrix( pos, recvmatrices[direction] );
        }
    }
#endif

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
};

#include "Communicator.tpp"
