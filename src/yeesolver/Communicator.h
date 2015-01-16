#pragma once

#include <cstdlib>   // malloc
#include <iostream>
#include <mpi.h>
#include "math/TVector.h"
#include "octree/Octree.h"
#include "simbox/SimulationBox.h"
#include "Directions.h"

#ifndef DEBUG_COMMUNICATOR
    #define DEBUG_COMMUNICATOR 1
#endif

template<typename T_OCTREE>
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

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int  processor_name_length;
    int  rank, worldsize;

    T_OCTREE & tree;
#if 1==0
    /* Because of reasons when accessing and direction conversion, these      *
     * also have an array element for the rank itself at [0] = (0,0,0) !!!    */
    int  nneighbors;
    int* neighbors;

    MPI_Request * sendrequests;
    MPI_Request * recvrequests;
    
    SimBox::CellMatrix * sendmatrices;
    SimBox::CellMatrix * recvmatrices;
#endif

#if 1==0
    void StartGuardUpdate( int timestep = 1 ) // Asynchron, returns status
    {
        /* Only cycle through direct neighbors first, see notes on paper !!!  */
        for ( int direction = 1; direction <= this->nneighbors; direction++ )
        {
            /* Limit send directions to very nearest neighbors, no diagonals  */
            /* if ( getAxis(direction) == -1 ) {
                recvrequests[direction] = MPI_REQUEST_NULL;
                sendrequests[direction] = MPI_REQUEST_NULL;
                continue;
            } */

            VecI size = this->getBorderSizeInDirection    (direction);
            VecI pos  = this->getBorderPositionInDirection(direction);

            sendmatrices[direction] = simbox.t[timestep]->cells.getPartialMatrix( pos, size );
            MPI_Isend( sendmatrices[direction].data, sendmatrices[direction].getSize().product()
                       * sizeof(SimBox::CellMatrix::Datatype), MPI_CHAR,
                       neighbors[direction], 53+direction, communicator, &(sendrequests[direction]) );

            #if DEBUG_COMMUNICATOR >= 1
                terr << "[Rank " << this->rank << " in ComBox] Send to Direction ";
                terr << getDirectionVector<T_DIM>( direction );
                terr << "(=lin:" << direction << "=Rank:" << neighbors[direction] << ")";
                terr << " => pos: " << pos << " size: " << size << endl << flush;
                terr << endl << "Sent Matrix: " << endl;
                    {const SimBox::CellMatrix & m = sendmatrices[direction];
                    VecI size = m.getSize();
                    VecI ind(0);
                    for ( ind[Y_AXIS]=0; ind[Y_AXIS]<size[Y_AXIS]; ind[Y_AXIS]++) {
                        for ( ind[X_AXIS]=0; ind[X_AXIS]<size[X_AXIS]; ind[X_AXIS]++)
                            terr << m[ ind ].value << " ";
                        terr << endl;
                    }
                    terr << endl;}
                terr << endl << flush;
            #endif

            int opdir = getOppositeDirection<T_DIM>( direction );
            recvmatrices[opdir] = sendmatrices[direction]; // copy because we want same size
            MPI_Irecv( recvmatrices[opdir].data, recvmatrices[opdir].getSize().product() *
                       sizeof(SimBox::CellMatrix::Datatype), MPI_CHAR,
                       neighbors[opdir], 53+direction, communicator, &(recvrequests[opdir]) );
        }
    }

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

    /******************************* Constructor ******************************/
    OctreeCommunicator(void) {
        std::cerr << "Communicator needs to be given an octree handle when constructed!\n";
    }

    OctreeCommunicator( T_OCTREE & tree ) : rank(0), worldsize(1), tree(tree) {
        /* Initialize MPI */
        /*MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
        MPI_Get_processor_name(processor_name, &processor_name_length);*/
#if 1==0
        /* TODO */
        nneighbors = 0;
        terr << "Number of neighboring processes     : " << nneighbors;
        terr << "Border size to send ( per process ) : ";
        
        this->nneighbors = nneighbors;
        this->neighbors    = (int*)         malloc( sizeof(int)        *(nneighbors+1) );
        this->recvrequests = (MPI_Request*) malloc( sizeof(MPI_Request)*(nneighbors+1) );
        this->sendrequests = (MPI_Request*) malloc( sizeof(MPI_Request)*(nneighbors+1) );
        this->sendmatrices = (SimBox::CellMatrix*) malloc( sizeof(SimBox::CellMatrix)*(nneighbors+1) );
        this->recvmatrices = (SimBox::CellMatrix*) malloc( sizeof(SimBox::CellMatrix)*(nneighbors+1) );
        if ( neighbors == NULL or recvrequests == NULL or sendrequests == NULL )
            terr << "Couldn't allocate Memory!";
        for (int i=0; i<=nneighbors; i++) {
            sendrequests[i] = MPI_REQUEST_NULL;
            recvrequests[i] = MPI_REQUEST_NULL;
        }
#endif
    }

    /******************************** Destructor ******************************/
    ~OctreeCommunicator() {
#if 1==0
        free( this->neighbors );
        free( this->sendrequests );
        free( this->recvrequests );
        free( this->sendmatrices );
        free( this->recvmatrices );
#endif
    }
};
