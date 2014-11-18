#pragma once

#include <mpi.h>
#include <cstdlib>  // malloc
#include <iostream>  // malloc
#include "Vector.h"
#include "SimulationBox.h"

using namespace std;


#define DEBUG_COMMUNICATOR 1



template<int T_DIMENSION>
class CommTopo{
private:
      /* Make Copy- and Assignment-Constructors unavailable, to prevent       *
       * copies. If there is an error at                                      *
       *   Singleton var = Singleton::getInstance                             *
       * use a reference!                                                     *
       *   Singleton & var = Singleton::getInstance                           */
      CommTopo(const CommTopo&);                // Don't Implement
      CommTopo& operator=(const CommTopo&);     // Don't implement

public:
    static CommTopo& getInstance(void)
    {
        static CommTopo instance; // Guaranteed to be destroyed.
        return instance;
    }
    ~CommTopo() {
        free( this->neighbors );
        free( this->sendrequests );
        free( this->recvrequests );
        free( this->sendmatrices );
        free( this->recvmatrices );
    };

    /**************************************************************************
     * 3 bits are not enough, they give only 8 possibilities. That wouldn't   *
     * include diagonal neighbors of all sorts (edges, corners in 3D)         *
     * Therefore we need tertiary math:                                       *
     *   0=>-1 from position, 1=>stay at pos, 2=>+1 from position along axis  *
     * So the axes in 3D become:                                              *
     *   axis 0: -1=left  , 0=stay, +1=right                                  *
     *   axis 1: -3=bottom, 0=stay, +3=top    (bottom+right <! stay)          *
     *   axis 2: -9=back  , 0=stay, +9=front                                  *
     * Difference between two states of axis n must be smaller than range     *
     * accessible by alle axes < n                                            *
     * We get the following structure:                                        *
     *                                                       back:            *
     *                                    stay:         --------------        *
     *                  front:        --------------   | 11 | 12 | 13 |       *
     *              --------------   |  2 |  3 |  4 |  |--------------|       *
     *       top:  | -7 | -6 | -5 |  |--------------|  |  8 |  9 | 10 |       *
     *             |--------------|  | -1 |  0 |  1 |  |--------------|       *
     *      stay:  |-10 | -9 | -8 |  |--------------|  |  5 |  6 |  7 |       *
     *             |--------------|  | -4 | -3 | -2 |   --------------        *
     *    bottom:  |-13 |-12 |-11 |   --------------                          *
     *              --------------                                            *
     *               ^    ^    ^                                              *
     *             left  stay right                                           *
     * As can be seen arbitrary diagonal neighbors can be easily derived by   *
     * adding instead von anding the values, like it would be the case for    *
     * a bitmask:                                                             *
     *   left + top = -1-+3 = 2   =>  correct                                 *
     * by shifting with +13 we can map this to positive numbers               *
     * These numbers although nicely to look at, proove problematic (tried to *
     * access array[-1]. It is better to have only positive values, which is  *
     * possible if we say: negative value -> 2*abs( negative value ). We get  *
     * the following structure:                                               *
     *   axis 0: 2=left  , 0=stay, +1=right                                   *
     *   axis 1: 6=bottom, 0=stay, +3=top    (bottom+right <! stay)           *
     *   axis 2: 18=back , 0=stay, +9=front                                   *
     *                                                       back:            *
     *                                    stay:         --------------        *
     *                  front:        --------------   | 23 | 21 | 22 |       *
     *              --------------   |  5 |  3 |  4 |  |--------------|       *
     *       top:  | 14 | 12 | 13 |  |--------------|  | 20 | 18 | 19 |       *
     *             |--------------|  |  2 |  0 |  1 |  |--------------|       *
     *      stay:  | 11 |  9 | 10 |  |--------------|  | 26 | 24 | 25 |       *
     *             |--------------|  |  8 |  6 |  7 |   --------------        *
     *    bottom:  | 17 | 15 | 16 |   --------------                          *
     *              --------------                                            *
     *               ^    ^    ^                                              *
     *             left  stay right                                           *
     **************************************************************************/

    /* These are the same access-specifiers like in SimulationBox.h -> merge? */
    static const int X_AXIS = 0;
    static const int Y_AXIS = 1;
    static const int Z_AXIS = 2;

    /* This basically is a 3x3 Matrix => gut use getLinearIndex from          *
     * BaseMatrix to calculate these directions on the fly. It should anyway  *
     * not be necessary for N Dimensions to use LEFT,RIGHT,...                */
    static const int LEFT   = 2;
    static const int RIGHT  = 1;
    static const int BOTTOM = 6;
    static const int TOP    = 3;
    static const int BACK   = 18;
    static const int FRONT  = 9;
    
    /* These should be changed to using BaseMatrix::getVectorIndex.           *
     * Problem: Will have to check for (0,0,0,...) at every for loop !!!      *
     *   - as all the values would be positive we also should shift it by the *
     *     position of the center, to get directions relative to it           */
    VecI getDirectionVector( int direction ) {
        VecI vec(0);
        assert( direction >= 1 and direction <= nneighbors );
        /* This may not be needed for 2D, but it also doesn't make the        *
         * results wrong                                                      */
        switch( direction % 3 ) {
            case 1: vec[X_AXIS] = +1; break;
            case 2: vec[X_AXIS] = -1; break;
        }
        direction /= 3;
        switch( direction % 3 ) {
            case 1: vec[Y_AXIS] = +1; break;
            case 2: vec[Y_AXIS] = -1; break;
        }
        direction /= 3;
        switch( direction % 3 ) {
            case 1: vec[Z_AXIS] = +1; break;
            case 2: vec[Z_AXIS] = -1; break;
        }
        return vec;
    }
    
    /* reverse of the above getDirectionVector function */
    int getDirectionNumber( const VecI v ) {
        int direction = 0;
        int prevrange = 1;
        for (int i=0; i<T_DIMENSION; i++) {
            int value = v[i];
            if (value==-1) value = 2;
            direction += value * prevrange;
            prevrange *= 3;
        }
        assert( direction >=1 and direction <=this->nneighbors );
        assert( getDirectionVector( direction ) == v );
        return direction;
    }
    
    int getOppositeDirection( int direction ) {
        return getDirectionNumber( getDirectionVector( direction )*(-1) );
    }
    
    /* returns the axis corresponding to a direction:                         *
     *   RIGHT,LEFT -> 0, TOP,BOTTOM -> 1, FRONT,BACK -> 2, ...               *
     * returns -1 if e.g. RIGHT+LEFT                                          */
    int getAxis( int direction ) {
        VecI v( getDirectionVector( direction ) );
        int axis = -1;
        for (int i=0; i<T_DIMENSION; i++)
            if ( v[i] != 0 ) {
                if (axis == -1)
                    axis = i;
                else
                    /* 2nd axis found -> on diagonal axes return -1 */
                    return -1;
            }
        return axis;
    }

    static const int dim = T_DIMENSION;

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int  processor_name_length;
    int  rank, worldsize;
    MPI_Comm  communicator;

    int  periodic[T_DIMENSION];
    int  nthreads[T_DIMENSION];
    int  coords  [T_DIMENSION];

    /* Because of reasons when accessing and direction conversion, these      *
     * also have an array element for the rank itself at [0] = (0,0,0) !!!    */
    int  nneighbors;
    int* neighbors;
    
    MPI_Request * sendrequests;
    MPI_Request * recvrequests;

    SimBox::CellMatrix * sendmatrices;
    SimBox::CellMatrix * recvmatrices;

    void StartGuardUpdate( void ) // Asynchron, returns status
    {
        /* Only cycle through direct neighbors first, see notes on paper !!!  */
        for ( int direction = 1; direction <= this->nneighbors; direction++ ) {
            int axis = getAxis(direction);
            if ( axis == -1 ) {
                recvrequests[direction] = MPI_REQUEST_NULL;
                sendrequests[direction] = MPI_REQUEST_NULL;
                continue;
            }
            
            SimBox & simbox = SimBox::getInstance();
            VecI size( simbox.localcells );
            size[axis] = simbox.guardsize;
            VecI pos ( simbox.guardsize );
            pos[axis] += simbox.localcells[axis] - 1; // we don't want a size, but a position!
            
            #if DEBUG_COMMUNICATOR >= 1
                cout << "[Rank " << this->rank << " in ComBox] Send to Direction ";
                getDirectionVector( direction ).Print();
                cout << "(=lin:" << direction << "=Rank:" << neighbors[direction] << ")";
                cout << " => pos: "; pos.Print();
                cout << " size: "; size.Print();
                cout << endl << flush;
            #endif
            
            sendmatrices[direction] = simbox.cells.getPartialMatrix( pos, size );
            MPI_Isend( sendmatrices[direction].data, sendmatrices[direction].getSize().product()
                       * sizeof(SimBox::CellMatrix::Datatype), MPI_CHAR,
                       neighbors[direction], 53+direction, communicator, &(sendrequests[direction]) );
            
            int opdir = getOppositeDirection( direction );
            #if DEBUG_COMMUNICATOR >= 2
                cout << "[Rank " << this->rank << " in ComBox] Recv from Direction ";
                getDirectionVector( opdir ).Print();
                cout << "(=lin:" << opdir << "=Rank:" << neighbors[opdir] << ")";
                cout << endl << flush;
            #endif
            
            recvmatrices[opdir] = sendmatrices[direction]; // copy because we want same size
            MPI_Irecv( recvmatrices[opdir].data, recvmatrices[opdir].getSize().product() *
                       sizeof(SimBox::CellMatrix::Datatype), MPI_CHAR,
                       neighbors[opdir], 53+opdir, communicator, &(recvrequests[opdir]) );
        }
    }

    /* Waits for MPI_Recv_Requests to finish and then copies the received     *
     * partial matrix from the buffer to the simulation matrix. Because of    *
     * complicated stride in all directions it can't be communicated directly */
    void FinishGuardUpdate( void ) {
        while(true) {
            int index;
            MPI_Waitany(nneighbors, &(recvrequests[1]), &index, MPI_STATUSES_IGNORE);
            /* Because we don't let MPI_Waitany watch the 0th request, which  *
             * would be a recv from the process itself, we get the index      *
             * shifted by one back !                                          */
            int direction = index + 1;
            #if DEBUG_COMMUNICATOR >= 2
                cout << "Got Result from Direction: " << direction << " Index: " << index << endl;
            #endif
            
            /* this index is returned if all recvrequests are MPI_REQUEST_NULL*/
            if ( index == MPI_UNDEFINED )
                break;
            assert( direction >= 1 and direction <= this->nneighbors );
            /* delete request (neccessary?) for MPI_UNDEFINED */
            this->recvrequests[direction] = MPI_REQUEST_NULL;
            
            /* if we receive from left, we can use the code for send to right *
             * for calculating pos, size !                                    */
            SimBox & simbox = SimBox::getInstance();
            int axis = getAxis(direction);
            if ( axis == -1 ) {
                cout << "Receiving from diagonal neighbor not yet implemented!" << endl;
            }
            assert( axis >= 0 and axis < T_DIMENSION);
            
            VecI size( simbox.localcells );
            size[axis] = simbox.guardsize;
            VecI pos ( simbox.guardsize );
            pos[axis] += simbox.localcells[axis];
            
            #if DEBUG_COMMUNICATOR >= 1
                cout << "[Rank " << this->rank << " in ComBox] Recv from Direction ";
                getDirectionVector( direction ).Print();
                cout << "(=lin:" << direction << "=Rank:" << neighbors[direction] << ")";
                cout << " => pos: "; pos.Print();
                cout << " size: "; size.Print();
                cout << endl << "Received Matrix: " << endl;
                    {const SimBox::CellMatrix & m = recvmatrices[direction];
                    VecI size = m.getSize();
                    VecI ind(0);
                    for ( ind[Y_AXIS]=0; ind[Y_AXIS]<size[Y_AXIS]; ind[Y_AXIS]++) {
                        for ( ind[X_AXIS]=0; ind[X_AXIS]<size[X_AXIS]; ind[X_AXIS]++)
                            cout << m[ ind ].value << " ";
                        cout << endl;
                    }
                    cout << endl;}
                cout << endl << flush;
            #endif
            
#if 1==1
            //simbox.cells[pos] = (recvmatrices[direction])[VecI(0)];
            const SimBox::CellMatrix & m = recvmatrices[direction];
            for ( int i=0; i<m.getSize().product(); i++ ) {
                //assert( pos + m.getVectorIndex(i) < this->size );
                simbox.cells[ pos + m.getVectorIndex(i) ] = m[i];
            }
            //simbox.cells.insertMatrix( pos, recvmatrices[direction] );
#endif
        }
        
        #if 1 == 0
            // MPI_Waitall( nneighbors, sendrequests, MPI_STATUSES_IGNORE );
            // sends shouldn't matter for us ...
            // MPI_Waitall( nneighbors, recvrequests, MPI_STATUSES_IGNORE );
            // waitall makes problem if rest of array are no valid requests :(
            // waitall isn't performant anyway. Better use waitany and analyze
            // whether the packet is from the left or right or so on rank and then
            // dynamicall copy buffer to simbox-matrix
            MPI_Wait( &(recvrequests[LEFT]), MPI_STATUSES_IGNORE );

            /* Copy Buffer Matrices into simBox Matrix */
            VecI pos;

            /* Receive from left */
            pos [X_AXIS] = 0;
            pos [Y_AXIS] = simbox.guardsize;
            simbox.cells.insertMatrix( pos, recvmatrices[LEFT] );

            /* Receive from right */
            MPI_Wait( &(recvrequests[RIGHT]), MPI_STATUSES_IGNORE );
            pos [X_AXIS] = simbox.guardsize + simbox.localcells[X_AXIS];
            pos [Y_AXIS] = simbox.guardsize;
            simbox.cells.insertMatrix( pos, recvmatrices[RIGHT] ); 
        #endif 
    };

    
#if DEBUG_COMMUNICATOR >= 1
    void Print( void ) {
        cout << "Printing from Communicator on rank " << this->rank;
        cout << " with cartesian coordinates: "; VecI( this->coords ).Print();
        cout << " My neighbors are: " << endl;
        
        for (int direction = 1; direction <= this->nneighbors; direction++) {
            (VecI(this->coords) + getDirectionVector( direction )).Print();
            cout << "(lin=" << direction << ")";
            cout << " -> Rank: " << this->neighbors[direction] << endl << flush;
        }
        
        /* Wait a bit till everything is flushed out, to not get scrambled output */
        double t0 = MPI_Wtime();
        while( MPI_Wtime() - t0 < 0.1 ) {}
    }
#endif
    
private:

    /* Constructor */
    CommTopo(void) {
        for (int i=0; i<T_DIMENSION; i++) {
            this->periodic[i] = true;
            this->nthreads[i] = 0;
        }

        int nneighbors = 1;
        for (int i=0; i<T_DIMENSION; i++)
            nneighbors *= 3;
        nneighbors--;
        cout << "Number of Neighbors: " << nneighbors;
        this->nneighbors = nneighbors;
        this->neighbors    = (int*)         malloc( sizeof(int)        *(nneighbors+1) );
        this->recvrequests = (MPI_Request*) malloc( sizeof(MPI_Request)*(nneighbors+1) );
        this->sendrequests = (MPI_Request*) malloc( sizeof(MPI_Request)*(nneighbors+1) );
        this->sendmatrices = (SimBox::CellMatrix*) malloc( sizeof(SimBox::CellMatrix)*(nneighbors+1) );
        this->recvmatrices = (SimBox::CellMatrix*) malloc( sizeof(SimBox::CellMatrix)*(nneighbors+1) );
        if ( neighbors == NULL or recvrequests == NULL or sendrequests == NULL )
            cerr << "Couldn't allocate Memory!";
        for (int i=0; i<=nneighbors; i++) {
            sendrequests[i] = MPI_REQUEST_NULL;
            recvrequests[i] = MPI_REQUEST_NULL;
        }
            
        /* Initialize MPI */
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
        MPI_Get_processor_name(processor_name, &processor_name_length);

        /* Create 2D-Torus Topology and get coords */
        MPI_Dims_create( worldsize, T_DIMENSION, nthreads );
        MPI_Cart_create( MPI_COMM_WORLD, T_DIMENSION, nthreads, periodic,
                         true /* reorder */ , &communicator);
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_Cart_coords( communicator, rank, T_DIMENSION, coords );

        /* Get Neighbors from Topology */
        for (int direction = 1; direction <= this->nneighbors; direction++) {
            /* get rank from cartesian coords -> enables diagonal neighbors ! */
            MPI_Cart_rank( this->communicator, ( VecI(this->coords) + 
                getDirectionVector( direction ) ), &(neighbors[direction]));
        }
    }

};
