/*

rm Main.exe; mpic++ Main.cpp -o Main.exe -Wall -std=c++0x; mpirun -n 4 ./Main.exe

ToDo:
  - Do more MPI Error Management
  - Watch Out for Matrices indices (row,col) vs. (x,y)
  - Run Test-Kernel (Add 4 nearest neighbors up and mod 10 ) to show how it is
    done with the iterator
  - Test every method / function singelhandendly for correctness !
        -> Done for Vector, BaseMatrix
  - Test in 3D ... (need better output)
  - Use Functions of 3x3 Matrix-class for Direction conversion
  - Document everything better
  - Beende den template-Wahn -.-. Vor allem Bei Vector und Basematrix mit T_DIM!

Testcase 1:

[Rank 0]                             [Rank 1]
MPI-Coords: (0,0)                    MPI-Coords: (1,0)
Global Coordinates: (0,0)            Global Coordinates: (5,0)
My Raw Matrix with Guard(size=1) is  My Raw Matrix with Guard(size=1) is
0 5 5 5 5 5 0                        0 5 5 5 5 5 0
5 6 7 7 6 5 4                        5 4 3 3 4 5 6
5 7 9 9 7 5 3                        5 3 1 1 3 5 7
5 7 9 9 7 5 3                        5 3 1 1 3 5 7
5 6 7 7 6 5 6                        5 4 3 3 4 5 4
5 5 5 5 5 5 5                        5 5 5 5 5 5 5
0 4 3 3 4 5 0                        0 6 7 7 6 5 0


Testcase 2:
Global Coordinates: (0,0)            Global Coordinates: (5,0)
My Raw Matrix with Guard(size=1) is  My Raw Matrix with Guard(size=1) is
0 0 0 0 0 0 0                        0 0 0 0 0 0 0
0 6 7 7 6 5 0                        0 4 3 3 4 5 0
0 7 9 9 7 5 0                        0 3 1 1 3 5 0
0 7 9 9 7 5 0                        0 3 1 1 3 5 0
0 6 7 7 6 5 0                        0 4 3 3 4 5 0
0 5 5 5 5 5 0                        0 5 5 5 5 5 0
0 0 0 0 0 0 0                        0 0 0 0 0 0 0

[Rank 1]                             [Rank 1]
MPI-Coords: (1,0)                    MPI-Coords: (1,0)
Global Coordinates: (5,0)            Global Coordinates: (5,0)
My Raw Matrix with Guard(size=1) is  My Raw Matrix with Guard(size=1) is
5 5 5 5 5 5 5                        5 5 5 5 5 5 5
5 6 7 7 6 5 4                        5 4 3 3 4 5 6
5 7 9 9 7 5 3                        5 3 1 1 3 5 7
5 7 9 9 7 5 3                        5 3 1 1 3 5 7
5 6 7 7 6 5 4                        5 4 3 3 4 5 6
5 5 5 5 5 5 5                        5 5 5 5 5 5 5
5 6 7 7 6 5 4                        5 4 3 3 4 5 6


Testcase 3:

Global Coordinates: (0,0)            Global Coordinates: (5,0)
My Raw Matrix with Guard(size=2) is  My Raw Matrix with Guard(size=2) is
0 0 0 0 0 0 0 0 0                    0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0                    0 0 0 0 0 0 0 0 0
0 0 9 9 7 5 3 0 0                    0 0 1 1 3 5 7 0 0
0 0 9 9 7 5 3 0 0                    0 0 1 1 3 5 7 0 0
0 0 7 7 6 5 4 0 0                    0 0 3 3 4 5 6 0 0
0 0 5 5 5 5 5 0 0                    0 0 5 5 5 5 5 0 0
0 0 3 3 4 5 6 0 0                    0 0 7 7 6 5 4 0 0
0 0 0 0 0 0 0 0 0                    0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0                    0 0 0 0 0 0 0 0 0

[Rank 0]                             [Rank 1]
MPI-Coords: (0,0)                    MPI-Coords: (1,0)
Global Coordinates: (0,0)            Global Coordinates: (5,0)
My Raw Matrix with Guard(size=2) is  My Raw Matrix with Guard(size=2) is
5 5 5 5 5 5 5 5 5                    5 5 5 5 5 5 5 5 5
5 4 3 3 4 5 6 7 7                    5 6 7 7 6 5 4 3 3
5 7 9 9 7 5 3 1 1                    5 3 1 1 3 5 7 9 9
5 7 9 9 7 5 3 1 1                    5 3 1 1 3 5 7 9 9
5 6 7 7 6 5 4 3 3                    5 4 3 3 4 5 6 7 7
5 5 5 5 5 5 5 5 5                    5 5 5 5 5 5 5 5 5
5 4 3 3 4 5 6 7 7                    5 6 7 7 6 5 4 3 3
5 7 9 9 7 5 3 1 1                    5 3 1 1 3 5 7 9 9
5 7 9 9 7 5 3 1 1                    5 3 1 1 3 5 7 9 9


*/


#ifndef M_PI
#   define M_PI 3.14159265358979323846
#endif

#define DEBUG 1
#define GUARDSIZE 1
#define SIMDIM 2

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc
#include <mpi.h>
#include <time.h>

#include "teestream/TeeStream.h"

#include "math/TVector.h"
typedef Vec<double, SIMDIM> VecD;
typedef Vec<int   , SIMDIM> VecI;

#include "Cell.h"
#include "simbox/SimulationBox.h"
typedef SimulationBox::SimulationBox<SIMDIM,CellData> SimBox; // needed by Communicator.h
typedef typename SimBox::CellMatrix SimMatrix;
#include "Communicator.h"

using namespace std;


int main( int argc, char **argv )
{
    SimBox simBox;
    CommTopo<SIMDIM> comBox(simBox);
    tout.Open( string(""), comBox.rank );
    srand( clock() * comBox.rank );

    /* init( VecD abspos, VecI localcells, VecD cellsize, int guardsize, int bufferpages ) */
    VecD abspos;
    abspos[0] = comBox.coords[0];
    abspos[1] = comBox.coords[1];
    simBox.init( abspos*5, 5, 1, GUARDSIZE, 3 );

    SimBox::IteratorType it = simBox.getIterator( SimulationBox::CORE + SimulationBox::BORDER );
    VecD globcoords = simBox.getGlobalPosition( it );

    for ( it=it.begin(); it!=it.end(); ++it ) {
        VecD globcoords = simBox.getGlobalPosition( it );
        (simBox.t[0]->cells)[it.icell].value = (5 +
          int( 5*sin(2.*M_PI * globcoords[0]/2. ) *
                 sin(2.*M_PI * globcoords[1]/2. )
          ) ) % 10;
        // simBox[it].value = globcoords[1];
    }

    // Set guard to 0 (not really necessary, only looks better)
    for ( it=simBox.getIterator( SimulationBox::GUARD ).begin(); it!=it.end(); ++it )
        simBox.t[0]->cells[it.icell].value = 0;

    const int rank      = comBox.rank;
    const int worldsize = comBox.worldsize;
    const MPI_Comm commTorus = comBox.communicator;

#if DEBUG == 1
    if (rank == 0)
        tout << "Torus Size: (" << VecI(comBox.nthreads) << endl;

    {MPI_Barrier( commTorus );
    int beacon;
    if ( rank == 0 )
    {
        comBox.Print();
        tout << "Global Coordinates: " << globcoords << endl;
        simBox.PrintValues();
        /* Let other ranks print their Matrices */
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
    }
    else
    {
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
        comBox.Print();
        tout << "Global Coordinates: " << globcoords << endl;
        simBox.PrintValues();
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
    }
    /* Wait a bit till everything is flushed out, to not get scrambled output */
    double t0 = MPI_Wtime();
    while( MPI_Wtime() - t0 < 0.5 ) {}
    MPI_Barrier( commTorus );}

    /* Test getPartialMatrix */
    if (rank == 0) {
        tout << "==============" << endl
             << "Update Guards!" << endl
             << "==============" << endl << flush;
        int size[SIMDIM] = {3,4};
        int  pos[SIMDIM] = {0,3};
        SimMatrix mat = simBox.t[0]->cells.getPartialMatrix( VecI(pos), VecI(size) );
        tout << "PartialMatrix at " << VecI(pos) << " of size " << VecI(size)
             << " in matrix of rank 0: " << endl;
        for (int i=0; i<size[0]; i++) {
            for (int j=0; j<size[1]; j++) {
                int ind[SIMDIM] = {i,j};
                tout << mat[VecI(ind)].value << " ";
            }
            tout << endl << flush;
        }
    }

    /* Wait a bit till everything is flushed out, to not get scrambled output */
    double t0 = MPI_Wtime();
    while( MPI_Wtime() - t0 < 0.5 ) {}
    MPI_Barrier( commTorus );
#endif

    comBox.StartGuardUpdate ( 0 ); // Asynchron, returns status
    comBox.FinishGuardUpdate( 0 );
    MPI_Barrier( commTorus );

#if DEBUG == 1
    {MPI_Barrier( commTorus );
    int beacon;
    if ( rank == 0 )
    {
        tout << "[Rank " << rank << "]" << endl
             << "MPI-Coords: " << VecD( comBox.coords ) << endl
             << "Global Coordinates: " << globcoords << endl;
        simBox.PrintValues();
        /* Let other ranks print their Matrices */
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
    }
    else
    {
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
        tout << "[Rank " << rank << "]" << endl
             << "MPI-Coords: " << VecD( comBox.coords ) << endl
             << "Global Coordinates: " << globcoords << endl;
        simBox.PrintValues();
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
    }
    MPI_Barrier( commTorus );}
#endif

    tout << "==============================================" << endl
         << "Delete Guard to demonstrate asynchronous send!" << endl
         << "==============================================" << endl << flush;
    for ( SimBox::IteratorType it=simBox.getIterator( SimulationBox::GUARD ).begin(); it!=it.end(); ++it )
        (simBox.t[0]->cells)[it.icell].value = 0;
    simBox.PrintValues();

    tout << "========================================" << endl
         << "Add all 4 (2D) neighbors to cell on Core" << endl
         << "========================================" << endl << flush;
    simBox.copyCurrentToPriorTimestep();
    comBox.StartGuardUpdate( 1 ); //  Sends Border data from last timestep (t[1])
    for ( SimBox::IteratorType it=simBox.getIterator( SimulationBox::CORE ).begin(); it!=it.end(); ++it ) {
        (simBox.t[0]->cells)[it.icell].value +=
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( RIGHT  ) ].value +
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( LEFT   ) ].value +
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( TOP    ) ].value +
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( BOTTOM ) ].value;
    }
    comBox.FinishGuardUpdate( 1 ); // Saves received data into t[1], because we only need prior time step
    for ( SimBox::IteratorType it=simBox.getIterator( SimulationBox::BORDER ).begin(); it!=it.end(); ++it ) {
        //VecI neighbor = it.icell + getDirectionVector( RIGHT, SIMDIM );
        (simBox.t[0]->cells)[it.icell].value +=
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( RIGHT  ) ].value +
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( LEFT   ) ].value +
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( TOP    ) ].value +
            (simBox.t[1]->cells)[ it.icell + getDirectionVector<SIMDIM>( BOTTOM ) ].value;
    }
    simBox.PrintValues();
    
    
    tout << "==========================" << endl
         << "Set Border,Core and Guard!" << endl
         << "==========================" << endl << flush;
    for ( SimBox::IteratorType it=simBox.getIterator( SimulationBox::CORE ).begin(); it!=it.end(); ++it )
        (simBox.t[0]->cells)[it.icell].value = 1;
    for ( SimBox::IteratorType it=simBox.getIterator( SimulationBox::BORDER ).begin(); it!=it.end(); ++it )
        (simBox.t[0]->cells)[it.icell].value = 8;
    for ( SimBox::IteratorType it=simBox.getIterator( SimulationBox::GUARD ).begin(); it!=it.end(); ++it )
        (simBox.t[0]->cells)[it.icell].value = 0;

    simBox.PrintValues();

    MPI_Finalize(); // doesn't work in destructor :S

}






/*
Done:
  - Use inArea of SimBox from Iterator, instead of it's own one
        => have global InArea and this->inAres which gives some extra parameters
           from its object
  - include per rank output with tout (extra file)
  - look which methods we could set to static in the classes ( shouldn't depend
    on variables inside the class, only on the template parameters ! )
        => none really ...
  - something is wrong with matrix copy assignment !
        const const SimBox::CellMatrix & m = recvmatrices[direction];
        SimBox::CellMatrix m = recvmatrices[direction];
        SimBox::CellMatrix m( recvmatrices[direction] );
    only first one does work -> called Destructor instead of just free(...)
       => not good
  - Test for only two matrices (send to itself) -> workds
  - Implement Communication (will have to link ComBox with SimBox :(, maybe
    even as global vars. Only hope it can be done without cyclic referencing
  - Find better encoding for the directions (get rid of CENTER / OFFSET ! )
*/