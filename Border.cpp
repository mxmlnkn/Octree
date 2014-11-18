/*

rm border.exe; mpic++ border.cpp -o border.exe -Wall -std=c++0x; mpirun -n 4 ./border.exe

ToDo:
  - Implement Communication (will have to link ComBox with SimBox :(, maybe
    even as global vars. Only hope it can be done without cyclic referencing
  - Do more MPI Error Management
  - Watch Out for Matrices indicey (row,col) vs. (x,y) (Remove Dirty hack!)
  - Run Test-Kernel (Add 4 nearest neighbors up and mod 10 ) to show how it is
    done with the iterator
  - Test every method / function singelhandendly for correctness !
  - Test for only two matrices (send to itself)
  - Test in 3D ... (need better output)
  - include per rank output with tout (extra file)
  - something is wrong with matrix copy assignment !
        const const SimBox::CellMatrix & m = recvmatrices[direction];
        SimBox::CellMatrix m = recvmatrices[direction];
        SimBox::CellMatrix m( recvmatrices[direction] );
    only first one does work ...
  - look which methods we could set to static in the classes ( shouldn't depend
    on varibales inside the class, only on the template parameters ! )
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

#include "Vector.h"
typedef Vec<double, SIMDIM> VecD;
typedef Vec<int   , SIMDIM> VecI;

#include "SimulationBox.h"
typedef SimulationBox<SIMDIM,GUARDSIZE> SimBox; // neded by Communicator.h

#include "Communicator.h"


using namespace std;




int main( int argc, char **argv )
{
    CommTopo<SIMDIM> & comBox = CommTopo<SIMDIM>::getInstance();
    SimBox & simBox = SimBox::getInstance();

    /* Init( VecD globsize, VecI globalcells, VecI localcells, int mpicoords[T_DIMENSION] ) */
    simBox.Init( 10, 10, 5, comBox.coords );

    SimBox::Iterator it = simBox.getIterator( simBox.CORE + simBox.BORDER );
    VecD globcoords = simBox.getGlobalPosition( it );
    
    for ( it=it.begin(); it!=it.end(); ++it ) {
        VecD globsize   = simBox.globsize;
        VecD globcoords = simBox.getGlobalPosition( it );
        simBox[it].value = (5 +
          int( 5*sin(2.*M_PI * globcoords[0]/globsize[0] ) *
                 sin(2.*M_PI * globcoords[1]/globsize[1] )
          ) ) % 10;
        // simBox[it].value = globcoords[1];
    }
    
    const int rank      = comBox.rank;
    const int worldsize = comBox.worldsize;
    const MPI_Comm commTorus = comBox.communicator;

#if DEBUG == 1
    {if (rank == 0) {
        cout << "Torus Size: (";
        for (int j=0; j<SIMDIM-1; j++)
            cout << comBox.nthreads[j] << ",";
        cout << comBox.nthreads[SIMDIM-1] << ")" << endl;
    }

    MPI_Barrier( commTorus );
    int beacon;
    if ( rank == 0 )
    {
        comBox.Print();
        cout << "Global Coordinates: "; globcoords.Print(); cout << endl;
        simBox.PrintValues();
        /* Let other ranks print their Matrices */
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
    }
    else
    {
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
        comBox.Print();
        cout << "Global Coordinates: "; globcoords.Print(); cout << endl;
        simBox.PrintValues();
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
    }
    /* Wait a bit till everything is flushed out, to not get scrambled output */
    double t0 = MPI_Wtime();
    while( MPI_Wtime() - t0 < 0.5 ) {}
    MPI_Barrier( commTorus );}
    
    /* Test getPartialMatrix */
    if (rank == 0) {
        cout << "==============" << endl;
        cout << "Update Guards!" << endl;
        cout << "==============" << endl << flush;
        int size[SIMDIM] = {3,4};
        int  pos[SIMDIM] = {0,3};
        BaseMatrix<CellData,SIMDIM> mat = simBox.cells.getPartialMatrix( VecI(pos), VecI(size) );
        cout << "PartialMatrix at "; VecI(pos).Print();
        cout << " of size "; VecI(size).Print();
        cout << " in matrix of rank 0: " << endl;
        for (int i=0; i<size[0]; i++) {
            for (int j=0; j<size[1]; j++) {
                int ind[SIMDIM] = {i,j};
                cout << mat[VecI(ind)].value << " ";
            }
            cout << endl << flush;
        }
    }
    
    /* Wait a bit till everything is flushed out, to not get scrambled output */
    double t0 = MPI_Wtime();
    while( MPI_Wtime() - t0 < 0.5 ) {}
    MPI_Barrier( commTorus );
#endif


    comBox.StartGuardUpdate(); // Asynchron, returns status
    comBox.FinishGuardUpdate();
    MPI_Barrier( commTorus );

#if DEBUG == 1
    {MPI_Barrier( commTorus );
    int beacon;
    if ( rank == 0 )
    {
        cout << "[Rank " << rank << "]" << endl;
        cout << "MPI-Coords: "; VecD( comBox.coords ).Print(); cout << endl;
        cout << "Global Coordinates: "; globcoords.Print(); cout << endl;
        simBox.PrintValues();
        /* Let other ranks print their Matrices */
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
    }
    else
    {
        MPI_Recv( &beacon, 1, MPI_INT, (rank-1)%worldsize, beacon, commTorus, MPI_STATUS_IGNORE );
        cout << "[Rank " << rank << "]" << endl;
        cout << "MPI-Coords: "; VecD( comBox.coords ).Print(); cout << endl;
        cout << "Global Coordinates: "; globcoords.Print(); cout << endl;
        simBox.PrintValues();
        MPI_Send( &beacon, 1, MPI_INT, (rank+1)%worldsize, beacon, commTorus );
    }
    MPI_Barrier( commTorus );}
#endif
    
    MPI_Finalize(); // doesn't work in destructor :S
}

/*
Done:
  - Find better encoding for the directions (get rid of CENTER / OFFSET ! )
*/