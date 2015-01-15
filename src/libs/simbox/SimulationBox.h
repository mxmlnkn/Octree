#pragma once

#include <cstdlib>
#include <iostream>
#include "math/TVector.h"
#include "math/TBaseMatrix.h"
#include "simbox/SimulationBoxDefines.h"
#include "simbox/SimulationBoxIterator.h"

#define DEBUG_SIMBOX 1
#define TIMESTEPS_NEEDED_FOR_CALCULATION 2


namespace SimulationBox {

/**************************************************************************
 * The 'Halo' in Musubi is here being called Guard for consistency with   *
 * PicOnGPU. Musbu doesn't differ between Core and Border. In PicOnGPU    *
 * these two are differentiated in order to do calculations on the Core,  *
 * while sending the data from the Border to the respective Guard of      *
 * processes.                                                             *
 * The Guard must be as wide as we need neighbor cells for one calcula-   *
 * tion in one cell. The Border is as wide as the Guard by definition.    *
 * Therefore cells in the Core can be calculated by only including the    *
 * annexing Core- or Bordercells!                                         *
 *                                             Rank 1                     *
 *                                        -----------------               *
 *                                       | G   G G G G   G |              *
 *                   Rank 0              |    ---------    |              *
 *              -----------------        | G | B B B B | G |              *
 *             | G   G G G G   G |       | G | B C C B | G |              *
 *             |    ---------    |       | G | B C C B | G |              *
 *             | G | B B B B | G |       | G | B B B B | G |              *
 *             | G | B C C B | G |       |    ---------    |              *
 *             | G | B C C B | G |       | G   G G G G   G |              *
 *             | G | B B B B | G |        -----------------               *
 *             |    ---------    |                                        *
 *             | G   G G G G   G |                                        *
 *              -----------------                                         *
 * In Order to do calculations on these areas easily, iterators and       *
 * relative cell-directions need to be included                           *
 * The internal data matrix will begin with the upper left corner like    *
 * any normal matrix                                                      *
 **************************************************************************/
bool InArea( const VecI & pos, const int & area, const VecI & localcells, int guardsize );


template<int T_DIM, typename T_CELLDATA>
class SimulationBox {
public:
    SimulationBox( void );

    /* Definitions and Aliases */
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    typedef BaseMatrix<T_CELLDATA,T_DIM> CellMatrix;
    typedef Iterator<T_DIM> IteratorType;
    static const int dim       = T_DIM;

    /* Current Time on which Calculations are to be done or where done are    *
     * held in t[0]. t[1] is the prior time step and so on. TimeData struct   *
	 * introduced so that it looks better: t[0]->cells[0] instead of          *
	 * problematic t[0][0]                                                    */
private:
    int ntimesteps;
    typedef struct TimeDataStruct {
        CellMatrix cells;
    } TimeData;
public:
    /* holds array to different versions of matrix, e.g. different times      *
     * array of pointers is being used instead of simple array in order to    *
     * cyclically swap the time pointers, instead of the whole data itself !  */
    TimeData ** t;

    int guardsize;
    VecD cellsize;
    VecI localcells; // cells in this SimulationBox/Matrix/Octree cell
    VecD abspos;  // of lower (or better upper ???) left point (without Guard)

    ~SimulationBox(void);

    bool inArea( const VecI & pos, const int & area ) const;
    IteratorType getIterator( const int area ) const;

    void init( VecD abspos, VecI localcells, VecD cellsize, int guardsize, int bufferpages );
    #if DEBUG_SIMBOX >= 1
        void PrintValues( int timestep = 0 ) const;
    #endif

    /* returns global position of cell by using global offset of the matrix,
     * */
    VecD getGlobalPosition ( const IteratorType it ) const;
    void copyCurrentToPriorTimestep( void );
};

} // namespace SimulationBox

#include "SimulationBox.tpp"