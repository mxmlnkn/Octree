#pragma once

#include "yeesolver/YeeCell.h"
#include "simbox/SimulationBox.h"
#include "math/TVector.h"

namespace YeeSolver {

typedef SimulationBox::SimulationBox<SIMDIM,YeeCell> OctCell;
typedef Vec<int,2> VecI;
const int X = 0;
const int Y = 0;
const int Z = 0;

void CalcH( OctCell & simBox, int tnext, int tcur, int area ) {
    typename OctCell::CellMatrix CellMatrix;
    typename OctCell::IteratorType itnext = simBox.getIterator( tnext, area );
    typename OctCell::IteratorType itcur  = simBox.getIterator( tcur , area );
    for ( itnext=itnext.begin(); itnext != itnext.end(); ++itnext ) {
        /* we misuse the iterator as a fancy access operator. icell is the    *
         * only member changed when iterating, so copying that suffices.      *
         * copying the whole iterator would be wrong, because the shall       *
         * have different srcmats like given by getIterator( timestep, ... )  */
        itcur.icell = itnext.icell;
        
        VecI xnext(0); xnext[X]++;
        VecI ynext(0); ynext[Y]++;
        VecI znext(0); znext[Z]++;

        /* Update all H components */
        double nom = (1.0 + 0.5*itnext->sigmaM*DELTA_T / itnext->mu);
        double Da  = (1.0 - 0.5*itnext->sigmaM*DELTA_T / itnext->mu) / nom;
        double Dbx = DELTA_T / ( itnext->mu * CELL_SIZE_X ) / nom;
        double Dby = DELTA_T / ( itnext->mu * CELL_SIZE_Y ) / nom;
        double Dbz = DELTA_T / ( itnext->mu * CELL_SIZE_Z ) / nom;

        #if DEBUG_MAIN_YEE >= 100
            if ( (itcur+ynext)->E != itcur->E )
                tout << "!!!! E[" << ynext << "]=" << (itcur+ynext)->E << " != E[" << pos << "]=" << itcur->E << "\n";
        #endif

        itnext->H[X] = Da * itcur->H[X] +
             Dbz*( (itcur+znext)->E[Y] - itcur->E[Y] )
            -Dby*( (itcur+ynext)->E[Z] - itcur->E[Z] );
        itnext->H[Y] = Da * itcur->H[Y] +
             Dbx*( (itcur+xnext)->E[Z] - itcur->E[Z] )
            -Dbz*( (itcur+znext)->E[X] - itcur->E[X] );
        itnext->H[Z] = Da * itcur->H[Z] +
             Dby*( (itcur+ynext)->E[X] - itcur->E[X] )
            -Dbx*( (itcur+xnext)->E[Y] - itcur->E[Y] );
    }

}

void CalcE( OctCell & simBox, int tnext, int tcur, int area ) {
    typename OctCell::CellMatrix CellMatrix;
    typename OctCell::IteratorType itnext = simBox.getIterator( tnext, area );
    typename OctCell::IteratorType itcur  = simBox.getIterator( tcur , area );
    for ( itnext=itnext.begin(); itnext != itnext.end(); ++itnext ) {
        itcur.icell = itnext.icell;
        
        VecI xprev(0); xprev[X]--;
        VecI yprev(0); yprev[Y]--;
        VecI zprev(0); zprev[Z]--;

		/* Now update all E components */
        double nom = (1.0 + 0.5*itnext->sigmaE*DELTA_T / itnext->epsilon);
        double Ca  = (1.0 - 0.5*itnext->sigmaE*DELTA_T / itnext->epsilon) / nom;
        double Cbx = DELTA_T / ( itnext->epsilon * CELL_SIZE_X ) / nom;
        double Cby = DELTA_T / ( itnext->epsilon * CELL_SIZE_Y ) / nom;
        double Cbz = DELTA_T / ( itnext->epsilon * CELL_SIZE_Z ) / nom;

        itnext->E[X] = Ca * itcur->E[X] +
             Cby*( itnext->H[Z] - (itnext+yprev)->H[Z] )
            -Cbz*( itnext->H[Y] - (itnext+zprev)->H[Y] );
        itnext->E[Y] = Ca * itcur->E[Y] +
             Cbz*( itnext->H[X] - (itnext+zprev)->H[X] )
            -Cbx*( itnext->H[Z] - (itnext+xprev)->H[Z] );
        itnext->E[Z] = Ca * itcur->E[Z] +
             Cbx*( itnext->H[Y] - (itnext+xprev)->H[Y] )
            -Cby*( itnext->H[X] - (itnext+yprev)->H[X] );
    }
}

} // namespace YeeSolver
