#pragma once

#include "yeesolver/YeeCell.h"
#include "simbox/SimulationBox.h"
#include "math/TVector.h"

namespace YeeSolver {

typedef SimulationBox::SimulationBox<SIMDIM,YeeCell> OctCell;
typedef Vec<int,SIMDIM> VecI;
const int X = 0;
const int Y = 1;
const int Z = 2;

void CalcH( OctCell & simbox, int tnext, int tcur, int area, double deltaT ) {
    typename OctCell::CellMatrix CellMatrix;
    typename OctCell::IteratorType itnext = simbox.getIterator( tnext, area );
    typename OctCell::IteratorType itcur  = simbox.getIterator( tcur , area );

    Vec<double,3> cellsize( (SIMDIM>=1) ? simbox.cellsize[0] : INF ,
                            (SIMDIM>=2) ? simbox.cellsize[1] : INF ,
                            (SIMDIM>=3) ? simbox.cellsize[2] : INF );

    for ( itnext=itnext.begin(); itnext != itnext.end(); ++itnext ) {
        /* we misuse the iterator as a fancy access operator. icell is the    *
         * only member changed when iterating, so copying that suffices.      *
         * copying the whole iterator would be wrong, because the shall       *
         * have different srcmats like given by getIterator( timestep, ... )  */
        itcur.icell = itnext.icell;

        VecI xnext(0); if (SIMDIM>=1) xnext[X]++;
        VecI ynext(0); if (SIMDIM>=2) ynext[Y]++;
        VecI znext(0); if (SIMDIM>=3) znext[Z]++;

        /* Update all H components */
        double nom = (1.0 + 0.5*itnext->sigmaM*deltaT / itnext->mu);
        double Da  = (1.0 - 0.5*itnext->sigmaM*deltaT / itnext->mu) / nom;
        double Dbx = deltaT / ( itnext->mu * cellsize[0] ) / nom;
        double Dby = deltaT / ( itnext->mu * cellsize[1] ) / nom;
        double Dbz = deltaT / ( itnext->mu * cellsize[2] ) / nom;

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

void CalcE( OctCell & simbox, int tnext, int tcur, int area, double deltaT ) {
    typename OctCell::CellMatrix CellMatrix;
    typename OctCell::IteratorType itnext = simbox.getIterator( tnext, area );
    typename OctCell::IteratorType itcur  = simbox.getIterator( tcur , area );

    Vec<double,3> cellsize( (SIMDIM>=1) ? simbox.cellsize[0] : INF ,
                            (SIMDIM>=2) ? simbox.cellsize[1] : INF ,
                            (SIMDIM>=3) ? simbox.cellsize[2] : INF );

    for ( itnext=itnext.begin(); itnext != itnext.end(); ++itnext ) {
        itcur.icell = itnext.icell;

        VecI xprev(0); if (SIMDIM>=1) xprev[X]--;
        VecI yprev(0); if (SIMDIM>=2) yprev[Y]--;
        VecI zprev(0); if (SIMDIM>=3) zprev[Z]--;

		/* Now update all E components */
        double nom = (1.0 + 0.5*itnext->sigmaE*deltaT / itnext->epsilon);
        double Ca  = (1.0 - 0.5*itnext->sigmaE*deltaT / itnext->epsilon) / nom;
        double Cbx = deltaT / ( itnext->epsilon * cellsize[0] ) / nom;
        double Cby = deltaT / ( itnext->epsilon * cellsize[1] ) / nom;
        double Cbz = deltaT / ( itnext->epsilon * cellsize[2] ) / nom;

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
