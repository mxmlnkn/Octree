#pragma once

#include "YeeCell.h"
#include "simbox/SimulationBox.h"
#include "math/TVector.h"

namespace YeeSolver {

typedef SimulationBox::SimulationBox<SIMDIM,YeeCell> OctCell;
typedef Vec<int,SIMDIM> VecI;
const int X = 0;
const int Y = 1;
const int Z = 2;

void CalcH( OctCell & simbox, int tnext, int tcur, int area, double dt ) {
    typename OctCell::CellMatrix CellMatrix;
    typename OctCell::IteratorType itnext = simbox.getIterator( tnext, area );
    typename OctCell::IteratorType itcur  = simbox.getIterator( tcur , area );

    Vec<double,3> cellsize( (SIMDIM>=1) ? simbox.cellsize[0] : INF ,
                            (SIMDIM>=2) ? simbox.cellsize[1] : INF ,
                            (SIMDIM>=3) ? simbox.cellsize[2] : INF );
    // ToDo: Check for Courant Criterium (now as we can calculate c_M :3) !!!

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
        double nom = (1.0 + 0.5*itnext->sigmaM*dt / itnext->mu);
        double Da  = (1.0 - 0.5*itnext->sigmaM*dt / itnext->mu) / nom;
        if ( itnext->sigmaM > sqrt(DBL_MAX) )
            Da = (1.0/itnext->sigmaM - 0.5 / itnext->mu) / (1.0/itnext->sigmaM + 0.5*dt / itnext->mu);
        double Dbx = dt / ( itnext->mu * cellsize[0] ) / nom;
        double Dby = dt / ( itnext->mu * cellsize[1] ) / nom;
        double Dbz = dt / ( itnext->mu * cellsize[2] ) / nom;

        #if DEBUG_MAIN_YEE >= 100
            if ( (itcur+ynext)->E != itcur->E )
                tout << "!!!! E[" << ynext << "]=" << (itcur+ynext)->E << " != E[" << pos << "]=" << itcur->E << "\n";
        #endif

        itnext->H[X] = Da * itnext->H[X] +
             Dbz*( (itnext+znext)->E[Y] - itnext->E[Y] )
            -Dby*( (itnext+ynext)->E[Z] - itnext->E[Z] )
            ;
        itnext->H[Y] = Da * itnext->H[Y] +
             Dbx*( (itnext+xnext)->E[Z] - itnext->E[Z] )
            -Dbz*( (itnext+znext)->E[X] - itnext->E[X] )
            ;
        itnext->H[Z] = Da * itnext->H[Z] +
             Dby*( (itnext+ynext)->E[X] - itnext->E[X] )
            -Dbx*( (itnext+xnext)->E[Y] - itnext->E[Y] )
            ;
    }

}

void CalcE( OctCell & simbox, int tnext, int tcur, int area, double dt ) {
    typename OctCell::CellMatrix CellMatrix;
    typename OctCell::IteratorType itnext = simbox.getIterator( tnext, area );
    typename OctCell::IteratorType itcur  = simbox.getIterator( tcur , area );

    Vec<double,3> cellsize( (SIMDIM>=1) ? simbox.cellsize[0] : INF ,
                            (SIMDIM>=2) ? simbox.cellsize[1] : INF ,
                            (SIMDIM>=3) ? simbox.cellsize[2] : INF );
    // ToDo: Check for Courant Criterium (now as we can calculate c_M :3) !!!

    for ( itnext=itnext.begin(); itnext != itnext.end(); ++itnext ) {
        itcur.icell = itnext.icell;

        VecI xprev(0); if (SIMDIM>=1) xprev[X]--;
        VecI yprev(0); if (SIMDIM>=2) yprev[Y]--;
        VecI zprev(0); if (SIMDIM>=3) zprev[Z]--;

		/* Now update all E components */
        double nom = (1.0 + 0.5*itnext->sigmaE*dt / itnext->epsilon);
        double Ca  = (1.0 - 0.5*itnext->sigmaE*dt / itnext->epsilon) / nom;
        if ( itnext->sigmaE > sqrt(DBL_MAX) )
            Ca = (1.0/itnext->sigmaE - 0.5 / itnext->epsilon) / (1.0/itnext->sigmaE + 0.5*dt / itnext->epsilon);
        double Cbx = dt / ( itnext->epsilon * cellsize[0] ) / nom;
        double Cby = dt / ( itnext->epsilon * cellsize[1] ) / nom;
        double Cbz = dt / ( itnext->epsilon * cellsize[2] ) / nom;

        itnext->E[X] = Ca * itnext->E[X] +
             Cby*( itnext->H[Z] - (itnext+yprev)->H[Z] )
            -Cbz*( itnext->H[Y] - (itnext+zprev)->H[Y] )
            ;
        itnext->E[Y] = Ca * itnext->E[Y] +
             Cbz*( itnext->H[X] - (itnext+zprev)->H[X] )
            -Cbx*( itnext->H[Z] - (itnext+xprev)->H[Z] )
            ;
        itnext->E[Z] = Ca * itnext->E[Z] +
             Cbx*( itnext->H[Y] - (itnext+xprev)->H[Y] )
            -Cby*( itnext->H[X] - (itnext+yprev)->H[X] )
            ;
    }
}

} // namespace YeeSolver
