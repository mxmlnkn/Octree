/* SIMULATION_SETUP 101 to 163 ***********************************/
bool leftabsorber   = curpos[X] - 0 < ABSORBING_BORDER_THICKNESS;
bool rightabsorber  = tree.center[X] + tree.size[X]/2 - curpos[X] < ABSORBING_BORDER_THICKNESS;
bool bottomabsorber = curpos[Y] - 0 < ABSORBING_BORDER_THICKNESS;
bool topabsorber    = tree.center[Y] + tree.size[Y]/2 - curpos[Y] < ABSORBING_BORDER_THICKNESS;

/** Result 009: absorbing Material on right side included ********/
if ( ( ABSORBER_SIDE[0] and leftabsorber )
or ( ( ABSORBER_SIDE[1] or /* Legacy */ CONTAINS(SIMULATION_SETUP,1) ) and rightabsorber )
or ( ABSORBER_SIDE[2] and bottomabsorber )
or ( ABSORBER_SIDE[3] and topabsorber ) ) {
    /* For INF instead of 2e8 it diverges! -> characteristic *
     * wave length for exponential decay -> 0 for INF ?      */
    itm->sigmaE  = ABSORBER_STRENGTH;
    itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
}


/*************** Result 014: broken total reflexion ***************/
/* two glass plates with small vacuum/air slit inbetween          */
/******************************************************************/
if ( CONTAINS(SIMULATION_SETUP,2) ) {
    const double n  = 1.33; // = sqrt( eps_r * mue_r )
    if ( curpos[X] < 0.333*tree.size[X]
    or   curpos[X] > 0.350*tree.size[X]  )
        itm->epsilon = EPS0 * n*n;
}
/********** Perfectly reflecting barrier with one slit ************/
if ( CONTAINS(SIMULATION_SETUP,3) ) {
    const double wy = LAMBDA;
    if ( curpos[X] > LAMBDA and curpos[X] < 2*LAMBDA )
    if ( curpos[Y] > 0.5 - wy/2 and curpos[Y] < 0.5 + wy/2 ) {
            itm->epsilon  = INF;//2*EPS0;
            //itm->mu       = INF;
    }
}
/************************* Circular Lense *************************/
if ( CONTAINS(SIMULATION_SETUP,4) ) {
    const double n  = 1.33; // = sqrt( eps_r * mue_r )
    if ( (curpos - SPHERICAL_LENSE_CENTER).norm() < SPHERICAL_LENSE_RADIUS )
        itm->epsilon = EPS0 * n*n;
    if ( (curpos[X] - 0 < ABSORBING_BORDER_THICKNESS) or
         (tree.center[X] + tree.size[X]/2 - curpos[X] < ABSORBING_BORDER_THICKNESS) or
         (tree.center[Y] + tree.size[Y]/2 - curpos[Y] < ABSORBING_BORDER_THICKNESS) or
         (curpos[Y] - 0 < ABSORBING_BORDER_THICKNESS)
       )
    {
        /* For INF instead of 2e8 it diverges! -> characteristic *
         * wave length for exponential decay -> 0 for INF ?       */
        itm->sigmaE  = ABSORBER_STRENGTH;
        itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
    }
    /* waveguide absorbers for WAVE_SPAWN_SETUP 3 */
    if ( ( (curpos >= SPAWN_POS + VecD( 0.0, SPAWN_AREA_SIZE[1] ) ) and
           (curpos <  SPAWN_POS + VecD( 0.0, SPAWN_AREA_SIZE[1] ) + WAVE_GUIDE_SIZE ) ) or
         ( (curpos >= SPAWN_POS - VecD( 0.0, WAVE_GUIDE_SIZE[1] ) ) and
           (curpos < SPAWN_POS + VecD( WAVE_GUIDE_SIZE[0], 0.0 ) ) )
       )
   {
        itm->sigmaE  = ABSORBER_STRENGTH;
        itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
   }
}
/************************ Parabolic Mirror ************************/
if ( CONTAINS(SIMULATION_SETUP,5) ) {
    const double xcalc = pow( curpos[1] - MIRROR_CENTER[1], 2.0 ) / ( 4.0*FOCAL_LENGTH );
    if ( curpos[0] < xcalc + MIRROR_CENTER[0] ) {
        itm->sigmaE  = INF;
        itm->sigmaM  = INF / SPEED_OF_LIGHT * MUE0/EPS0;
   }
}

/********* Using Lens to make point source to plane wave **********/
if ( CONTAINS(SIMULATION_SETUP,7) ) {
    /**************************************************************
     *           y                                                *
     * Absorber  ^    Air/Vacuum                                  *
     *         __|______         o...Source = M + R = (x0,y0)     *
     *       -   |  -   --                                        *
     *      - Air|   -    \                                       *
     *      |    o    | nG |                                      *
     *      -    |   -    /                                       *
     *       - __|__-___--                                        *
     *           |                                                *
     *           |---------------> x                              *
     **************************************************************/
    /* We want the center to be exactly the center there of the  *
     * wave which will be spawned! That's why we search for a    *
     * wanted pos and get back the center of the cell containing *
     * that position                                             */
    static bool firstCall = true;
    if ( firstCall ) {
        SPAWN_POS = VecD( ABSORBING_BORDER_THICKNESS + SPHERICAL_SCREEN_RADIUS,
                          ABSORBING_BORDER_THICKNESS + SPHERICAL_LENSE_CENTER[1]
                          + 0.84*SPHERICAL_LENSE_RADIUS );
        firstCall = false;
        /* if found cell, then send new corrected position to all other      *
         * processes. Problem here is, that the other processes don't know   *
         * at runtime, which process has had success / to receive from       */
        VecD * Ms = new VecD[combox.worldsize];
        VecD found = VecD(-8192);
        combox.findCell( SPAWN_POS, &found );
        MPI_Allgather( &found, sizeof(VecD), MPI_BYTE, Ms, sizeof(VecD), MPI_BYTE, MPI_COMM_WORLD );
        for ( int i = 0; i < combox.worldsize; ++i )
            if ( Ms[i] != VecD(-8192) ) {
                SPHERICAL_SCREEN_CENTER = Ms[i];
                SPAWN_POS               = Ms[i];
                break;
            }
        delete[] Ms;
        tout << "Corrected SPAWN_POS to: " << SPAWN_POS << "\n";
    }
    const double nVacuum = 1.0;
    const double nLense  = 1.33;
    const double e    = nVacuum / nLense; /* < 1 */
    const double & b  = ELLIPTIC_LENSE_SEMI_MINOR_AXIS;
    const double & R  = SPHERICAL_SCREEN_RADIUS;
    const double a    = b / sqrt( 1 - e*e );
    const double & x  = curpos[0];
    const double & y  = curpos[1];
    const double & x0 = SPHERICAL_SCREEN_CENTER[0];
    const double & y0 = SPHERICAL_SCREEN_CENTER[1];
    const double r    = (curpos - SPHERICAL_SCREEN_CENTER).norm();
    const double phi  = atan( (y-y0)/(x-x0) );
    double xEllipseRight = x0 + e*a + a*sqrt( 1 - pow( (y-y0)/b ,2 ) ); /* can be nan, which is good */
    double xEllipseLeft  = x0 + e*a - a*sqrt( 1 - pow( (y-y0)/b ,2 ) );
    bool isLense      = r >= R and x <= xEllipseRight and x >= xEllipseLeft;
    bool isSpawnGuard = r >= R and r <= R + ABSORBING_BORDER_THICKNESS
                        and ( fabs(phi) > M_PI/6 or x < x0 );
    /* comparison with NaN above will be always false! */
    if ( isSpawnGuard ) {
       itm->sigmaE  = ABSORBER_STRENGTH;
       itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
    }
    if ( isLense ) {
       itm->epsilon = EPS0 * nLense*nLense;
    }
}
