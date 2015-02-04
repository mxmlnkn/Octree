/********************** Gaussian Wave Pulse ***********************/
if ( CONTAINS(WAVE_SPAWN_SETUP,51) or CONTAINS(WAVE_SPAWN_SETUP,52) ) {
    /**************************************************************
     * Sine plane Wave going to Direction alpha and beginning     *
     * line going through pos0  y                                 *
     *               ^                                            *
     *               | \     e.g. p0 ( line includes p0! )        *
     *               |  --  /                                     *
     *               |    \     alpha                             *
     *               |     --  /                                  *
     *               |_______\__________ x                        *
     **************************************************************/
    double alpha    = 0 / 360. * 2*M_PI; //0.9*asin(1./n); // radian
    double NUMERICAL_SLOWING = 0.996797;
    double n                 = 1.0;
    double cM                = SPEED_OF_LIGHT*NUMERICAL_SLOWING/n;
    double T        = LAMBDA / cM;
    double kx       = 2*M_PI/LAMBDA * cos(alpha);
    double ky       = 2*M_PI/LAMBDA * sin(alpha);
    SPAWN_POS       = VecD( 0.1*SIM_SIZE[0], 0.75*SIM_SIZE[1] );
    SPAWN_AREA_SIZE = VecD( SIM_SIZE[0]/32, SIM_SIZE[1]/64 );
    assert( SPAWN_POS[1]+SPAWN_AREA_SIZE[1] < tree.center[1] + tree.size[1] );
    static bool firstIteration = true;
    if (firstIteration == true) {
        firstIteration = false;
        tout << "Spawn Position: " << SPAWN_POS << ", Size of Pulse: " << SPAWN_AREA_SIZE << "\n";
    }
    if ( CONTAINS(WAVE_SPAWN_SETUP,51) ) {
        itm->H[Y] = -1.0 / ( cM * MUE0 ) *
            TIME_SPAWN_FUNCTIONS::sinewave2d(
            T, 0.0, kx, curpos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
            CELL_SIZE[0]), ky, curpos[1]-SPAWN_POS[1] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
        itm->E[Z] = TIME_SPAWN_FUNCTIONS::sinewave2d(
            T, 0.5*DELTA_T, kx, curpos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
            CELL_SIZE[0]), ky, curpos[1]-SPAWN_POS[1] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
    }
    if ( CONTAINS(WAVE_SPAWN_SETUP,52) ) {
        /* B = µH = E/c */
        itm->H[Y] = - TIME_SPAWN_FUNCTIONS::GaussianBeam(
            curpos[0] - SPAWN_POS[0], curpos[1] - SPAWN_POS[1], 0,
            SPAWN_AREA_SIZE[1], LAMBDA, T, 1 ) / ( cM * MUE0 ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] )*
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
        /* E is half a time step later ! */
        itm->E[Z] = TIME_SPAWN_FUNCTIONS::GaussianBeam( curpos[0] - SPAWN_POS[0],
            curpos[1] - SPAWN_POS[1], 0.5, SPAWN_AREA_SIZE[1], LAMBDA, T, 1 ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] )*
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curpos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
    }
    /* Cut Off spawned wave to not dirty the whole simulation space */
    if ( std::abs(itm->E[Z]) < 0.001) {
        itm->E[Z] = 0;
        itm->H[Y] = 0;
    }
}
/********* Using Lens to make point source to plane wave **********/
if ( CONTAINS(WAVE_SPAWN_SETUP,7) ) {
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
                          + 0.7*SPHERICAL_LENSE_RADIUS );
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
