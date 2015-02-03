/********************** Gaussian Wave Pulse ***********************/
if ( WAVE_SPAWN_SETUP % 10 == 5 ) {
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
    if ( WAVE_SPAWN_SETUP == 51 ) {
        itm->H[Y] = -1.0 / ( cM * MUE0 ) *
            TIME_SPAWN_FUNCTIONS::sinewave2d(
            T, 0.0, kx, curPos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
            CELL_SIZE[0]), ky, curPos[1]-SPAWN_POS[1] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
        itm->E[Z] = TIME_SPAWN_FUNCTIONS::sinewave2d(
            T, 0.5*DELTA_T, kx, curPos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] -
            CELL_SIZE[0]), ky, curPos[1]-SPAWN_POS[1] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
    }
    if ( WAVE_SPAWN_SETUP == 52 ) {
        /* B = µH = E/c */
        itm->H[Y] = - TIME_SPAWN_FUNCTIONS::GaussianBeam(
            curPos[0] - SPAWN_POS[0], curPos[1] - SPAWN_POS[1], 0,
            SPAWN_AREA_SIZE[1], LAMBDA, T, 1 ) / ( cM * MUE0 ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] )*
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
        /* E is half a time step later ! */
        itm->E[Z] = TIME_SPAWN_FUNCTIONS::GaussianBeam( curPos[0] - SPAWN_POS[0],
            curPos[1] - SPAWN_POS[1], 0.5, SPAWN_AREA_SIZE[1], LAMBDA, T, 1 ) *
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[0], SPAWN_POS[0], SPAWN_AREA_SIZE[0] )*
            TIME_SPAWN_FUNCTIONS::GaussNotNormed( curPos[1], SPAWN_POS[1], SPAWN_AREA_SIZE[1] );
    }
    /* Cut Off spawned wave to not dirty the whole simulation space */
    if ( std::abs(itm->E[Z]) < 0.001) {
        itm->E[Z] = 0;
        itm->H[Y] = 0;
    }
}
/********* Using Lens to make point source to plane wave **********/
if ( WAVE_SPAWN_SETUP % 10 == 7 ) {
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
         firstCall = false;
         R = 2*LAMBDA;
         VecD guess = VecD(ABSORBING_BORDER_THICKNESS + R, SIM_SIZE[1]/2 );
         /* if found cell, then send new corrected position to all other      *
          * processes. Problem here is, that the other processes don't know   *
          * at runtime, which process has had success / to receive from       */
         VecD * Ms = new VecD[combox.worldsize];
         M = VecD(-8192);
         combox.findCell( guess, &M );
         MPI_Allgather( &M, sizeof(VecD), MPI_BYTE, Ms, sizeof(VecD), MPI_BYTE, MPI_COMM_WORLD );
         for ( int i = 0; i < combox.worldsize; ++i )
             if ( Ms[i] != VecD(-8192) ) {
                 SPAWN_POS = Ms[i];
                 M         = Ms[i];
                 break;
             }
         delete[] Ms;
         tout << "Corrected SPAWN_POS to: " << SPAWN_POS << "\n";
     }
     const double nVacuum = 1.0;
     const double nLense  = 1.33;
     const double e    = nVacuum / nLense;
     const double & b  = 4*R;
     const double a    = b / sqrt( 1 - e*e );
     const double & x  = curPos[0];
     const double & y  = curPos[1];
     const double & x0 = M[0];
     const double & y0 = M[1];
     const double r    = (curPos-M).norm();
     const double phi  = atan( (y-y0)/(x-x0) );
     double xLeftCirc  = fabs(y-y0) < R ? x0 - R*sqrt( 1 - pow( (y-y0)/R, 2 ) ) : x0;
     double xRightCirc = fabs(y-y0) < R ? x0 + R*sqrt( 1 - pow( (y-y0)/R, 2 ) ) : x0;
     double xEllipseRight = fabs(y-y0) < b ? x0 + e*a + a*sqrt( 1 - pow( (y-y0)/b ,2 ) ) : x0;
     double xEllipseLeft  = fabs(y-y0) < b ? x0 + e*a - a*sqrt( 1 - pow( (y-y0)/b ,2 ) ) : x0;
     bool isAbsorber   = x < x0 and x < xLeftCirc;
     bool isLense      = x >= xRightCirc and x <= xEllipseRight and x >= xEllipseLeft;
     bool isSpawnGuard = r >= R and r <= R + 2*ABSORBING_BORDER_THICKNESS and fabs(phi) > M_PI/6;
     /* comparison with NaN above will be always false! */
     if ( isAbsorber or isSpawnGuard ) {
        itm->sigmaE  = ABSORBER_STRENGTH;
        itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
     }
     if ( isLense ) {
        itm->epsilon = EPS0 * nLense*nLense;
     }
     /* Set AbsorberBorders */
     if ( (curPos[X] - 0 < ABSORBING_BORDER_THICKNESS) or
          (tree.center[X] + tree.size[X]/2 - curPos[X] < ABSORBING_BORDER_THICKNESS) or
          (tree.center[Y] + tree.size[Y]/2 - curPos[Y] < ABSORBING_BORDER_THICKNESS) or
          (curPos[Y] - 0 < ABSORBING_BORDER_THICKNESS)
        )
     {
         itm->sigmaE = ABSORBER_STRENGTH;
         itm->sigmaM = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
     }
}
