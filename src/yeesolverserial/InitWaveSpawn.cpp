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

