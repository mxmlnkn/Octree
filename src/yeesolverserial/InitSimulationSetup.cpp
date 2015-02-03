/********** Result 009: absorbing Material on right side **********/
if ( SIMULATION_SETUP == 1 ) {
    if ( curPos[0] > 0.2 ) {
        /* For INF instead of 2e8 it diverges! -> characteristic *
         * wave length for exponential decay -> 0 for INF ?       */
        itm->sigmaE  = ABSORBER_STRENGTH;
        itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
    }
}
/*************** Result 014: broken total reflexion ***************/
/* two glass plates with small vacuum/air slit inbetween          */
/******************************************************************/
if ( SIMULATION_SETUP == 2 ) {
    const double n  = 1.33; // = sqrt( eps_r * mue_r )
    if ( curPos[X] < 0.333*tree.size[X]
    or   curPos[X] > 0.350*tree.size[X]  )
        itm->epsilon = EPS0 * n*n;
}
/********** Perfectly reflecting barrier with one slit ************/
if ( SIMULATION_SETUP == 3 ) {
    const double wy = LAMBDA;
    if ( curPos[X] > LAMBDA and curPos[X] < 2*LAMBDA )
    if ( curPos[Y] > 0.5 - wy/2 and curPos[Y] < 0.5 + wy/2 ) {
            itm->epsilon  = INF;//2*EPS0;
            //itm->mu       = INF;
    }
}
/************************* Circular Lense *************************/
if ( SIMULATION_SETUP == 4 ) {
    const double n  = 1.33; // = sqrt( eps_r * mue_r )
    if ( (curPos - M).norm() < R )
        itm->epsilon = EPS0 * n*n;
    if ( (curPos[X] - 0 < ABSORBING_BORDER_THICKNESS) or
         (tree.center[X] + tree.size[X]/2 - curPos[X] < ABSORBING_BORDER_THICKNESS) or
         (tree.center[Y] + tree.size[Y]/2 - curPos[Y] < ABSORBING_BORDER_THICKNESS) or
         (curPos[Y] - 0 < ABSORBING_BORDER_THICKNESS)
       )
    {
        /* For INF instead of 2e8 it diverges! -> characteristic *
         * wave length for exponential decay -> 0 for INF ?       */
        itm->sigmaE  = ABSORBER_STRENGTH;
        itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
    }
    /* waveguide absorbers for WAVE_SPAWN_SETUP 3 */
    if ( ( (curPos >= SPAWN_POS + VecD( 0.0, SPAWN_AREA_SIZE[1] ) ) and
           (curPos <  SPAWN_POS + VecD( 0.0, SPAWN_AREA_SIZE[1] ) + WAVE_GUIDE_SIZE ) ) or
         ( (curPos >= SPAWN_POS - VecD( 0.0, WAVE_GUIDE_SIZE[1] ) ) and
           (curPos < SPAWN_POS + VecD( WAVE_GUIDE_SIZE[0], 0.0 ) ) )
       )
   {
        itm->sigmaE  = ABSORBER_STRENGTH;
        itm->sigmaM  = ABSORBER_STRENGTH / SPEED_OF_LIGHT * MUE0/EPS0;
   }
}
/************************ Parabolic Mirror ************************/
if ( SIMULATION_SETUP == 5 ) {
    const double xcalc = pow( curPos[1] - MIRROR_CENTER[1], 2.0 ) / ( 4.0*FOCAL_LENGTH );
    if ( curPos[0] < xcalc + MIRROR_CENTER[0] ) {
        itm->sigmaE  = INF;
        itm->sigmaM  = INF / SPEED_OF_LIGHT * MUE0/EPS0;
   }
}
