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
