if ( WAVE_SPAWN_SETUP == 1 or WAVE_SPAWN_SETUP == 2 or WAVE_SPAWN_SETUP == 7 ) {
    /* Function Generator on Cell in the Center */
    #if 1 == 0
    OctreeType::Node * node = tree.FindLeafContainingPos( SPAWN_POS );
    if ( ((OctreeCommType::CommData*)node->data[OctreeCommType::COMM_HEADER_INDEX])->rank == combox.rank ) {
        OctCell & simbox = *((OctCell*)node->data[OctreeCommType::CELL_DATA_INDEX]);
        #if DEBUG_MAIN_YEE >= 100
            tout << "Find position " << SPAWN_POS << " in tree sized " << tree.size << " positioned at center " << tree.center << " returned OctreeCell at " << node->center << " => Searching for position in simbox of that node with abspos " << simbox.abspos << ", localcells " << simbox.localcells << " and cellsize " << simbox.cellsize << "\n";
        #endif
        VecI targetIndex = simbox.findCellContaining( SPAWN_POS );
        simbox.t[0]->cells[targetIndex].E[Z] = 20*t_spawn_func( t * DELTA_T_SI );
        #if DEBUG_MAIN_YEE >=100
            tout << "Write to source in cell " << targetIndex << " in node at " << node->center << "\n";
        #endif
    }
    #endif 
    YeeCell * cell = combox.findCell( SPAWN_POS );
    if ( cell != NULL )
        cell->E[2] = 20*t_spawn_func( t * DELTA_T_SI );
}
if ( WAVE_SPAWN_SETUP == 2 ) {
    /* shield function generator in one direction */
    OctreeType::Node * node = tree.FindLeafContainingPos( SPAWN_POS );
    if ( ((OctreeCommType::CommData*)node->data[OctreeCommType::COMM_HEADER_INDEX])->rank == combox.rank ) {
        OctCell & simbox = *((OctCell*)node->data[OctreeCommType::CELL_DATA_INDEX]);
        VecI targetIndex = simbox.findCellContaining( SPAWN_POS );
        /*simbox.t[0]->cells[targetIndex + VecI(8,-4)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,-3)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,-2)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,-1)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,0)   ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,1)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,2)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,3)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(8,4)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,-4)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,-3)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,-2)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,-1)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,0)   ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,1)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,2)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,3)  ].E[Z] = 0;
        simbox.t[0]->cells[targetIndex + VecI(9,4)  ].E[Z] = 0;*/
    }
}
if ( WAVE_SPAWN_SETUP == 3 ) {
    /**********************************************************************
     * Sine plane Wave going to Direction alpha and beginning line going  *
     * through pos0  y                                                    *
     *               ^                                                    *
     *               | \     e.g. p0 ( line includes p0! )                *
     *               |  --  /                                             *
     *               |    \     alpha                                     *
     *               |     --  /                                          *
     *               |_______\__________ x                                *
     **********************************************************************/
    //double const n = 1.0;
    double alpha   = 0 / 360. * 2*M_PI; //0.9*asin(1./n); // radian
    //VecI pos0(0); pos0[X] = 6*lambda;
    double T0x     = LAMBDA / SPEED_OF_LIGHT;
    double T0y     = LAMBDA / SPEED_OF_LIGHT;
    double kx      = 2*M_PI/LAMBDA * cos(alpha);
    double ky      = 2*M_PI/LAMBDA * sin(alpha);
    double width   = SPAWN_AREA_SIZE[1];
    VecD pos       = SPAWN_POS;// + SPAWN_AREA_SIZE - 1;
    assert( SPAWN_POS[1]+SPAWN_AREA_SIZE[1] < tree.center[1] + tree.size[1] );
    
    VecD mincellsize = CELL_SIZE / pow( 2.0, combox.maxLevel - combox.minLevel );
    for (double x = SPAWN_POS[0]; x < SPAWN_POS[0] + SPAWN_AREA_SIZE[0]; x += mincellsize[0] )
    for (double y = SPAWN_POS[1]; y < SPAWN_POS[1] + SPAWN_AREA_SIZE[1]; y += mincellsize[1] ) {
        pos[0] = x; pos[1] = y;
        VecD cellpos;
        YeeCell * cell = combox.findCell( pos, &cellpos );
        if ( cell == NULL )
            continue;
        cell->E[2] = TIME_SPAWN_FUNCTIONS::sinewave2d( T0x, t * DELTA_T, kx,
                        pos[0] - (SPAWN_POS[0] + SPAWN_AREA_SIZE[0] - 
                        CELL_SIZE[0]), ky, pos[1]-SPAWN_POS[1] );
    }
}
