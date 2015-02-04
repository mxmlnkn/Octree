    /* NUMBER_OF_PARTICLES_PER_CELL, BOUNDARY_CONDITION, SPECIES, PNG_INTERVAL  -> Watch out for dependent Variables in Parameters.cpp!!! :  NUMBER_OF_EONS_PER_CELL, NUMBER_OF_IONS_PER_CELL       */

    while ( true ) {
        static struct option long_options[] = {
            {"timesteps"       , required_argument, 0, 't'},
            {"init-refinement" , required_argument, 0, 'i'},
            {"max-refinement"  , required_argument, 0, 'm'},
            {"number-of-cells" , required_argument, 0, 'n'},
            {"octree-setup"    , required_argument, 0, 'o'},
            {"simulation-setup", required_argument, 0, 's'},
            {"png-interval"    , required_argument, 0, 'p'},
            {"wave-source"     , required_argument, 0, 'w'},
            {"absorber"        , required_argument, 0, 'a'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "t:i:m:n:o:s:p:w:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
            case 't':
                NUMBER_OF_STEPS = atoi(optarg);
                break;
            case 'i':
                INITIAL_OCTREE_REFINEMENT = atoi(optarg);
                break;
            case 'm':
                MAX_OCTREE_REFINEMENT = atoi(optarg);
                break;
            case 'n':
                for (int i=0; i<SIMDIM; i++) {
                    assert( argv[optind-1+i][0] != '-' );
                    NUMBER_OF_CELLS[i]  = atoi(argv[optind-1+i]);
                }
                optind += SIMDIM-1; // extra arguments taken
                NUMBER_OF_CELLS_X   = NUMBER_OF_CELLS[0];
                if (SIMDIM > 1)
                    NUMBER_OF_CELLS_Y = NUMBER_OF_CELLS[1];
                else
                    NUMBER_OF_CELLS_Y = 1;
                if (SIMDIM > 2)
                    NUMBER_OF_CELLS_Z = NUMBER_OF_CELLS[2];
                else
                    NUMBER_OF_CELLS_Z = 1;
                NUMBER_OF_PARTICLES = NUMBER_OF_PARTICLES_PER_CELL *
                    NUMBER_OF_CELLS_X * NUMBER_OF_CELLS_Y * NUMBER_OF_CELLS_Z;
                SIM_SIZE = Vec<double,SIMDIM>( NUMBER_OF_CELLS ) * CELL_SIZE;
                SPAWN_POS              = SIM_SIZE / 2.0;
                SPHERICAL_LENSE_CENTER = VecD( 0.5*SIM_SIZE[0], 0.5*SIM_SIZE[1] );    // center of circle
                SPHERICAL_LENSE_RADIUS = 0.4*SIM_SIZE[1];
                MIRROR_WIDTH           = SIM_SIZE[1] - 2*ABSORBING_BORDER_THICKNESS;
                MIRROR_LENGTH          = SIM_SIZE[0]/8;
                break;
            case 'o':
                OCTREE_SETUP.push_back( atoi(optarg) );
                break;
            case 's':
                SIMULATION_SETUP.push_back( atoi(optarg) );
                /* 100 + direction:                                           *
                 * direction = 000000b         => 111111b = 63 activates      *
                 *             |||||+- left       absorber everywhere         *
                 *             ||||+-- right                                  *
                 *             |||+--- bottom                                 *
                 *             ||+---- top                                    *
                 *             |+----- back                                   *
                 *             +------ front                                  */
                 {int tsimsetup = atoi(optarg);
                if ( tsimsetup >= 101 and tsimsetup <= 163 ) {
                    for ( int bit = 0; bit < 2*SIMDIM; ++bit ) {
                        if ( (tsimsetup-100) & (1 << bit) )
                            ABSORBER_SIDE[bit] = true;
                    }
                 }}
                break;
            case 'p':
                PNG_INTERVAL = atoi(optarg);
                break;
            case 'w':
                WAVE_SPAWN_SETUP.push_back( atoi(optarg) );
                if ( CONTAINS(WAVE_SPAWN_SETUP, 3) ) {
                if ( SIMDIM == 2 and argv[optind+0][0] != '-' and argv[optind+1][0] != '-' )
                {
                    SPAWN_POS = Vec<double,SIMDIM>( atoi(argv[optind+0]), atoi(argv[optind+1]) );
                    optind += 2;
                }
                if ( SIMDIM == 2 and argv[optind+0][0] != '-' and argv[optind+1][0] != '-' )
                {
                    SPAWN_AREA_SIZE = Vec<double,SIMDIM>( atoi(argv[optind+0]), atoi(argv[optind+1]) );
                    optind += 2;
                }
                }
                break;
            case 'a':
                ABSORBER_STRENGTH = atoi(optarg);
                if ( argv[optind-2][0] != '-' ) {
                    ABSORBING_BORDER_THICKNESS = atoi(argv[optind-2]);
                    optind += 1;
                }
                break;
            default:
                abort();
        }
    }

if ( CONTAINS(WAVE_SPAWN_SETUP,7) or CONTAINS(OCTREE_SETUP,9) ) {
    SPAWN_POS = VecD( ABSORBING_BORDER_THICKNESS + SPHERICAL_SCREEN_RADIUS,
                      ABSORBING_BORDER_THICKNESS + SPHERICAL_LENSE_CENTER[1]
                      + 0.7*SPHERICAL_LENSE_RADIUS );
    SPHERICAL_SCREEN_CENTER = SPAWN_POS;
}
