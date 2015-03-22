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
            {"help"            , no_argument      , 0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "t:i:m:n:o:s:p:w:a:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
            case 'h':
                std::cout << "-t --timestep\n";
                std::cout << "   Set number of timesteps to calculate\n";
                std::cout << "-i --init-refinement\n";
                std::cout << "   Set initial homogenous refinement depth for Octree or Quadtree. 0 means no refinement. Only used of '-o 1' specified.\n";
                std::cout << "-m --mac-refinement\n";
                std::cout << "   Set maximum Refinement. If smaller or equal than number given to -i, then no complex octree refinement is being applied. See -o to specify what to refine further than homogenous\n";
                std::cout << "-n --number-of-cells\n";
                std::cout << "   Specify number of cells. Takes 1,2 or 3 arguments depending on SIMDIM in Parameters.cpp\n";
                std::cout << "-o --octree-setup\n";
                std::cout << "   Takes one argument, but this option can be specified arbitrarily often to combine setups\n";
                std::cout << "   1 - Refine homogenously up to level specified by -i\n";
                std::cout << "   6 - Refine sphere and wave spawning area\n";
                std::cout << "   7 - Refine cells at one point\n";
                std::cout << "   8 - Refine wave spawning area\n";
                std::cout << "   9 - Refine cells containing elliptic lense\n";
                std::cout << "-s simulation-setup\n";
                std::cout << "   Takes one argument, but this option can be specified arbitrarily often to combine setups\n";
                std::cout << "   100 + direction:                                     \n";
                std::cout << "   direction = 000000b         => 111111b = 63 activates\n";
                std::cout << "               |||||+- left       absorber everywhere   \n";
                std::cout << "               ||||+-- right                            \n";
                std::cout << "               |||+--- bottom                           \n";
                std::cout << "               ||+---- top                              \n";
                std::cout << "               |+----- back                             \n";
                std::cout << "               +------ front                            \n";
                std::cout << "   2 - Result 014: broken total reflexion\n";
                std::cout << "   3 - Perfectly reflecting barrier with one slit\n";
                std::cout << "   4 - Circular Lense\n";
                std::cout << "   5 - Parabolic Mirror\n";
                std::cout << "   7 - Elliptic Lense\n";
                std::cout << "-p png-interval\n";
                std::cout << "   Every timestep % <argument> a png is being written\n";
                std::cout << "-w --wave-source\n";
                std::cout << "    1 - Sine Wave over all timesteps at SPAWN_POS with wavelength LAMBDA_SI\n";
                std::cout << "    2 - old shield function generator in one direction around SPAWN_POS consisting of a reflecting wall\n";
                std::cout << "    3 - Plane wave slanted by ALPHA\n";
                std::cout << "   51 - Initializes a gaussian wave packet over many cells only at the first time step\n";
                std::cout << "   52 - similar to 51\n";
                std::cout << "-a --absorber\n";
                std::cout << "   Specify absorber thickness in number of cells\n";
                std::cout << "-h --help\n";
                std::cout << "   Print this message\n";
                exit(0);
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
                if ( CONTAINS(WAVE_SPAWN_SETUP, 3) and optind+1 < argc ) {
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
                ABSORBER_STRENGTH = atof(optarg);
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
