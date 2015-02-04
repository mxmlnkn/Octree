/**************************************************************************/
/* (1) Setup Octree Refinement ********************************************/
/**************************************************************************/

/********* refine all cells to initial homogenous min-Refinement **********/
if ( CONTAINS(OCTREE_SETUP, 1) ) {
    for ( int lvl=0; lvl<INITIAL_OCTREE_REFINEMENT; lvl++) {
        for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it )
            if ( it->IsLeaf() and it->getLevel()==lvl ) it->GrowUp();
    }
}
/*********************** Refine sphere boundaries *************************/
if ( CONTAINS(OCTREE_SETUP, 6) ) {
    assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
    const VecD   & M = SPHERICAL_LENSE_CENTER;
    const double & R = SPHERICAL_LENSE_RADIUS;
    for ( int lvl=INITIAL_OCTREE_REFINEMENT; lvl<MAX_OCTREE_REFINEMENT; lvl++) {
        /* Get all circle angles, where it intersects with a cell border */
        std::list<double> lphi;
        std::list<double>::iterator it;
        VecD cellsize = globSize / pow(2,lvl);
        double linexmin = ceil ( (M[0]-R)/cellsize[0] ) * cellsize[0];
        double linexmax = floor( (M[0]+R)/cellsize[0] ) * cellsize[0];
        for ( double linex=linexmin; linex<=linexmax; linex += cellsize[0] ) {
            /* acos in [0,2*pi] */
            double phi = acos( (linex-M[0])/R );
            lphi.push_back( phi );
            /* also add value mirrored at y-axis to stack */
            lphi.push_back( 2*M_PI-phi );
        }
        double lineymin = ceil ( (M[1]-R)/cellsize[1] ) * cellsize[1];
        double lineymax = floor( (M[1]+R)/cellsize[1] ) * cellsize[1];
        for ( double liney=lineymin; liney<=lineymax; liney += cellsize[1] ) {
            /* asin in [-pi,pi] */
            double phi = asin( (liney-M[1])/R );
            lphi.push_back( phi < 0 ? 2*M_PI+phi : phi );
            /* also add value mirrored at x-axis to stack */
            lphi.push_back( M_PI-phi );
        }
        lphi.sort();

        /* Echo all found angles */
        #if DEBUG_MAIN_YEE >= 100
            tout << "Angle list contains:";
            for (it=lphi.begin(); it!=lphi.end(); ++it)
                tout << ' ' << *it;
            tout << '\n';
        #endif

        /* Grow up all cells, with which the circle intersects. Find them by  *
         * using an angle between to successive circle intersection angles    */
        for (it=lphi.begin(); it!=lphi.end(); ++it) {
            VecD pos(0);
            std::list<double>::iterator itnext = it;
            double phi;
            if ( ++itnext == lphi.end() ) {
                itnext = lphi.begin();
                phi = 0.5 * (2*M_PI + *itnext + *it);
            } else
                phi = 0.5 * (*itnext + *it);
            pos[0] = M[0] + R*cos(phi);
            pos[1] = M[1] + R*sin(phi);
            OctreeType::Node * node = tree.FindLeafContainingPos(pos);
            if ( node->getLevel() == lvl )
                node->GrowUp();
        }
    }
}
/******* Just grow up one point in upper left area (minimal setup) ********/
if ( CONTAINS(OCTREE_SETUP,7) ) {
    tree.root->GrowUp();
    tree.FindLeafContainingPos( 0.75*globSize )->GrowUp();
}
/******* Grow up Spawning Area ********/
if ( (CONTAINS(OCTREE_SETUP, 6) or CONTAINS(OCTREE_SETUP, 8)) and false ) {
    double minCellSizeY = tree.size[1] / pow(2.,MAX_OCTREE_REFINEMENT);
    for ( VecD pos = SPAWN_POS; pos < SPAWN_POS + VecD(SPAWN_AREA_SIZE); pos[1] += minCellSizeY ) {
        OctreeType::Node * curNode = tree.FindLeafContainingPos(pos);
        if ( curNode->getLevel() < MAX_OCTREE_REFINEMENT ) {
            curNode->GrowUp();
            pos[1] -= minCellSizeY;
        }
    }
}
/**************************************************************************/
if ( CONTAINS(OCTREE_SETUP,9) ) {
    assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
    for ( int lvl=INITIAL_OCTREE_REFINEMENT; lvl<MAX_OCTREE_REFINEMENT; lvl++)
        for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it ) {
            if ( not it->IsLeaf() or it->getLevel() != lvl )
                continue;
            int cellsinside  = 0;
            int cellsoutside = 0;
            /* if the cell itself or one of it's neighbors is inside lense,   *
             * then refine the cell                                           */
            for ( int i=0; i < pow(3,SIMDIM) - 1; ++i ) {
                VecI dir = getDirectionVector<SIMDIM>(i);
                typename OctreeType::Node * neighbor = it->getNeighbor(dir,VecI(0));
                if ( neighbor == NULL )
                    continue;
                VecD curpos = tree.toGlobalCoords( neighbor->center );
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
                double xRightCirc = fabs(y-y0) <= R ? x0 + R*sqrt( 1 - pow( (y-y0)/R, 2 ) ) : x0;
                double xEllipseRight = fabs(y-y0) <= b ? x0 + e*a + a*sqrt( 1 - pow( (y-y0)/b ,2 ) ) : x0;
                double xEllipseLeft  = fabs(y-y0) <= b ? x0 + e*a - a*sqrt( 1 - pow( (y-y0)/b ,2 ) ) : x0;
                bool insideLense = x >= xRightCirc and x <= xEllipseRight and x >= xEllipseLeft;

                if ( insideLense )
                    cellsinside++;
                else
                    cellsoutside++;
            }
            if ( cellsinside > 0 and cellsoutside > 0 )
                it->GrowUp();
        }
}

tout << "Tree-Integrity: " << tree.CheckIntegrity() << "\n\n";
