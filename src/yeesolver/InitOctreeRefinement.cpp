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
    const VecD   & M = SPHERICAL_LENSE_CENTER;
    const double & R = SPHERICAL_LENSE_RADIUS;
    assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
    for ( int lvl=INITIAL_OCTREE_REFINEMENT; lvl<MAX_OCTREE_REFINEMENT; lvl++) {
        std::stack<OctreeType::Node*> torefine;
        for ( OctreeType::iterator it = tree.begin(); it != tree.end(); ++it )
            if ( it->IsLeaf() and ( M - it->center*tree.size ).norm() <= R ) {
                bool oneneighboroutside = false;
                for ( int i = 0; i < pow(3,SIMDIM); ++i ) {
                    OctreeType::Node * neighbor =
                        it->getNeighbor( getDirectionVector<SIMDIM>(i), VecI(0) );
                    if ( neighbor == NULL ) {
                        oneneighboroutside = true;
                        continue;
                    }
                    if ( ( M - neighbor->center*tree.size ).norm() > R ) {
                        oneneighboroutside = true;
                        break;
                    }
                }
                if ( oneneighboroutside )
                    torefine.push( &(*it) );
            }
        while ( not torefine.empty() ) {
            torefine.top()->GrowUp();
            torefine.pop();
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
            for ( int i=0; i < pow(3,SIMDIM); ++i ) {
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
