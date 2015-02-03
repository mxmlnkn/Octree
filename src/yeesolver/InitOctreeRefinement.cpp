    /**************************************************************************/
    /* (1) Setup Octree Refinement ********************************************/
    /**************************************************************************/

    /********* refine all cells to initial homogenous min-Refinement **********/
    if ( OCTREE_SETUP == 7 or OCTREE_SETUP == 6 or OCTREE_SETUP == 9 ) {
        for ( int lvl=0; lvl<INITIAL_OCTREE_REFINEMENT; lvl++) {
            for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it )
                if ( it->IsLeaf() and it->getLevel()==lvl ) it->GrowUp();
        }
    }
    /*********************** Refine certain boundaries ************************/
    if ( OCTREE_SETUP == 6 ) {
        assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
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
    if ( OCTREE_SETUP == 7 ) {
        tree.root->GrowUp();
        tree.FindLeafContainingPos( 0.75*globSize )->GrowUp();
    }
    /******* Grow up Spawning Area ********/
    if ( (OCTREE_SETUP == 8 or OCTREE_SETUP == 6) and false ) {
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
    if ( OCTREE_SETUP == 7 ) {
        assert( MAX_OCTREE_REFINEMENT >= INITIAL_OCTREE_REFINEMENT );
        for ( int lvl=INITIAL_OCTREE_REFINEMENT; lvl<MAX_OCTREE_REFINEMENT; lvl++)
            for ( OctreeType::iterator it=tree.begin(); it != tree.end(); ++it ) {
                bool insideLense = false;
                if ( it->IsLeaf() and it->getLevel()==lvl and insideLense )
                    it->GrowUp();
            }
    }
    
    tout << "Tree-Integrity: " << tree.CheckIntegrity() << "\n\n";