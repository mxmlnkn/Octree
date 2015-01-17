#pragma once

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
const int OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::COMM_HEADER_INDEX = 0;

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
const int OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::CELL_DATA_INDEX = 1;

/********************************* Constructor ********************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::OctreeCommunicator
( T_OCTREE & tree, const VecI periodic )
: rank(0), worldsize(1), tree(tree), periodic(periodic), 
  minLevel(0), maxLevel(0), comDataPtr(NULL),
  NLeaves(0), cellsPerOctreeCell(0), guardsize(0), timestepbuffers(0),
  neighbors(NULL), recvrequests(NULL), sendrequests(NULL)
{
    /* Initialize MPI */
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &processor_name_length);
    
    this->neighbors    = new NeighborData[worldsize];
#if 1==0
    this->recvrequests = new MPI_Request[nneighbors+1];
    this->sendrequests = new MPI_Request[nneighbors+1];
    this->sendmatrices = new SimBox::CellMatrix[nneighbors+1];
    this->recvmatrices = new SimBox::CellMatrix[nneighbors+1];
    if ( neighbors == NULL or recvrequests == NULL or sendrequests == NULL )
        terr << "Couldn't allocate Memory!";
    for (int i=0; i<=nneighbors; i++) {
        sendrequests[i] = MPI_REQUEST_NULL;
        recvrequests[i] = MPI_REQUEST_NULL;
    }
#endif
}

/********************************** Destructor ********************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::~OctreeCommunicator() {
    for ( typename T_OCTREE::iterator it=tree.begin(); it!=tree.end(); ++it ) {
        if ( it->data.size() >= 2 ) {
            int curRank = ((CommData*)it->data[COMM_HEADER_INDEX])->rank;
            assert( curRank == this->rank );
            delete ((OctCell*)it->data[CELL_DATA_INDEX]);
        }
    }

    if ( comDataPtr != NULL )
        delete[] comDataPtr;

    delete[] this->neighbors;
#if 1==0
    delete[] this->sendrequests;
    delete[] this->recvrequests;
    delete[] this->sendmatrices;
    delete[] this->recvmatrices;
#endif
}

/********************************* initCommData *******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::initCommData
( VecI cellsPerOctreeCell, int guardsize, int timesteps ) {
    /* number of cells rises exponentially. Noone will use 2^255 Cells! */
    this->minLevel = 255;
    this->maxLevel = 0;
    for ( typename T_OCTREE::iterator it = tree.begin(); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) {
            int curLevel = it->getLevel();
            assert( curLevel <= 255 );
            if ( minLevel > curLevel )
                minLevel = curLevel;
            assert( curLevel >= 0 );
            if ( maxLevel < curLevel )
                maxLevel = curLevel;
        }
    
    this->cellsPerOctreeCell = cellsPerOctreeCell;
    this->guardsize = guardsize;
    this->timestepbuffers = timesteps;
    /* allocate data which stores assigned ranks and other communication      *
     * information to which pointers will given to octree. And default it to  *
     * the last rank                                                          */
    this->NLeaves = tree.root->countLeaves();
    this->comDataPtr = new CommData[NLeaves];
    /* Insert Pointers to array elements of allocated chunk into octree and   *
     * assign weighting and caclulate total weighting and optimalCosts        */
    int dataInserted = 0;
    this->totalCosts = 0;
    for ( typename T_OCTREE::iterator it = tree.begin(); it != tree.end(); ++it ) {
        if ( it->IsLeaf() ) {
            assert( dataInserted < NLeaves );
            assert( it->data.empty() );
            it->data.push_back( &(comDataPtr[dataInserted]) );
            ++dataInserted;

            assert( COMM_HEADER_INDEX == 0 );
            ((CommData*)it->data[COMM_HEADER_INDEX])->weighting = 1. / it->size.min();

            this->totalCosts += ((CommData*)it->data[COMM_HEADER_INDEX])->weighting;
        }
    }
    this->optimalCosts = totalCosts / double(worldsize);
    tout << "Total Costs: " << totalCosts << " => Optimal Costs: " << optimalCosts << std::endl;
}

/******************************* distributeCells ******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::distributeCells(int ordering) {
    this->distOrdering = ordering;
    /* Assign cells to all the processes and allocate memory */
    double cumulativeCosts = 0;
    int curRank = 0;
    for ( typename T_OCTREE::iterator it=tree.begin(ordering); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) {
            cumulativeCosts += 1. / it->getSize().min();
            if ( cumulativeCosts >= optimalCosts and curRank != worldsize-1) {
                cumulativeCosts = 1. / it->getSize().min();
                curRank++;
            }
            ((CommData*)it->data[COMM_HEADER_INDEX])->rank = curRank;

            /* Allocate memory for celldata only if cell belongs to us */
            if ( curRank == this->rank ) {
                NOwnLeaves++;
                assert(it->data.size() == CELL_DATA_INDEX);

                it->data.push_back(
                    new SimulationBox::SimulationBox<SIMDIM,T_CELLTYPE> (
                        tree.toGlobalCoords(it->center - 0.5*it->size), /* abspos, lower left */
                        cellsPerOctreeCell, /* localcells */
                        tree.size / pow( 2, it->getLevel() ) / cellsPerOctreeCell, /* cellsize */
                        this->guardsize, this->timestepbuffers )
                );
            }
        }

    /* Count cells to transmit per neighbor process */
    for ( typename T_OCTREE::iterator it=tree.begin(ordering); it != tree.end(); ++it )
        if ( it->IsLeaf() ) {
            terr << "Count Cells to transmit at " << it->center << "\n";
            /* only look at neighboorhood of cells assigned to me */
            int curRank = ((CommData*)it->data[COMM_HEADER_INDEX])->rank;
            if ( curRank != this->rank )
                continue;

            /* Cycle through all neighboring octree nodes of same size */
            assert( T_DIM == 2 );
            int dirs[4] = {RIGHT,LEFT,TOP,BOTTOM};
            for (int i=0; i<compileTime::pow(2,T_DIM); ++i) {
                VecI vDir = getDirectionVector<T_DIM>(dirs[i]);
                typename T_OCTREE::Node & neighbor = *(it->getNeighbor(vDir,periodic));

                terr << "  Looking at neighbor at " << neighbor.center << "\n";

                /* Neighbors may not be leaves (e.g. if it points to a large  *
                 * node with smaller neighbors) , therefore we need to cycle  *
                 * through them again. If it is of same size or even larger   *
                 * and a leaf, then the for loop will end after one iteration *
                 * In Order to check whether the node or its childs are       *
                 * direct neighbors, wie test if getNeighbor in the opposite  *
                 * direction returns the current node (it)                    */
                VecI vOppositeDir = getDirectionVector<T_DIM>(
                                       getOppositeDirection<T_DIM>( dirs[i] ) );
                for ( typename T_OCTREE::iterator it2 = neighbor.begin(ordering);
                      it2 != neighbor.end(); ++it2 )
                if ( it2->IsLeaf() )
                if ( it->IsInside( it2->getNeighbor(vOppositeDir,periodic)->center )
                or   it2->getNeighbor(vOppositeDir,periodic)->IsInside( it->center ) ) {
                    terr << "    Child of neighbor at " << it2->center << "\n";

                    /* If neighbor also belongs to us, then no mpi necessary */
                    int neighborRank = ((CommData*)it2->data[COMM_HEADER_INDEX])->rank;
                    if ( neighborRank == this->rank )
                        continue;

                    /* If two cells have the same neighbor (because they are  *
                     * small, then that neighbor must not tagged in 'ToRecv'  *
                     * twice                                                  */
                    BorderData toRecv = { &(*it2), vOppositeDir };
                    std::list<BorderData> & lsr = neighbors[neighborRank].srecvmats;
                    if ( std::find( lsr.begin(), lsr.end(), toRecv) == lsr.end() )
                        neighbors[neighborRank].srecvmats.push_back(toRecv);

                    /* only send this, if not already send to neighbor thread *
                     * e.g. if it is a bigger node than it's neighbor, and    *
                     * more than one neighbors are on the same process. While *
                     * it is only necessary to recv one cell one time, it may *
                     * be necessary to send one cell to multiple processes    */
                    BorderData toSend = { &(*it), vDir };
                    std::list<BorderData> & lss = neighbors[neighborRank].ssendmats;
                    if ( std::find( lss.begin(), lss.end(), toSend) == lss.end() )
                        neighbors[neighborRank].ssendmats.push_back(toSend);
                }
            }
        }

    /* sort send and recv list in order to not scramble up the order when     *
     * communicating in one large bulk                                        */
    for (int i=0; i<this->worldsize; ++i) {
        neighbors[i].ssendmats.sort();
        neighbors[i].srecvmats.sort();
    }
}

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>:: StartGuardUpdate(int timestep)
{
    #if 1==0
    /* Only cycle through direct neighbors first, see notes on paper !!!  */
    for ( int direction = 1; direction <= this->nneighbors; direction++ )
    {
        /* Limit send directions to very nearest neighbors, no diagonals  */
        /* if ( getAxis(direction) == -1 ) {
            recvrequests[direction] = MPI_REQUEST_NULL;
            sendrequests[direction] = MPI_REQUEST_NULL;
            continue;
        } */

        VecI size = this->getBorderSizeInDirection    (direction);
        VecI pos  = this->getBorderPositionInDirection(direction);

        sendmatrices[direction] = simbox.t[timestep]->cells.getPartialMatrix( pos, size );
        MPI_Isend( sendmatrices[direction].data, sendmatrices[direction].getSize().product()
                   * sizeof(SimBox::CellMatrix::Datatype), MPI_CHAR,
                   neighbors[direction], 53+direction, communicator, &(sendrequests[direction]) );

        #if DEBUG_COMMUNICATOR >= 1
            terr << "[Rank " << this->rank << " in ComBox] Send to Direction ";
            terr << getDirectionVector<T_DIM>( direction );
            terr << "(=lin:" << direction << "=Rank:" << neighbors[direction] << ")";
            terr << " => pos: " << pos << " size: " << size << endl << flush;
            terr << endl << "Sent Matrix: " << endl;
                {const SimBox::CellMatrix & m = sendmatrices[direction];
                VecI size = m.getSize();
                VecI ind(0);
                for ( ind[Y_AXIS]=0; ind[Y_AXIS]<size[Y_AXIS]; ind[Y_AXIS]++) {
                    for ( ind[X_AXIS]=0; ind[X_AXIS]<size[X_AXIS]; ind[X_AXIS]++)
                        terr << m[ ind ].value << " ";
                    terr << endl;
                }
                terr << endl;}
            terr << endl << flush;
        #endif

        int opdir = getOppositeDirection<T_DIM>( direction );
        recvmatrices[opdir] = sendmatrices[direction]; // copy because we want same size
        MPI_Irecv( recvmatrices[opdir].data, recvmatrices[opdir].getSize().product() *
                   sizeof(SimBox::CellMatrix::Datatype), MPI_CHAR,
                   neighbors[opdir], 53+direction, communicator, &(recvrequests[opdir]) );
    }
    #endif
}
