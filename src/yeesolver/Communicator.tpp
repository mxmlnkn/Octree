#pragma once

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
const int OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::COMM_HEADER_INDEX = 0;

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
const int OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::CELL_DATA_INDEX = 1;

/******************************* Constructor ******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::OctreeCommunicator(void) {
    std::cerr << "Communicator needs to be given an octree handle when constructed!\n";
}

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::OctreeCommunicator( T_OCTREE & tree )
: rank(0), worldsize(1), tree(tree), comDataPtr(NULL), NLeaves(0)
{
    /* Initialize MPI */
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &processor_name_length);

#if 1==0
    /* TODO */
    nneighbors = 0;
    terr << "Number of neighboring processes     : " << nneighbors;
    terr << "Border size to send ( per process ) : ";

    this->nneighbors = nneighbors;
    this->neighbors    = (int*)         malloc( sizeof(int)        *(nneighbors+1) );
    this->recvrequests = (MPI_Request*) malloc( sizeof(MPI_Request)*(nneighbors+1) );
    this->sendrequests = (MPI_Request*) malloc( sizeof(MPI_Request)*(nneighbors+1) );
    this->sendmatrices = (SimBox::CellMatrix*) malloc( sizeof(SimBox::CellMatrix)*(nneighbors+1) );
    this->recvmatrices = (SimBox::CellMatrix*) malloc( sizeof(SimBox::CellMatrix)*(nneighbors+1) );
    if ( neighbors == NULL or recvrequests == NULL or sendrequests == NULL )
        terr << "Couldn't allocate Memory!";
    for (int i=0; i<=nneighbors; i++) {
        sendrequests[i] = MPI_REQUEST_NULL;
        recvrequests[i] = MPI_REQUEST_NULL;
    }
#endif
}

/******************************** Destructor ******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::~OctreeCommunicator() {
    if ( comDataPtr != NULL )
        free(comDataPtr);

    for ( typename T_OCTREE::iterator it=tree.begin(); it!=tree.end(); ++it ) {
        if ( it->data.size() >= 2 ) {
            assert( ((CommData*)it->data[COMM_HEADER_INDEX])->rank == this->rank );
            free( (T_CELLTYPE*) it->data[CELL_DATA_INDEX] );
        }
    }

#if 1==0
    free( this->neighbors );
    free( this->sendrequests );
    free( this->recvrequests );
    free( this->sendmatrices );
    free( this->recvmatrices );
#endif
}

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::initCommData(void) {
    /* allocate data which stores assigned ranks and other communication      *
     * information to which pointers will given to octree. And default it to  *
     * the last rank                                                          */
    this->NLeaves = tree.root->countLeaves();
    this->comDataPtr = (CommData*) malloc( sizeof(CommData)*NLeaves );
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

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::distributeCells(VecI cellsPerOctreeCell, int ordering) {
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
            /* Allocate memory for celldata */
            if ( curRank == this->rank ) {
                NOwnLeaves++;
                it->data.push_back( new BaseMatrix<T_CELLTYPE,T_DIM>(cellsPerOctreeCell) );
            }
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
