#pragma once

/******************************* Constructor ******************************/
template<typename T_OCTREE>
OctreeCommunicator<T_OCTREE>::OctreeCommunicator(void) {
    std::cerr << "Communicator needs to be given an octree handle when constructed!\n";
}

template<typename T_OCTREE>
OctreeCommunicator<T_OCTREE>::OctreeCommunicator( T_OCTREE & tree )
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
template<typename T_OCTREE>
OctreeCommunicator<T_OCTREE>::~OctreeCommunicator() {
    if ( comDataPtr != NULL )
        free(comDataPtr);
#if 1==0
    free( this->neighbors );
    free( this->sendrequests );
    free( this->recvrequests );
    free( this->sendmatrices );
    free( this->recvmatrices );
#endif
}

template<typename T_OCTREE>
void OctreeCommunicator<T_OCTREE>::initCommData(void) {
    /* allocate data which stores assigned ranks and other communication      *
     * information to which pointers will given to octree. And default it to  *
     * the last rank                                                          */
    this->NLeaves = tree.root->countLeaves();
    this->comDataPtr = (CommData*) malloc( sizeof(CommData)*NLeaves );
    /* Insert Pointers to array elements of allocated chunk into octree,      *
     * assign weighting and rank (default last rank) and caclulate total weig.*/
    int dataInserted = 0;
    this->totalCosts = 0;
    for ( typename T_OCTREE::iterator it = tree.begin(); it != tree.end(); ++it )
        if ( it->IsLeaf() ) {
            assert( dataInserted < NLeaves );
            assert( it->data.empty() );
            it->data.push_back( &(comDataPtr[dataInserted]) );
            ++dataInserted;

            assert( COMM_HEADER_INDEX == 0 );
            ((CommData*)it->data[COMM_HEADER_INDEX])->rank      = worldsize-1;
            ((CommData*)it->data[COMM_HEADER_INDEX])->weighting = 1. / it->size.min();

            this->totalCosts += ((CommData*)it->data[COMM_HEADER_INDEX])->weighting;
        }
    this->optimalCosts = totalCosts / double(worldsize);
    tout << "Total Costs: " << totalCosts << " => Optimal Costs: " << optimalCosts << std::endl;
}
