#pragma once

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
const int OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::COMM_HEADER_INDEX = 0;

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
const int OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::CELL_DATA_INDEX = 1;

    /* squash the cube (CORE+BORDER) we want to send to guardsize in each     *
     * direction we want to sent it to, e.g.                                  *
     *   guard = 1, core+border = 3 => begin with 3x3x3 cube.                 *
     *   send to LEFT or RIGHT => make it a 1x3x3 surface                     *
     *   send to TOP or BOTTOM => make it a 3x1x3 surface                     *
     *   send to LEFT+TOP      => make it a 1x1x3 line of cells (edge)        */
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
typename OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::VecI
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::getBorderSizeInDirection
( const int & direction ) const
{
    VecI size = cellsPerOctreeCell;
    VecI vDir( getDirectionVector<T_DIM>( direction ) );
    for (int i=0; i<T_DIM; ++i)
        if ( vDir[i] != 0 )
            size[i] = guardsize;
    return size;
}

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
typename OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::VecI
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::getBorderPositionInDirection
( const int & direction ) const
{
    VecI pos ( guardsize );
    VecI vDir( getDirectionVector<T_DIM>( direction ) );
    for ( int i=0; i<T_DIM; ++i )
        if ( vDir[i] == 1 )
            pos[i] += cellsPerOctreeCell[i] - guardsize;
    //assert( simbox.inArea( pos, SimulationBox::BORDER ) );
    return pos;
}

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
typename OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::VecI
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::getGuardSizeInDirection
( const int & direction ) const
{
    return this->getBorderSizeInDirection( direction );
}

template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
typename OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::VecI
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::getGuardPositionInDirection
( const int & direction ) const
{
    VecI pos = this->getBorderPositionInDirection( direction );
    VecI vDir( getDirectionVector<T_DIM>( direction ) );
    for ( int i=0; i<T_DIM; ++i ) {
        assert( abs(vDir[i]) <= 1 );
        pos[i] += vDir[i]*guardsize;
    }
    return pos;
}

/********************************* Constructor ********************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::OctreeCommunicator
( T_OCTREE & tree, const VecI periodic )
: rank(0), worldsize(1), tree(tree), periodic(periodic),
  minLevel(0), maxLevel(0), comDataPtr(NULL),
  NLeaves(0), cellsPerOctreeCell(0), guardsize(0), timestepbuffers(0),
  neighbors(NULL), sendrequests(NULL), recvrequests(NULL)
{
    /* Initialize MPI */
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &processor_name_length);

    this->neighbors    = new NeighborData[worldsize];
    this->sendrequests = new MPI_Request[worldsize];
    this->recvrequests = new MPI_Request[worldsize];
    if ( neighbors == NULL or sendrequests == NULL or recvrequests == NULL )
        terr << "Couldn't allocate Memory!";
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

    for (int i=0; i < this->worldsize; ++i) {
        if( neighbors[i].sendData != NULL )
            free(neighbors[i].sendData);
        if( neighbors[i].recvData != NULL )
            free(neighbors[i].recvData);
    }
    delete[] this->recvrequests;
    delete[] this->sendrequests;
    delete[] this->neighbors;
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

                /* Neighbors may not be leaves (e.g. if it points to a large  *
                 * node with smaller neighbors) , therefore we need to cycle  *
                 * through them again. If it is of same size or even larger   *
                 * and a leaf, then the for loop will end after one iteration *
                 * In Order to check whether the node or its children are     *
                 * direct neighbors, wie test if getNeighbor in the opposite  *
                 * direction returns the current node (it)                    */
                int oppositeDir = getOppositeDirection<T_DIM>( dirs[i] );
                VecI vOppositeDir = getDirectionVector<T_DIM>(oppositeDir);
                for ( typename T_OCTREE::iterator it2 = neighbor.begin(ordering);
                      it2 != neighbor.end(); ++it2 )
                if ( it2->IsLeaf() )
                if ( it->IsInside( it2->getNeighbor(vOppositeDir,periodic)->center )
                or   it2->getNeighbor(vOppositeDir,periodic)->IsInside( it->center ) )
                {
                    /* If neighbor also belongs to us, then no mpi necessary */
                    int neighborRank = ((CommData*)it2->data[COMM_HEADER_INDEX])->rank;
                    if ( neighborRank == this->rank )
                        continue;

                    /* If two cells have the same neighbor (because they are  *
                     * small, then that neighbor must not tagged in 'ToRecv'  *
                     * twice                                                  */
                    BorderData toRecv = { &(*it2), oppositeDir };
                    ToCommList & lsr = neighbors[neighborRank].sRecvData;
                    if ( std::find( lsr.begin(), lsr.end(), toRecv) == lsr.end() )
                        neighbors[neighborRank].sRecvData.push_back(toRecv);

                    /* only send this, if not already send to neighbor thread *
                     * e.g. if it is a bigger node than it's neighbor, and    *
                     * more than one neighbors are on the same process. While *
                     * it is only necessary to recv one cell one time, it may *
                     * be necessary to send one cell to multiple processes    */
                    BorderData toSend = { &(*it), dirs[i] };
                    ToCommList & lss = neighbors[neighborRank].sSendData;
                    if ( std::find( lss.begin(), lss.end(), toSend) == lss.end() )
                        neighbors[neighborRank].sSendData.push_back(toSend);
                }
            }
        }

    /* sort send and recv list in order to not scramble up the order when     *
     * communicating in one large bulk. Also allocate packed buffer          */
    for (int i=0; i<this->worldsize; ++i) {
        neighbors[i].sSendData.sort();
        neighbors[i].sRecvData.sort();

        /* Count how many doubles to allocate. If cellsPerOctreeCell is the   *
         * same in every dimension, it would be easy, but if not, then we     *
         * need to count how often every direction is accessed / communicated */
        int cellsToSend = 0;
        ToCommList & lss = neighbors[i].sSendData;
        for ( typename ToCommList::iterator it = lss.begin(); it != lss.end(); ++it ) {
            cellsToSend += guardsize * getBorderSizeInDirection(it->direction).product();
            /* old formula with which I want to compare only works for non-   *
             * diagonal directions                                            */
            if ( getDirectionVector<T_DIM>(it->direction).sum() == 1 )
                assert( getBorderSizeInDirection(it->direction).product() ==
                        cellsPerOctreeCell.product() /
                        cellsPerOctreeCell[ getAxis<T_DIM>(it->direction) ] );
        }
        neighbors[i].cellsToSend = cellsToSend;

        int cellsToRecv = 0;
        ToCommList & lsr = neighbors[i].sSendData;
        for ( typename ToCommList::iterator it = lsr.begin(); it != lsr.end(); ++it ) {
            cellsToRecv += guardsize * getBorderSizeInDirection(it->direction).product();
            if ( getDirectionVector<T_DIM>(it->direction).sum() == 1 )
                assert( getBorderSizeInDirection(it->direction).product() ==
                        cellsPerOctreeCell.product() /
                        cellsPerOctreeCell[ getAxis<T_DIM>(it->direction) ] );
        }
        neighbors[i].cellsToRecv = cellsToRecv;

        /* actually allocate counted bytes */
        if( neighbors[i].sendData != NULL )
            free(neighbors[i].sendData);
        if( neighbors[i].recvData != NULL )
            free(neighbors[i].recvData);
        neighbors[i].recvData = (T_CELLTYPE*) malloc( sizeof(T_CELLTYPE) * cellsToRecv );
        neighbors[i].sendData = (T_CELLTYPE*) malloc( sizeof(T_CELLTYPE) * cellsToSend );

        tout << "Send " << cellsToSend << " doubles to process " << i
             << "and receive " << cellsToRecv << " from it\n";
    }
}

/****************************** StartGuardUpdate ******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::StartGuardUpdate(int timestep)
{
    /* copy data needed order as dictated by sxxxxmats into buffers xxxxnats */
    for ( int neighborRank=0; neighborRank < this->worldsize; ++neighborRank ) {
        ToCommList & lss = neighbors[neighborRank].sSendData;
        if ( neighborRank == this->rank or lss.empty() ) {
            sendrequests[neighborRank] = MPI_REQUEST_NULL;
            recvrequests[neighborRank] = MPI_REQUEST_NULL;
            continue;
        }
        int cellsCopied = 0;
        for ( typename ToCommList::iterator it = lss.begin(); it != lss.end(); ++it ) {
            VecI size = this->getBorderSizeInDirection    (it->direction);
            VecI pos  = this->getBorderPositionInDirection(it->direction);

            CellMatrix tmp = ((OctCell*)it->node->data[CELL_DATA_INDEX])->t[timestep]->cells.getPartialMatrix( pos, size );
            memcpy( &(neighbors[neighborRank].sendData[cellsCopied]), tmp.data,
                    sizeof(T_CELLTYPE) * tmp.size.product()         );
            cellsCopied += tmp.size.product();
            assert( cellsCopied <= neighbors[neighborRank].cellsToSend );
        }
        MPI_Isend( neighbors[neighborRank].sendData, neighbors[neighborRank].cellsToSend * sizeof(T_CELLTYPE),
                   MPI_CHAR, neighborRank, 0, MPI_COMM_WORLD, &(sendrequests[neighborRank]) );
        MPI_Irecv( neighbors[neighborRank].recvData, neighbors[neighborRank].cellsToRecv * sizeof(T_CELLTYPE),
                   MPI_CHAR, neighborRank, 0, MPI_COMM_WORLD, &(recvrequests[neighborRank]) );

#if DEBUG_COMMUNICATOR >= 1
        /* Debug Output of recv- and send-list in comBox */
        tout << "I send to process " << neighborRank << ":\n";
        for ( typename ToCommList::iterator it = neighbors[neighborRank].sSendData.begin();
              it != neighbors[neighborRank].sSendData.end(); ++it ) {
            tout << "    ";
            switch ( it->direction ) {
                case LEFT  : tout << "Left";   break;
                case RIGHT : tout << "Right";  break;
                case BOTTOM: tout << "Bottom"; break;
                case TOP   : tout << "Top";    break;
            }
            tout << " side of node at " << it->node->center << "\n";
        }
        tout << "I receive from process " << neighborRank << ":\n";
        for ( typename ToCommList::iterator it = neighbors[neighborRank].sRecvData.begin();
              it != neighbors[neighborRank].sRecvData.end(); ++it ) {
            tout << "    ";
            switch ( it->direction ) {
                case LEFT  : tout << "Left";   break;
                case RIGHT : tout << "Right";  break;
                case BOTTOM: tout << "Bottom"; break;
                case TOP   : tout << "Top";    break;
            }
            tout << " side of node at " << it->node->center << "\n";
        }
#endif
    }
}

/***************************** FinishGuardUpdate ******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::FinishGuardUpdate( int timestep ) {
    while(true) {
        int rankFinished;
        /* MPI_Waitany ignores / and sets finished ones to MPI_REQUEST_NULL */
        MPI_Waitany( worldsize, recvrequests, &rankFinished, MPI_STATUSES_IGNORE);
        assert( rankFinished != this->rank );
        /* this index is returned if all recvrequests are MPI_REQUEST_NULL */
        if ( rankFinished == MPI_UNDEFINED )
            break;

        /* copy all buffer data into octree cells */
        ToCommList & lsr = neighbors[rankFinished].sRecvData;
        int cellsCopied = 0;
        for ( typename ToCommList::iterator it = lsr.begin(); it != lsr.end(); ++it )
        {
            VecI size = this->getBorderSizeInDirection    (it->direction);
            VecI pos  = this->getBorderPositionInDirection(it->direction);

            CellMatrix tmp(size);
            memcpy( tmp.data, &(neighbors[rankFinished].recvData[cellsCopied]),
                    sizeof(T_CELLTYPE) * tmp.size.product() );
            ((OctCell*)it->node->data[CELL_DATA_INDEX])->t[timestep]->cells.insertMatrix(pos,tmp);

            cellsCopied += tmp.size.product();
            assert( cellsCopied <= neighbors[rankFinished].cellsToRecv );
        }

        /* interpolate or copy Borders of neighboring cells, which are now    *
         * existent and up-to-date, even if the octree cells are on other     *
         * nodes                                                              */
        tout << "ToDo: Interpolation and Tests!\n";

        /*VecI size = getGuardSizeInDirection    (direction);
        VecI pos  = getGuardPositionInDirection(direction);
        simbox.t[timestep]->cells.insertMatrix( pos, recvmatrices[direction] );*/
    }
}

/********************************* PrintPNG ***********************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::PrintPNG(int timestep, const char * name )
{
    const int X = 0;
    const int Y = 1;

    assert( T_DIM == 2 );
    VecI sizepx = this->cellsPerOctreeCell * pow(2,this->maxLevel);
    time_t t = time(0);
    struct tm * now = localtime( &t );
    std::stringstream filenamepng;
    filenamepng << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
                << now->tm_mday << "_" << now->tm_hour << "-"
                << now->tm_min << "_";
    filenamepng << name << "_rank-" << this->rank << "_Ex.png";
    pngwriter image( sizepx[X],sizepx[Y], 1.0, filenamepng.str().c_str() );
    tout << "Create " << sizepx << "px sized png named " << filenamepng.str() << "\n";

    for ( typename T_OCTREE::iterator it=tree.begin(); it!=tree.end(); ++it )
    if ( it->IsLeaf() ) if ( it->data.size() >= 2 )
    {
        OctCell & data   = *((OctCell*)it->data[this->CELL_DATA_INDEX]);
        BaseMatrix<YeeCell,SIMDIM> & matrix = data.t[0]->cells;

        /* No Resize, if uniform octree cells. If not uniform, then scale up  *
         * smaller ones, by simpling filling the rest of the space            */
        int resizeFactor = pow( 2, this->maxLevel - it->getLevel() );

        /* shift internal octree coords from [-0.5,0.5] to [0,1.0] then       *
         *  then shift it->center to it->lower left corner : - it.size/2      *
         *  then get octree cell index in that level : / it.size              *
         *  then scale up to internal cells which will be pixels: *localcells *
         *  then scale up pixels in that level to maxLevel : * resizeFactor   */
        VecI abspos = ( it->center + tree.root->size/2 - it->size/2 ) / it->size;
        assert( fmod( (it->center + tree.size/2 - it->size/2)[X], it->size[X] ) == 0 );
        assert( fmod( (it->center + tree.size/2 - it->size/2)[Y], it->size[Y] ) == 0 );
        abspos *= this->cellsPerOctreeCell * resizeFactor;

        /* abspos member of SimulationBox is initialized bei Communicator.tpp *
         * with it->center - 0.5*it->size, meaning lower left corner with     *
         * internal units of OctreeNode ( no rounding errors should happen )  */
        typename OctCell::IteratorType it = data.getIterator( 0, SimulationBox::CORE + SimulationBox::BORDER );
        for ( it = it.begin(); it != it.end(); ++it ) {
            /* pngwriter begins counting pixels with 1 instead of 0 -.- */
            VecI pos = 1 + abspos + resizeFactor*(it.icell-it.guardsize);
            VecI posTo = pos + resizeFactor - 1;
            image.filledsquare( pos[X],pos[Y], posTo[X],posTo[Y],
                matrix[it.icell].E[0][X], 0.0, -matrix[it.icell].E[0][X]);
        }
    }
    image.close();
}
