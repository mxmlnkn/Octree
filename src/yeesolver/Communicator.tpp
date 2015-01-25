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
( T_OCTREE & p_tree, const VecI p_periodic )
: rank(0), worldsize(1), processor_name_length(0), tree(p_tree),
  periodic(p_periodic), minLevel(0), maxLevel(0), comDataPtr(NULL),
  distOrdering(Octree::Ordering::Morton), NLeaves(0), NOwnLeaves(0),
  cellsPerOctreeCell(0), guardsize(0), timestepbuffers(0),
  neighbors(NULL), sendrequests(NULL), recvrequests(NULL), timestepSent(-1)
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
        if ( it->data.size() >= 2 )
            delete ((OctCell*)it->data[CELL_DATA_INDEX]);
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
( VecI p_cellsPerOctreeCell, int p_guardsize, int timesteps ) {
    /* number of cells rises exponentially. Noone will use 2^255 Cells! */
    this->minLevel           = tree.getMinLevel();
    this->maxLevel           = tree.getMaxLevel();
    this->cellsPerOctreeCell = p_cellsPerOctreeCell;
    this->guardsize          = p_guardsize;
    this->timestepbuffers    = timesteps;
    /* allocate data which stores assigned ranks and other communication      *
     * information to which pointers will given to octree. And default it to  *
     * the last rank                                                          */
    this->NLeaves            = tree.root->countLeaves();
    this->comDataPtr         = new CommData[NLeaves];
    /* Insert Pointers to array elements of allocated chunk into octree and   *
     * assign weighting and caclulate total weighting and optimalCosts        */
    int dataInserted = 0;
    for ( typename T_OCTREE::iterator it = tree.begin(); it != tree.end(); ++it ) {
        if ( it->IsLeaf() ) {
            assert( dataInserted < NLeaves );
            assert( it->data.empty() );
            it->data.push_back( &(comDataPtr[dataInserted]) );
            ++dataInserted;

            assert( COMM_HEADER_INDEX == 0 );
            /* ToDo: Use weighting operator() !!! */
            ((CommData*)it->data[COMM_HEADER_INDEX])->weighting = 1;
            //((CommData*)it->data[COMM_HEADER_INDEX])->weighting = 1. / it->size.min();
            //((CommData*)it->data[COMM_HEADER_INDEX])->weighting = it->size.min() * pow( 2, tree.getMaxLevel() );
        }
    }
}

/******************************* distributeCells ******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::distributeCells(int ordering) {
    this->distOrdering = ordering;
    
    /* Calculate Costs optimal Costs */
    double totalCosts = 0;
    for ( typename T_OCTREE::iterator it = tree.begin(); it != tree.end(); ++it )
        if ( it->IsLeaf() )
            totalCosts += ((CommData*)it->data[COMM_HEADER_INDEX])->weighting;
    double optimalCosts = totalCosts / double(worldsize);
    tout << "Total Costs: " << totalCosts << " => Optimal Costs: " << optimalCosts << std::endl;
    
    /* Assign cells to all the processes and allocate memory */
    double cumulativeCosts = 0;
    int curRank = 0;
    for ( typename T_OCTREE::iterator it=tree.begin(ordering); it!=tree.end(); ++it )
        if ( it->IsLeaf() ) {
            CommData * curCommData = (CommData*)it->data[COMM_HEADER_INDEX];
            const double curWeighting = curCommData->weighting;
            cumulativeCosts += curWeighting;
            if ( cumulativeCosts >= optimalCosts and curRank != worldsize-1) {
                cumulativeCosts = curWeighting;
                curRank++;
            }
            curCommData->rank = curRank;

            /* Allocate memory for celldata only if cell belongs to us */
            if ( curRank == this->rank ) {
                NOwnLeaves++;
                assert(it->data.size() == CELL_DATA_INDEX);
                VecD cellsize = tree.size / pow( 2, it->getLevel() ) / VecD(cellsPerOctreeCell);
                assert( cellsize == VecD(CELL_SIZE) / pow( 2, it->getLevel() - tree.getMinLevel() ) );

                it->data.push_back(
                    new SimulationBox::SimulationBox<SIMDIM,T_CELLTYPE> (
                        tree.toGlobalCoords(it->center - 0.5*it->size), /* abspos, lower left */
                        cellsPerOctreeCell, /* localcells */
                        cellsize,
                        this->guardsize, this->timestepbuffers )
                );
            }
        }

    /* Count cells to transmit per neighbor process */
    for ( typename T_OCTREE::iterator it=tree.begin(ordering); it != tree.end(); ++it )
        if ( it->IsLeaf() ) {
            /* only look at neighboorhood of cells assigned to me */
            if ( ((CommData*)it->data[COMM_HEADER_INDEX])->rank != this->rank )
                continue;

            /* Cycle through all neighboring octree nodes of same size */
            assert( T_DIM == 2 );
            int dirs[ compileTime::pow(2,T_DIM) ] = {RIGHT,LEFT,TOP,BOTTOM};
            for (int i=0; i < compileTime::pow(2,T_DIM); ++i) {
                VecI vDir = getDirectionVector<T_DIM>(dirs[i]);
                typename T_OCTREE::Node & neighbor = *(it->getNeighbor(vDir,periodic));

                /* Neighbors may be not leaves (e.g. right neighbor of 3 is   *
                 * is parent of 1). In that case iterate over     +-----+-+-+ *
                 * it's children. 4 of 1 will with this strategy  |     |4|1| *
                 * iterate only one time, because  it2.begin()    |  5  +-+-+ *
                 * is a leaf. To check whether the node or its    |     |3|2| *
                 * children are direct neighbors, we test if      +-----+-+-+ *
                 * getNeighbor in the opposite direction returns the current  *
                 * node (it). Because this only works from small to big, also *
                 * test the other way around. E.g. left of 4 is 5, but right  *
                 * of 5 will yield parent of 4                                */
                int oppositeDir = getOppositeDirection<T_DIM>( dirs[i] );
                VecI vOppositeDir = getDirectionVector<T_DIM>(oppositeDir);
                for ( typename T_OCTREE::iterator it2 = neighbor.begin(ordering);
                      it2 != neighbor.end(); ++it2 ) if ( it2->IsLeaf() )
                if ( &neighbor == &(*it2) /* neighbor is leaf */ or
                     it2->getNeighbor(vOppositeDir,periodic) == &(*it) )
                {
                    /* Doesn't matter whether or not it also belongs to us   */
                    int neighborRank = ((CommData*)it2->data[COMM_HEADER_INDEX])->rank;

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
#if 1==0
        /* Allocate OctreeCells / SimBoxes we will receive Borders from */
        for ( typename ToCommList::iterator it = neighbors[i].sRecvData.begin();
              it != neighbors[i].sRecvData.end(); ++it ) {
            /* I don't receive from myself */
            assert( ((CommData*)it->node->data[COMM_HEADER_INDEX])->rank != this->rank );
            /* Because some cells recv more than one border they are multiple *
             * times in recv list                                             */
            if ( it->node->data.size() > CELL_DATA_INDEX)
                continue;
            it->node->data.push_back(
                new SimulationBox::SimulationBox<SIMDIM,T_CELLTYPE> (
                    tree.toGlobalCoords(it->node->center - 0.5*it->node->size),
                    cellsPerOctreeCell, /* localcells */
                    tree.size / pow( 2, it->node->getLevel() ) / cellsPerOctreeCell,
                    this->guardsize, this->timestepbuffers )
            );
        }
#endif
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
        ToCommList & lsr = neighbors[i].sRecvData;
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

        tout << "Send " << cellsToSend << " YeeCells to process " << i
             << " and receive " << cellsToRecv << " from it\n";
    }
}

/****************************** StartGuardUpdate ******************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::StartGuardUpdate(int timestep) {
    this->timestepSent = timestep;
    for ( int neighborRank=0; neighborRank < this->worldsize; ++neighborRank )
    {
        ToCommList & lss = neighbors[neighborRank].sSendData;
        if ( neighborRank == this->rank or lss.empty() ) {
            sendrequests[neighborRank] = MPI_REQUEST_NULL;
            recvrequests[neighborRank] = MPI_REQUEST_NULL;
            continue;
        }

        /* Don't pack and send data to ourselves */
        if (neighborRank == this->rank)
            continue;

        /* copy data needed order as dictated by s*mats into buffers *mats */
        int cellsCopied = 0;
        for ( typename ToCommList::iterator it = lss.begin(); it != lss.end(); ++it ) {
            VecI size = this->getBorderSizeInDirection    (it->direction);
            VecI pos  = this->getBorderPositionInDirection(it->direction);

            CellMatrix tmp = ((OctCell*)it->node->data[CELL_DATA_INDEX])->t[timestep]->cells.getPartialMatrix( pos, size );
            memcpy( &(neighbors[neighborRank].sendData[cellsCopied]), tmp.data,
                    sizeof(T_CELLTYPE) * tmp.size.product() );
            cellsCopied += tmp.size.product();
            assert( cellsCopied <= neighbors[neighborRank].cellsToSend );
        }

        /* Send data Buffer away */
        MPI_Isend( neighbors[neighborRank].sendData, neighbors[neighborRank].
                   cellsToSend * int(sizeof(T_CELLTYPE)), MPI_CHAR, neighborRank, 0,
                   MPI_COMM_WORLD, &(sendrequests[neighborRank]) );
        MPI_Irecv( neighbors[neighborRank].recvData, neighbors[neighborRank].
                   cellsToRecv * int(sizeof(T_CELLTYPE)), MPI_CHAR, neighborRank, 0,
                   MPI_COMM_WORLD, &(recvrequests[neighborRank]) );

#if DEBUG_COMMUNICATOR >= 10
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
    /* set process itself to MPI_Status_Done in order to also interpolate     *
     * in those cases                                                         */

    int rankFinished = -1;
    while(true) {
        if (rankFinished == -1)
            rankFinished = this->rank;
        else {
            /* MPI_Waitany ignores / and sets finished ones to MPI_REQUEST_NULL */
            MPI_Waitany( worldsize, recvrequests, &rankFinished, MPI_STATUSES_IGNORE);
            assert( rankFinished != this->rank );
        }

        /* this index is returned if all recvrequests are MPI_REQUEST_NULL */
        if ( rankFinished == MPI_UNDEFINED )
            break;

        /* only needed for rankFinished == this->rank */
        ToCommList & lss = neighbors[rankFinished].sSendData;
        typename ToCommList::iterator itSend = lss.begin();

        /* copy and interpolate all buffer data into octree cells */
        int cellsCopied = 0;
        ToCommList & lsr = neighbors[rankFinished].sRecvData;
        for ( typename ToCommList::iterator it = lsr.begin(); it != lsr.end(); ++it )
        {
            VecI posB  = this->getBorderPositionInDirection(it->direction);
            VecI sizeB = this->getBorderSizeInDirection    (it->direction);
            CellMatrix tmpB(sizeB);
            /* For borders we transmit to ourselves don't read recvData */
            if ( rankFinished == this->rank ) {
                VecI size = this->getBorderSizeInDirection    (itSend->direction);
                VecI pos  = this->getBorderPositionInDirection(itSend->direction);
                assert(timestepSent >= 0 and timestepSent < timestepbuffers);
                tmpB = ((OctCell*)itSend->node->data[CELL_DATA_INDEX])->
                        t[timestepSent]->cells.getPartialMatrix( pos, size );
            } else
                memcpy( tmpB.data, &(neighbors[rankFinished].recvData[cellsCopied]),
                        sizeof(T_CELLTYPE) * tmpB.size.product() );
            cellsCopied += tmpB.size.product();
            assert( cellsCopied <= neighbors[rankFinished].cellsToRecv );

            /* Copy the newly received borders into the correct guards. For   *
             * that we first need to calculate alle the guards we need to put *
             * the received data into. Do this by traversing over the         *
             * neighbor:                                       +-----+-----+  *
             *  - 1-7 is the traversal order and also the id   |     |4b|1a|  *
             *  - a-c is the assigned process id               |  5b |--+--|  *
             *  - process a gets borders from 4b,3b and 7b     |     |3b|2a|  *
             *  - border 2a to 1a and vice-versa is being      |-----+-----+  *
             *    interpolated (in this case just copied) in   |     |     |  *
             *    StartGuardUpdate, because those don't need   |  6c |  7c |  *
             *    to be communicated.                          |     |     |  *
             * 4b will check neighbors: 5b,3b,1a,7c only       +-----+-----+  *
             * 1a belongs to us (a) and needs to be updated with data 4b.     *
             * For 7c top border the neighbors will be: 5b,6c and the         *
             * parent of  4b,1a,... (actually twice: neighbor above and       *
             * below if  periodic ) In that case we need to travers all       *
             * children and check their neighbor in opposite direction is     *
             * 7c. For 7c top border this will yield 2a and 3b, latter is     *
             * ignored, because I am process a.                               */

            VecI vDir = getDirectionVector<T_DIM>(it->direction);
            typename T_OCTREE::Node * neighbor = it->node->getNeighbor( vDir, periodic );
            int oppositeDir = getOppositeDirection<T_DIM>(it->direction);
            VecI vOppositeDir = getDirectionVector<T_DIM>(oppositeDir);

            VecI posG  = this->getGuardPositionInDirection(oppositeDir);
            VecI sizeG = this->getGuardSizeInDirection    (oppositeDir);
            assert( sizeG == sizeB );
            for ( typename T_OCTREE::iterator itT = neighbor->begin(); itT != neighbor->end(); ++itT )
            if ( itT->IsLeaf() and ((CommData*)itT->data[COMM_HEADER_INDEX])->rank == this->rank )
                if ( neighbor == &(*itT) /* neighbor is leaf */
                or   itT->getNeighbor(vOppositeDir,periodic) == it->node )
            {
                /* E.g. send 5b left to 1a right, then size of 1a stays, but  *
                 * we only need the upper half of left border of 5b           */
                VecI posS, sizeS, posT, sizeT;
                VecD lowerleftSource = it->node->center - it->node->size / 2;
                VecD lowerleftTarget = itT->center - itT->size / 2;
                /* automatic conversion makes (-1,0) to (1,0), because    *
                 * bool(-1)=1. Apply calculated offset in all directions  *
                 * except the one, we want to send to / the neighbor is   */
                VecI mask = VecI(1) - VecI( (Vec<bool,T_DIM>)(vDir) );
                assert( mask.sum() == T_DIM-1 );
                if ( /* received data */ it->node->size > itT->size /* owned neighbor */ ) {
                    sizeS  = (VecI(1) - mask) * sizeB;
                    for (int i = 0; i < itT->size.dim; ++i) {
                        assert( fmod( it->node->size[i], itT->size[i] ) == 0.0 );
                        assert( mask[i] * sizeB[i] % int( it->node->size[i] / itT->size[i] ) == 0 );
                    }
                    sizeS += mask * sizeB / VecI( it->node->size  / itT->size );
                    sizeT  = sizeG;
                    VecI offset = ( lowerleftTarget - lowerleftSource ) /
                                  it->node->size * VecD(cellsPerOctreeCell);
                    posS  = posB + offset * mask;
                    posT  = posG;
                } else { /* this case includes same-size-case */
                    sizeS  = sizeB;
                    sizeT  = (VecI(1) - mask) * sizeG;
                    sizeT += mask * sizeG / VecI( itT->size / it->node->size );
                    VecI offset = ( lowerleftSource - lowerleftTarget ) /
                                  itT->size * VecD(cellsPerOctreeCell);
                    posS  = posB;
                    posT  = posG + offset * mask;;
                }
                CellMatrix tmpS(sizeS);
                CellMatrix tmpT(sizeT);
                assert(sizeS != VecI(0));
                assert(sizeT != VecI(0));

#if DEBUG_COMMUNICATOR >= 10
                /* Debug Output for interpolating borders to guards */
                tout << "Interpolate matrix at " << posS << " sized " << sizeS
                     << " cut out from ";
                switch ( it->direction ) {
                    case LEFT  : tout << "Left";   break;
                    case RIGHT : tout << "Right";  break;
                    case BOTTOM: tout << "Bottom"; break;
                    case TOP   : tout << "Top";    break;
                }
                tout << " Border of node at " << it->node->center << " sized "
                     << it->node->size << "\n -> to ";
                tout << "matrix at " << posT << " sized " << sizeT
                     << " cut out from ";
                switch ( oppositeDir ) {
                    case LEFT  : tout << "Left";   break;
                    case RIGHT : tout << "Right";  break;
                    case BOTTOM: tout << "Bottom"; break;
                    case TOP   : tout << "Top";    break;
                }
                tout << " Guard of node at " << itT->center << " sized "
                     << itT->size << "\n";
#endif

                /* Because we only received the border, not the whole cell    *
                 * matrix, posB must not be used to access data, just to      *
                 * calculate the offset above. Instead 'posB' is 'offset'     */
                tmpS = tmpB.getPartialMatrix( posS - posB, sizeS );
                if ( sizeS != sizeT )
                    tmpS.NearestResizeTo( tmpT );
                else
                    tmpT = tmpS;

                ((OctCell*)itT->data[CELL_DATA_INDEX])->t[timestep]->cells.insertMatrix(posT,tmpT);
            }
            if (rankFinished == this->rank) {
                assert(itSend != lss.end());
                ++itSend;
            }
        }
    } // MPI_Waitany
}

/********************************* PrintPNG ***********************************/
template<int T_DIM, typename T_OCTREE, typename T_CELLTYPE>
template<typename T_FUNC>
void OctreeCommunicator<T_DIM,T_OCTREE,T_CELLTYPE>::PrintPNG
(int timestep, const char * name, T_FUNC colorFunctor )
{
    assert( T_DIM == 2 );
    VecI sizepx = this->cellsPerOctreeCell * int(pow(2,this->maxLevel));
    time_t t = time(0);
    struct tm * now = localtime( &t );
    static std::stringstream s_timestamp;
    if ( s_timestamp.rdbuf()->in_avail() == 0 ) {
        s_timestamp << 1900 + now->tm_year << "-" << 1 + now->tm_mon 
                    << "-" << now->tm_mday << "_" << now->tm_hour << "-"
                    << now->tm_min << "_output";
        tout << "Creating directory " << boost::filesystem::absolute(s_timestamp.str()) << " with boost\n";
        boost::filesystem::create_directory( boost::filesystem::absolute(s_timestamp.str()) );
    }
    std::stringstream filenamepng;
    filenamepng << s_timestamp.str() << "/" << name << "_rank-" << this->rank << "_t" << timestep << ".png";
    pngwriter image( sizepx[0],sizepx[1], 1.0, filenamepng.str().c_str() );
#if DEBUG_COMMUNICATOR >= 10
    tout << "Create " << sizepx << "px sized png named " << filenamepng.str() << "\n";
#endif
    for ( typename T_OCTREE::iterator it=tree.begin(); it!=tree.end(); ++it )
    if ( it->IsLeaf() ) if ( it->data.size() >= 2 )
    {
        OctCell & data   = *((OctCell*)it->data[this->CELL_DATA_INDEX]);

        /* No Resize, if uniform octree cells. If not uniform, then scale up  *
         * smaller ones, by simpling filling the rest of the space            */
        int resizeFactor = (int) pow( 2, this->maxLevel - it->getLevel() );

        /* shift internal octree coords from [-0.5,0.5] to [0,1.0] then       *
         *  then shift it->center to it->lower left corner : - it.size/2      *
         *  then get octree cell index in that level : / it.size              *
         *  then scale up to internal cells which will be pixels: *localcells *
         *  then scale up pixels in that level to maxLevel : * resizeFactor   */
        VecI abspos = ( it->center + tree.root->size/2 - it->size/2 ) / it->size;
        for ( int i=0; i<T_DIM; ++i )
            assert( fmod( (it->center + tree.root->size/2 - it->size/2)[i], it->size[i] ) == 0 );
        abspos *= this->cellsPerOctreeCell * resizeFactor;

        /* abspos member of SimulationBox is initialized bei Communicator.tpp *
         * with it->center - 0.5*it->size, meaning lower left corner with     *
         * internal units of OctreeNode ( no rounding errors should happen )  */
        for ( typename OctCell::IteratorType itm = data.getIterator( timestep, 
              SimulationBox::CORE + SimulationBox::BORDER ).begin();
              itm != itm.end(); ++itm )
        {
            /* pngwriter begins counting pixels with 1 instead of 0 -.- */
            VecI pos = 1 + abspos + resizeFactor*( itm.icell - itm.guardsize );
            VecI posTo = pos + resizeFactor - 1;
            Vec<double,3> color = colorFunctor(*itm);
            image.filledsquare( pos[0],pos[1], posTo[0],posTo[1],
                                color[0], color[1], color[2] );
        }
    }
    image.close();
}
