#pragma once

struct CellIterator {
    const OctreeCommunicator & combox;
    typename T_OCTREE::iterator it;
    typename OctCell::IteratorType itm;
    int timestep, area, ordering;

    /******************************** Constructor *********************************/
    CellIterator ( const OctreeCommunicator & pcombox, int parea, int ptimestep, int pordering )
    : combox(pcombox), it(combox.tree.begin(pordering)), itm(),
      timestep(ptimestep), area(parea), ordering(pordering)
    {}

    CellIterator begin( void ) const {
        CellIterator res(*this);

        for ( res.it = combox.tree.begin( ordering ); res.it != res.it.end(); ++res.it ) {
            if ( res.it->IsLeaf() ) {
                CommData * commdata = (CommData*) res.it->data[COMM_HEADER_INDEX];
                if ( commdata->rank == combox.rank ) {
                    OctCell * data = (OctCell*) res.it->data[CELL_DATA_INDEX];
                    res.itm = data->getIterator( timestep, area ).begin();
                    break;
                }
            }
        }

        if ( res.it == res.it.end() ) {
            res.itm = res.itm.end();
        }

        return res;
    }
    /* return VecI(-1) which is a end marker. Will have to set it to that *
     * in operator++ for example                                          */
    CellIterator end( void ) const {
        CellIterator tmp( *this );
        tmp.it  = tmp.it.end();
        tmp.itm = tmp.itm.end();
        return tmp;
    }

    CellIterator & operator++( void ) {
        if ( it == it.end() )
            return *this;
        if ( itm != itm.end() ) {
            ++itm;
            if ( itm != itm.end() )
                return *this;
        }

        //tout << "Increment it\n";
        /* only iterate further over octree if matrix iterator reached end */
        for ( it = ++it; it != it.end(); ++it ) {
            if ( it->IsLeaf() ) {
                CommData * commdata = (CommData*) it->data[COMM_HEADER_INDEX];
                if ( commdata->rank == combox.rank ) {
                    OctCell * data = (OctCell*) it->data[CELL_DATA_INDEX];
                    itm = data->getIterator( timestep, area ).begin();
                    break;
                }
            }
        }

        /* can also be reached, if itm=itm.end() and it=it.end() */
        return *this;
    }

    bool operator==( const CellIterator & rhs ) {
        bool alland = true;
        alland = alland & ( this->timestep  == rhs.timestep  );
        alland = alland & ( this->area      == rhs.area      );
        alland = alland & ( this->ordering  == rhs.ordering  );
        alland = alland & ( this->it        == rhs.it        );
        alland = alland & ( this->itm       == rhs.itm       );
        return alland;
    }

    bool operator!=( const CellIterator & rhs ) {
        return !( (*this) == rhs );
    }

    T_CELLTYPE& operator*( void ) const {
        return *itm;
    }

    T_CELLTYPE* operator->( void ) const {
        return &(*itm);
    }

    VecD getGlobalPosition(void) const {
        OctCell * simbox = (OctCell*) it->data[CELL_DATA_INDEX];
        return simbox->getGlobalPosition( itm );
    }
};

