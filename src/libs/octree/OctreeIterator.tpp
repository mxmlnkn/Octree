#pragma once

#define DEBUG_OCTREE_ITERATOR 0


namespace Octree {

template<int T_DIM>
typename OctreeIterator<T_DIM>::orderingTableStruct
    OctreeIterator<T_DIM>::orderingTable;

template<>
typename OctreeIterator<2>::orderingTableStruct
    OctreeIterator<2>::orderingTable =
          {
            /* Gray-Code ordering */
            { /* ordering Table */
              { /* Row 0 meaning parent has orientation 0 */
                /* Columns. 0th column contains index and orientation */
                { 0,1,3,2 }, /* ordering if orientation is 0 */
                { 3,2,0,1 }  /* ordering if orientation is 1 */
              },
              /* orientation Table */
              { /* Row 0 meaning parent has orientation 0 */
                /* Columns. 0th column contains index and orientation */
                { 0,1,1,0 }, /* orientation for children if orientation is 0 */
                { 1,0,0,1 }  /* orientation for children if orientation is 1 */
              }
            },
            /* Hilbert ordering */
            { /* ordering Table */
              {
                { 0,1,3,2 },
                { 0,2,3,1 },
                { 3,1,0,2 },
                { 3,2,0,1 }
              },
              /* orientation Table */
              {
                { 1,0,0,2 },
                { 0,1,1,3 },
                { 3,2,2,0 },
                { 2,3,3,1 }
              }
            }
          };

template<>
typename OctreeIterator<3>::orderingTableStruct
    OctreeIterator<3>::orderingTable =
          {
            /* Gray-Code ordering */
            { /* ordering Table */
              { /* Row 0 meaning parent has orientation 0 */
                /* Columns. 0th column contains index and orientation */
                { 0,1,3,2, 6,7,5,4 }, /* ordering if orientation is 0 */
                { 5,4,6,7, 3,2,0,1 }, /* ordering if orientation is 1 */
                { 3,2,0,1, 5,4,6,7 }, /* ordering if orientation is 2 */
                { 6,7,5,4, 0,1,3,2 }  /* ordering if orientation is 3 */
              },
              /* orientation Table */
              { /* Row 0 meaning parent has orientation 0 */
                /* Columns. 0th column contains index and orientation */
                { 0,1,2,3, 3,2,1,0 }, /* orientation for children if orientation is 0 */
                { 1,0,3,2, 2,3,0,1 }, /* orientation for children if orientation is 1 */
                { 2,3,0,1, 1,0,3,2 }, /* orientation for children if orientation is 2 */
                { 3,2,1,0, 0,1,2,3 }  /* orientation for children if orientation is 3 */
              }
            },
            /* Hilbert ordering */
            { /* ordering Table */
              {
                {0,1,3,2, 6,7,5,4 }, /* orientation  0 */
                {0,4,6,2, 3,7,5,1 }, /* orientation  1 */
                {0,1,5,4, 6,7,3,2 }, /* orientation  2 */
                {5,1,0,4, 6,2,3,7 }, /* orientation  3 */
                {3,7,6,2, 0,4,5,1 }, /* orientation  4 */
                {6,7,3,2, 0,1,5,4 }, /* orientation  5 */
                {5,1,3,7, 6,2,0,4 }, /* orientation  6 */
                {0,4,5,1, 3,7,6,2 }, /* orientation  7 */
                {5,4,0,1, 3,2,6,7 }, /* orientation  8 */
                {5,4,6,7, 3,2,0,1 }, /* orientation  9 */
                {0,2,3,1, 5,7,6,4 }, /* orientation 10 */
                {6,4,0,2, 3,1,5,7 }, /* orientation 11 */
                {5,7,3,1, 0,2,6,4 }, /* orientation 12 */
                {3,7,5,1, 0,4,6,2 }, /* orientation 13 */
                {6,4,5,7, 3,1,0,2 }, /* orientation 14 */
                {0,2,6,4, 5,7,3,1 }, /* orientation 15 */
                {6,2,0,4, 5,1,3,7 }, /* orientation 16 */
                {6,2,3,7, 5,1,0,4 }, /* orientation 17 */
                {3,2,0,1, 5,4,6,7 }, /* orientation 18 */
                {6,7,5,4, 0,1,3,2 }, /* orientation 19 */
                {5,7,6,4, 0,2,3,1 }, /* orientation 20 */
                {3,2,6,7, 5,4,0,1 }, /* orientation 21 */
                {3,1,0,2, 6,4,5,7 }, /* orientation 22 */
                {3,1,5,7, 6,4,0,2 }  /* orientation 23 */
              },
              /* orientation Table */
              {
                {1 ,2 ,0 ,3 , 4 ,0 ,5 ,6 },
                {0 ,7 ,1 ,8 , 5 ,1 ,4 ,9 },
                {15,0 ,2 ,22, 20,2 ,19,23},
                {20,6 ,3 ,23, 15,3 ,16,22},
                {22,13,4 ,12, 11,4 ,1 ,20},
                {11,19,5 ,20, 22,5 ,0 ,12},
                {9 ,3 ,6 ,2 , 21,6 ,17,0 },
                {10,1 ,7 ,11, 12,7 ,13,14},
                {12,9 ,8 ,14, 10,8 ,18,11},
                {6 ,8 ,9 ,7 , 17,9 ,21,1 },
                {7 ,15,10,16, 13,10,12,17},
                {5 ,14,11,9 , 0 ,11,22,8 },
                {8 ,20,12,19, 18,12,10,5 },
                {18,4 ,13,5 , 8 ,13,7 ,19},
                {17,11,14,1 , 6 ,14,23,7 },
                {2 ,10,15,18, 19,15,20,21},
                {19,17,16,21, 2 ,16,3 ,18},
                {14,16,17,15, 23,17,6 ,10},
                {13,21,18,17, 7 ,18,8 ,16},
                {16,5 ,19,4 , 3 ,19,2 ,13},
                {3 ,12,20,13, 16,20,15,4 },
                {23,18,21,10, 14,21,9 ,15},
                {4 ,23,22,6 , 1 ,22,11,3 },
                {21,22,23,0 , 9 ,23,14,2 }
              }
            }
          };

/* used by begin() */
template<int T_DIM>
inline OctreeIterator<T_DIM>::OctreeIterator( void )
 : root(NULL), todo(), done(), curpos(-1), ordering(-1), traversallist() {}

template<int T_DIM>
OctreeIterator<T_DIM>::OctreeIterator( Node * proot, int p_ordering )
 : root(proot), todo(),done(), curpos(-1), ordering(p_ordering), traversallist()
{
    if ( ordering == Ordering::Rows ) {
        this->curpos = 0;
        int maxdepth = root->getMaxLevel() - root->getLevel();
        VecD index = 0.5;
        VecD t_pos = root->center - root->size/2 + index*root->size/pow(2,maxdepth);
    } else {
        ToDoData toBeStored = { /* child index */ 0, proot, /* orientation of root */ 0 };
        todo.push( toBeStored );
    }
}

template<int T_DIM>
inline OctreeIterator<T_DIM>::~OctreeIterator( void ) {}

/****************************** Copy Constructor ******************************/
template<int T_DIM>
inline OctreeIterator<T_DIM>::OctreeIterator( const OctreeIterator & src )
: root(src.root), todo(src.todo), done(src.done), curpos(src.curpos),
  ordering(src.ordering), traversallist(src.traversallist)
{}

/**************************** Assignment Operator *****************************/
template<int T_DIM>
inline OctreeIterator<T_DIM> & OctreeIterator<T_DIM>::operator=
( const OctreeIterator & src ) {
    this->root          = src.root;
    this->todo          = src.todo;
    this->done          = src.done;
    this->curpos        = src.curpos;
    this->ordering      = src.ordering;
    this->traversallist = src.traversallist;
    return *this;
}

/***************************** Increment Operator *****************************/
template<int T_DIM>
OctreeIterator<T_DIM>& OctreeIterator<T_DIM>::operator++( void ) {
    if ( this->ordering == Ordering::Rows ) {
        this->curpos++;
        if ( this->curpos >= int(traversallist.size()) )
            this->curpos = -1;
        return *this;
    }
/* The iterator position is defined by what is currently at the top of the    *
 * todo-stack. Therefore we try to add exactly one element to the stack top.  *
 * If the stack top is a leaf, then that is the last iteration we were at, so *
 * pop it and look at next element, which should be a parent. That parent     *
 * also has save an ichild index specifying how much children we already      *
 * pushed to the todo-stack from this parent. If all children already were    *
 * pushed, then pop this parent. When pushing children no difference is being *
 * whether the child is a leaf or a parent !                                  */
    while( !todo.empty() ) {
        const Node * currentNode = todo.top().node;
        if ( currentNode->IsLeaf() ) {
            todo.pop();
        } else if ( todo.top().ichild < todo.top().node->nchildren ) {
/* The first child being pushed is the root with an initial orientation of 0. *
 * From that orientation we can derive the orientation and ordering for the   *
 * childs. So instead of accessing children[0],children[1],... we access:     *
 *   children[ orderingTable.Hilbert.ordering[0 <- root orientation ][0]      *
 *   children[ orderingTable.Hilbert.ordering[0 <- root orientation ][1]      *
 * The orientation of that child 0,1,... then is given by:                    *
 *   children[ orderingTable.Hilbert.orientation[0][0,1,...]                  */
            int childIndex       = -1;
            int childOrientation = -1;
            switch ( this->ordering ) {
                case Ordering::Morton:
                    childIndex = todo.top().ichild;
                    childOrientation = todo.top().orientation;
                    break;
                case Ordering::GrayCode:
                    childIndex = this->orderingTable.GrayCode.ordering
                        [todo.top().orientation][todo.top().ichild];
                    childOrientation = orderingTable.GrayCode.orientation
                        [todo.top().orientation][todo.top().ichild];
                    break;
                case Ordering::Hilbert:
                    childIndex = this->orderingTable.Hilbert.ordering
                        [todo.top().orientation][todo.top().ichild];
                    childOrientation = orderingTable.Hilbert.orientation
                        [todo.top().orientation][todo.top().ichild];
                    break;
            }
            const Node * childPtr = currentNode->getChildPtr( childIndex );
            assert( childPtr != NULL );
            /* child index we currently are at is first to begin with: 0 */
            ToDoData toBeStored = { 0, const_cast<Node*>( childPtr ),
                                    childOrientation                 };
            todo.top().ichild++;
            todo.push(toBeStored);
            return *this;
        } else /* all children of parents already pushed once */
            todo.pop();
    }
    return *this;
}

template<int T_DIM>
inline OctreeIterator<T_DIM> OctreeIterator<T_DIM>::operator++( int ) {
    OctreeIterator tmp( *this );
	++(*this);
	return tmp;
}

template<int T_DIM>
inline bool OctreeIterator<T_DIM>::operator==( const OctreeIterator & src ) {
    /* bool allequal = true;
    allequal = allequal and ( this->todo     == src.todo     );
    allequal = allequal and ( this->done     == src.done     );
    allequal = allequal and ( this->curpos   == src.curpos   );
    allequal = allequal and ( this->ordering == src.ordering );
    return allequal; */
    return not (*this != src);
}

template<int T_DIM>
inline bool OctreeIterator<T_DIM>::operator!=( const OctreeIterator & src ) {
    if (ordering == Ordering::Rows)
        return this->curpos != src.curpos;
    return this->todo.size() != src.todo.size() or
           this->curpos != src.curpos;
}

template<int T_DIM>
inline typename OctreeIterator<T_DIM>::Node & OctreeIterator<T_DIM>::operator*( void ) const {
    if ( this->ordering == Ordering::Rows ) {
        return *(traversallist[curpos]);
    } else
        return *(this->todo.top().node);
}

template<int T_DIM>
inline typename OctreeIterator<T_DIM>::Node * OctreeIterator<T_DIM>::operator->( void ) const {
    if ( this->ordering == Ordering::Rows ) {
        #ifndef NDEBUG
            bool toassert = curpos < int(traversallist.size());
            if ( !toassert )
                tout << "curpos=" << curpos << " >= " << traversallist.size()
                     << "=traversallist.size()\n";
            assert( curpos < int(traversallist.size()) );
        #endif
        assert( curpos < int(traversallist.size()) );
        return traversallist[curpos];
    } else
        return this->todo.top().node;
}

/* return iterator with empty stack */
template<int T_DIM>
inline OctreeIterator<T_DIM> OctreeIterator<T_DIM>::begin( void )
{
    OctreeIterator<T_DIM> tmp = *this;
    if ( tmp.traversallist.size() == 0 and tmp.ordering == Ordering::Rows ) {
        int maxdepth   = tmp.root->getMaxLevel() - tmp.root->getLevel();
        VecI t_size    = int(pow(2,maxdepth));
        #if DEBUG_OCTREE_ITERATOR > 0
            int nleaves    = tmp.root->countLeaves();
            int ileave    = 0;
            clock_t tstart = clock();
            clock_t ttmp   = clock();
        #endif

        for ( int i=0; i < t_size.product(); i++ ) {
            VecD index = VecD(ConvertLinearToVectorIndex<T_DIM>(i,t_size)) + 0.5;
            Node * found = tmp.root->FindLeafContainingPos( tmp.root->center -
                tmp.root->size/2 + index * tmp.root->size/pow(2,maxdepth) );
            /* Search in list to avoid traversing one node two times */
            if ( std::find( tmp.done.begin(), tmp.done.end(), found) == tmp.done.end() ) {
                tmp.done.push_back(found);
                tmp.traversallist.push_back(found);
                tmp.curpos = i;
                #if DEBUG_OCTREE_ITERATOR > 0
                    if ( found->IsLeaf() )
                            ileave++;
                #endif
            }
            #if DEBUG_OCTREE_ITERATOR > 0
                if ( clock() - ttmp > CLOCKS_PER_SEC ) {
                    ttmp = clock();
                    tout << "[" << (double) 100.0 * ileave / nleaves
                         << "%] Currently at node at " << tmp.done.back()->center << "\n";
                }
            #endif
        }
        #if DEBUG_OCTREE_ITERATOR > 0
            tout << "Iterator took " << clock() - tstart << " in total.\\n";
        #endif
    }
    if ( tmp.traversallist.size() > 0 )
        tmp.curpos = 0;
    else
        tmp.curpos = -1;
    return tmp;
}

/* return iterator with empty stack */
template<int T_DIM>
inline OctreeIterator<T_DIM> OctreeIterator<T_DIM>::end( void ) const {
    OctreeIterator it;
    return it;
}

} // namespace Octree
