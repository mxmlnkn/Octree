#pragma once

#define DEBUG_OCTREE_ITERATOR 0


namespace Octree {

template<int T_DIM>
typename Node<T_DIM>::iterator::orderingTableStruct
    Node<T_DIM>::iterator::orderingTable =
          {
            /* Gray-Code ordering */
            { /* ordering Table */
              { /* Row 0 meaning parent has orientation 0 */
                /* Columns. 0th column contains index and orientation */
                { 0,1,3,2 },
                { 3,2,0,1 }
              },
              /* orientation Table */
              { /* Row 0 meaning parent has orientation 0 */
                /* Columns. 0th column contains index and orientation */
                { 0,1,1,0 },
                { 1,0,0,1 }
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

/* used by begin() */
template<int T_DIM>
Node<T_DIM>::iterator::iterator( void )
 : todo(), done(), curpos(-1), ordering(-1) {}

template<int T_DIM>
Node<T_DIM>::iterator::iterator( Node * root, int p_ordering )
 : todo(), done(), curpos(-1), ordering(p_ordering)
{
    if ( ordering == Ordering::Rows ) {
        this->curpos = 0;
        int maxdepth = root->getMaxLevel()-root->getLevel();
        VecD index = 0.5;
        VecD t_pos = root->center - root->size/2 + index*root->size/pow(2,maxdepth);
        done.push_back( root );
        done.push_back( root->FindLeafContainingPos( t_pos ) );
    } else {
        tododata toBeStored = { /* child index */ 0, root, /* orientation of root */ 0 };
        todo.push( toBeStored );
    }
}

template<int T_DIM>
Node<T_DIM>::iterator::~iterator( void ) {}

/****************************** Copy Constructor ******************************/
template<int T_DIM>
Node<T_DIM>::iterator::iterator( const iterator & src ) 
: todo(src.todo), done(src.done), curpos(src.curpos), ordering(src.ordering)
{}

/**************************** Assignment Operator *****************************/
template<int T_DIM>
typename Node<T_DIM>::iterator & Node<T_DIM>::iterator::operator=
( const iterator & src ) {
    this->todo     = src.todo;
    this->done     = src.done;
    this->curpos   = src.curpos;
    this->ordering = src.ordering;
    return *this;
}

/***************************** Increment Operator *****************************/
template<int T_DIM>
typename Node<T_DIM>::iterator& Node<T_DIM>::iterator::operator++( void ) {
    assert( T_DIM == 2 or ordering == Ordering::Morton );
    if ( this->ordering == Ordering::Rows ) {
        Node * root  = done.front();
        int maxdepth = root->getMaxLevel() - root->getLevel();
        VecI t_size  = int(pow(2,maxdepth));
        int i = this->curpos + 1;
        /* set to end(), if nothing valid found it won't be changed back */
        this->curpos = -1;
        for ( ; i < t_size.product(); i++ ) {
            VecD index = VecD(ConvertLinearToVectorIndex<T_DIM>(i,t_size)) + 0.5;
            Node * found = root->FindLeafContainingPos( root->center - 
                           root->size/2 + index*root->size/pow(2,maxdepth) );
            if ( std::find( done.begin(), done.end(), found) == done.end() ) {
                done.push_back(found);
                this->curpos = i;
                break;
            }
        }
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
            tododata toBeStored = { 0, const_cast<Node*>( childPtr ),
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
typename Node<T_DIM>::iterator Node<T_DIM>::iterator::operator++( int ) {
    iterator tmp( *this );
	++(*this);
	return tmp;
}

template<int T_DIM>
bool Node<T_DIM>::iterator::operator==( const iterator & src ) {
    /* bool allequal = true;
    allequal = allequal and ( this->todo     == src.todo     );
    allequal = allequal and ( this->done     == src.done     );
    allequal = allequal and ( this->curpos   == src.curpos   );
    allequal = allequal and ( this->ordering == src.ordering );
    return allequal; */
    return not (*this != src);
}

template<int T_DIM>
bool Node<T_DIM>::iterator::operator!=( const iterator & src ) {
    if (ordering == Ordering::Rows) {
    }
    return this->todo.size() != src.todo.size() or 
           this->curpos != src.curpos;
}

template<int T_DIM>
typename Node<T_DIM>::Node & Node<T_DIM>::iterator::operator*( void ) const {
    if ( this->ordering == Ordering::Rows )
        return *(this->done.back());
    else
        return *(this->todo.top().node);
}

template<int T_DIM>
typename Node<T_DIM>::Node * Node<T_DIM>::iterator::operator->( void ) const {
    if ( this->ordering == Ordering::Rows )
        return this->done.back();
    else
        return this->todo.top().node;
}

/* return iterator with empty stack */
template<int T_DIM>
typename Node<T_DIM>::iterator Node<T_DIM>::iterator::end( void ) {
    iterator it;
    return it;
}

} // namespace Octree
