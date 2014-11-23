#pragma once

#define DEBUG_OCTREE 0


namespace Octree {


/* Should only be used in conjunction with operator= !! */
template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::Octree( void )
: center( VecD(0) ), size( VecD(0) ), root( NULL ) {}

/* Copy Constructor */
template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::Octree( const Octree & src ) {
    this->center  = src.center;
    this->size    = src.size;
    this->root    = new Node( NULL, VecD(0), VecD(1) );
    *(this->root) = *(src.root);
}

template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::~Octree( void ) {
    delete this->root;
}

/* Copy Assignment */
template<typename T_DTYPE, int T_DIM>
typename Octree<T_DTYPE,T_DIM>::Octree& Octree<T_DTYPE,T_DIM>::operator=(const Octree & src)
{
    this->center = src.center;
    this->size   = src.size;
    if (this->root != NULL) {
        delete this->root;
        this->root = new Node( NULL, VecD(0), VecD(1) );
    }
    *(this->root) = *(src.root);
    this->root->parent = NULL;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::Octree( const VecD center, const VecD size )
: center( center ), size( size ) {
/* root Node has no parent, therefore it's parent is set to NULL !            */
    this->root = new Node( NULL, VecD(0), VecD(1) );
}

template<typename T_DTYPE, int T_DIM>
const typename Octree<T_DTYPE,T_DIM>::Node * Octree<T_DTYPE,T_DIM>::GetNodePtr
( const VecD center ) const
{
    const Node * tmp = this->root;
    while( ! tmp->IsLeaf() ) {
        if ( tmp->center == center )
            return tmp;
        VecI direction = 2*VecI( center.GreaterThan( tmp->center )) - 1;
        tmp = tmp->children[ tmp->ConvertDirectionToNumber( direction ) ];
    }
    if ( tmp->center == center )
        return tmp;
    else
        return NULL;
}

/* Very slow, shouldn't be used inside a simulation ! In order to traverse    *
 * this "non-recursively" we simulate a stack with <stack>. In the stack the  *
 * next todo Node combined with the next to do child-number is being stored.  *
 * The child number will be incremented after completely traversing that      *
 * child. If child number is higher than can be indexed we are at the end and *
 * that corresponding node can be popped from the stack                       */
template<typename T_DTYPE, int T_DIM>
typename Octree<T_DTYPE,T_DIM>::VecD Octree<T_DTYPE,T_DIM>::FindData
( T_DTYPE * const object ) const
{
    typedef struct{ int ichild; const Node* node; } tododata;
    std::stack<tododata> todo;
    tododata tmp = { 0, this->root };
    todo.push( tmp );
    while( !todo.empty() ) {
        const Node * current = todo.top().node;
        if ( current->IsLeaf() ) {
            typename Node::Datalist::const_iterator it = current->data.begin();
            while ( it != current->data.end() ) {
                if ( it->object == object )
                    return this->size * it->pos;
                ++it;
            }
            todo.pop();
        } else {
            if ( todo.top().ichild < current->nchildren ) {
                tmp.node   = current->children[ todo.top().ichild++ ];
                tmp.ichild = 0;
                todo.push(tmp);
            } else
                todo.pop();
        }
    }
    return VecD( NAN );
}

template<typename T_DTYPE, int T_DIM>
void Octree<T_DTYPE,T_DIM>::InsertData( const VecD pos, T_DTYPE * const object ) {
    this->root->FindLeafContainingPos( pos / this->size )->
                InsertData( pos / this->size, object );
}

template<typename T_DTYPE, int T_DIM>
bool Octree<T_DTYPE,T_DIM>::RemoveData( const VecD pos, T_DTYPE * const object )
{
    return this->root->FindLeafContainingPos( pos / this->size )->RemoveData( pos / this->size, object );
}

template<typename T_DTYPE, int T_DIM>
bool Octree<T_DTYPE,T_DIM>::MoveData( const VecD pos, T_DTYPE * const object,
const VecD newpos ) {
    Node * tmp = this->root->FindLeafContainingPos( pos / this->size );
/* If the new position is still inside the same octant, then we don't need to *
 * relocate it and we can simply search for the correct data in the list and  *
 * change the position                                                        */
    typename Node::Datalist::iterator it = tmp->data.begin();
    while ( it != tmp->data.end() ) {
        if ( it->object == object )
            break;
        ++it;
    }
    if (it == tmp->data.end())
        return false;
    if ( tmp->IsInside( newpos / this->size ) )
        it->pos = newpos / this->size;
/* If the new position will change the parent leaf of that object, then it's  *
 * being relocated by removing and inserting it anew (could be optimized)     */
    else {
        RemoveData( pos   , object );
        InsertData( newpos, object );
    }
    return true;
}


template<typename T_DTYPE, int T_DIM>
bool Octree<T_DTYPE,T_DIM>::CheckIntegrity( void ) {
    bool foundError = false;
    typename Octree::Octree<T_DTYPE,T_DIM>::iterator it = this->begin();
    while( it != this->end() ) {
        if ( it->IsLeaf() ) {
            typename Node::Datalist::iterator dataIt = it->data.begin();
            while ( dataIt != it->data.end() ) {
                if ( ! it->IsInside( dataIt->pos ) ) {
                    foundError = true;
                    std::cerr << "Node at " << this->size * it->center << " sized " << this->size * it->size << " contains child it can't manage: " << dataIt->pos << " -> " << *(dataIt->object) << " !\n";
                }
                ++dataIt;
            }
        }
        
        /* Count children and check their parents */
        int nchilds = 0;
        for (int i=0; i < it->nchildren; ++i) {
            if (it->children[i] != NULL) {
                nchilds++;
                if ( it->children[i]->parent != &(*it) ) {
                    foundError = true;
                    std::cerr << "Parent at " << this->size * it->center << " kidnapped the child at " << it->children[i]->center << "!\n";
                }
            }
        }
            
        if (nchilds > 0 and it->IsLeaf()) {
            foundError = true;
            std::cerr << "Node at " << this->size * it->center << " sized " << this->size * it->size << " contains childs even though it is a leaf node !\n";
        }
        if (nchilds < it->nchildren and !it->IsLeaf()) {
            foundError = true;
            std::cerr << "Parent at " << this->size * it->center << " sized " << this->size * it->size << " lost one of its children !\n";
        }
        ++it;
    }
    return !foundError;
}



} // namespace Octree

/******************************** Debug output ********************************/

template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out, Octree::Octree<T_DTYPE,T_DIM>& tree ) {
    typedef typename Octree::Octree<T_DTYPE,T_DIM>::Node Node;
    typename Octree::Octree<T_DTYPE,T_DIM>::iterator it = tree.begin();
    while( it != tree.end() ) {
        if ( it->IsLeaf() ) {
            for ( int i=0; i < it->getLevel(); ++i )
                out << "  ";
            out << "Leaf at "    << tree.size * it->center
                << " with size " << tree.size * it->size << std::endl;
            typename Node::Datalist::iterator dataIt = it->data.begin();
            while ( dataIt != it->data.end() ) {
                for ( int i=0; i < it->getLevel()+1; ++i )
                    out << "  ";
                out << tree.size * dataIt->pos << " -> " << *(dataIt->object) << "\n";
                ++dataIt;
            }
        } else {
            for ( int i=0; i < it->getLevel(); ++i )
                out << "  ";
            out << "Parent at "  << tree.size * it->center
                << " with size " << tree.size * it->size << std::endl;
        }
        ++it;
    }
    return out;
}

