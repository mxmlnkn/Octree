#pragma once

#define DEBUG_OCTREE 0


namespace Octree {

/****************************** Copy Constructor ******************************/
template<int T_DIM>
Octree<T_DIM>::Octree( const Octree & src )
: center(src.center), size(src.size), root( new Node( NULL, VecD(0), VecD(1) ) )
{
    *(this->root) = *(src.root);
}

template<int T_DIM>
Octree<T_DIM>::~Octree( void ) {
    delete this->root;
}

/****************************** Copy Assignment *******************************/
template<int T_DIM>
typename Octree<T_DIM>::Octree& Octree<T_DIM>::operator=(const Octree & src)
{
    this->center = src.center;
    this->size   = src.size;
    if (this->root != NULL)
        delete this->root;
    this->root = new Node( NULL, VecD(0), VecD(1) );
    /* copy root of src to newly allocated space of this object's root */
    *(this->root) = *(src.root);
    this->root->parent = NULL;
    return *this;
}

/********************************* Constructor ********************************/
template<int T_DIM>
Octree<T_DIM>::Octree( const VecD p_center, const VecD p_size )
: center( p_center ), size( p_size ), root( NULL ) {
/* root Node has no parent, therefore it's parent is set to NULL !            */
    this->root = new Node( NULL, VecD(0), VecD(1) );
}

template<int T_DIM>
Octree<T_DIM>::Octree( void )
: center( VecD(0) ), size( VecD(1) ), root( NULL ) {
    this->root = new Node( NULL, VecD(0), VecD(1) );
}


/********************************* GetNodePtr *********************************/
template<int T_DIM>
const typename Octree<T_DIM>::Node * Octree<T_DIM>::GetNodePtr
( const VecD p_center ) const
{
    Node * tmp = root;
    while( ! tmp->IsLeaf() )
    {
        if ( tmp->center == p_center )
            return tmp;
        VecD direction = 2*VecI( p_center.GreaterOrEqualThan( tmp->center )) - 1;
        tmp = tmp->children[ tmp->ConvertDirectionToNumber( direction ) ];
    }

    if ( tmp->center == p_center )
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
template<int T_DIM>
typename Octree<T_DIM>::VecD Octree<T_DIM>::FindData
( void * const object ) const {
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
                    return toGlobalCoords(it->pos);
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

/********************************* InsertData *********************************/
template<int T_DIM>
void Octree<T_DIM>::InsertData( const VecD pos, void * const object ) {
    const VecD intpos = toInternalCoords(pos);
    this->root->FindLeafContainingPos(intpos)->data.push_back(object);
}

/********************************* RemoveData *********************************/
template<int T_DIM>
bool Octree<T_DIM>::RemoveData( const VecD pos, const int index ) {
    typename Node::Datalist data = root->FindLeafContainingPos(toInternalCoords(pos))->data;
    if ( index >= data.size() )
        return false;
    data.erase( data.begin() + index );
    return true;
}

/*************************** FindLeafContainingPos ****************************/
template<int T_DIM>
typename Octree<T_DIM>::Node* Octree<T_DIM>::FindLeafContainingPos( const VecD & pos ) {
    return this->root->FindLeafContainingPos(toInternalCoords(pos));
}

template<int T_DIM>
int Octree<T_DIM>::getMinLevel( void ) {
    return root->getMinLevel();
}

template<int T_DIM>
int Octree<T_DIM>::getMaxLevel( void ) {
    return root->getMaxLevel();
}

/******************************* CheckIntegrity *******************************/
template<int T_DIM>
bool Octree<T_DIM>::CheckIntegrity( void ) const {
    bool foundError = false;
    typename Octree<T_DIM>::iterator it = this->begin();
    while( it != this->end() ) {
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

/*********************************** begin ************************************/
template<int T_DIM>
typename Octree<T_DIM>::Node::iterator Octree<T_DIM>::begin( int ordering ) const {
    typename Node::iterator it( this->root, ordering );
    return it.begin();
}

/************************************ end *************************************/
template<int T_DIM>
typename Octree<T_DIM>::Node::iterator Octree<T_DIM>::end( void ) const {
    typename Node::iterator it;
    return it.end();
}

template<int T_DIM>
typename Octree<T_DIM>::VecD Octree<T_DIM>::toInternalCoords
( const VecD pos ) const {
    return (pos - center) / size;
}

template<int T_DIM>
int Octree<T_DIM>::countLeaves(void) const {
    int leavesFound = 0;
    for ( typename Node::iterator it = this->begin(); it != this->end(); ++it )
        if ( it->IsLeaf() )
            ++leavesFound;
    return leavesFound;
}

template<int T_DIM>
typename Octree<T_DIM>::VecD Octree<T_DIM>::toGlobalCoords
( const VecD pos ) const {
    return center + pos*size;
}



} // namespace Octree

/******************************** Debug output ********************************/

template<int T_DIM>
std::ostream& operator<<( std::ostream& out, Octree::Octree<T_DIM>& tree ) {
    typedef typename Octree::Octree<T_DIM>::Node Node;
    typename Octree::Octree<T_DIM>::iterator it = tree.begin();
    while( it != tree.end() ) {
        if ( it->IsLeaf() ) {
            for ( int i=0; i < it->getLevel(); ++i )
                out << "  ";
            out << "Leaf at "    << tree.size * it->center
                << " with size " << tree.size * it->size << " contains pointers:\n";
            typename Node::Datalist::iterator dataIt = it->data.begin();
            while ( dataIt != it->data.end() ) {
                for ( int i=0; i < it->getLevel()+1; ++i )
                    out << "  ";
                out << *dataIt << "\n";
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

