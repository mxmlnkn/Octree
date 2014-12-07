#pragma once

#define DEBUG_OCTREE_NODE 8




namespace Octree {

/* This should only be used in conjunction with operator= !!                  */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::Node( void )
: parent( NULL ), center( VecD(0) ), size( VecD(0) ) {
    for (int i=0; i<nchildren; i++)
        this->children[i] = NULL;
}

/* Copy Constructor */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::Node(const Node & src) {
    assert(false);
}

/* Copy Assigment (Problematic because of const members :( )                  */
template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::Node& Node<T_DTYPE,T_DIM>::operator=(const Node & src)
{
    this->parent = NULL;
    this->center = src.center;
    this->size   = src.size;
    for (int i=0; i<this->nchildren; ++i) {
        if (this->children[i] != NULL)
            delete this->children[i];
        if (src.children[i] != NULL) {
              this->children[i]  = new class Node;
            *(this->children[i]) = *(src.children[i]);
              this->children[i]->parent = this;
        } else
            this->children[i] = NULL;
    }
    /* it breaks here ? */
    this->data = src.data; /* <list> has copy assignment operator */
    return *this;
}

/* Destructor */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::~Node( void ) {
    for (int i=0; i<nchildren; i++)
        if ( this->children[i] != NULL )
            delete this->children[i];
}

/* Constructor */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::Node(Node * const parent, VecD const cent, VecD const size)
: parent( parent ), center( cent ), size( size ) {
    for (int i=0; i<nchildren; i++)
        this->children[i] = NULL;
}

template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::IsInside ( const VecD & pos ) const {
    return (pos >= (this->center - 0.5*this->size)) and
           (pos <  (this->center + 0.5*this->size));
}

template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::IsLeaf( void ) const {
/* We assume the correctness of this class, i.e. all children are NULL or all *
 * children are not NULL                                                      */
    return children[0] == NULL;
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::VecI Node<T_DTYPE,T_DIM>::ConvertNumberToDirection
( const int number ) const {
    VecI direction;
    int tmp = number;
    for (int i=T_DIM-1; i>=0; --i) {
        direction[i] = ( (tmp & 1) == 0 ) ? -1 : 1;
        tmp = tmp >> 1;
    }
#if DEBUG_OCTREE_NODE >= 11
    std::cerr << "number: " << number << " => " << direction << std::endl;
#endif
    //assert( ConvertDirectionToNumber(direction) == number );
    return direction;
}

template<typename T_DTYPE, int T_DIM>
int Node<T_DTYPE,T_DIM>::ConvertDirectionToNumber( const VecI direction ) const {
    int tmp = 0;
    for (int i=0; i<T_DIM; ++i) {
        tmp = ( tmp << 1) | (1+direction[i])/2;
    }
#if DEBUG_OCTREE_NODE >= 11
    std::cerr << "direction: " << direction << " => " << tmp << std::endl;
#endif
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::Node * Node<T_DTYPE,T_DIM>::FindLeafContainingPos
( const VecD & pos ) {
/* Prone to rounding errors :(, except if worldsize is e.g. 1 and center is   *
 * 0, because in that case all sizes and new centers will be in the form of   *
 * 1/2^N which can be represented exactly with floating points!               *
 *   => Use this internally and overlay it with user-sizes !                  */
    bool insideNode = this->IsInside(pos);
    if (!insideNode) {
        std::cout << "pos: " << pos << " center: " << this->center << " size: "
                  << this->size << std::endl;
        assert( insideNode );
    }
    Node * tmp = this;
#if DEBUG_OCTREE_NODE >= 10
        std::cerr << "Try to find correct leaf for " << pos << std::endl;
#endif
    while( ! tmp->IsLeaf() ) {
/* E.g. center is (0,0) then pos =(0.3,-0.1) will with Greater() result in    *
 * (1,0) which in turn will be converted to (1,-1).                           */
        VecD direction = 2*VecI( pos.GreaterThan( tmp->center )) - 1;
        tmp = tmp->children[ ConvertDirectionToNumber( direction ) ];
#if DEBUG_OCTREE_NODE >= 10
        std::cerr << "Comparison with center " << tmp->center << " resulted in "
                  << pos.GreaterThan( tmp->center ) << " -> "
                  << 2*VecI( pos.GreaterThan( tmp->center )) - 1 << std::endl;
#endif
    }
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
void Node<T_DTYPE,T_DIM>::GrowUp( void )
{
#if DEBUG_OCTREE_NODE >= 10
    std::cerr << "Leaf at " << this->center << " growing Up!\n";
#endif
    for (int i=0; i<compileTime::pow(2,T_DIM); ++i) {
        VecD direction = VecD( ConvertNumberToDirection( i ) );
        this->children[i] = new class Node( this, this->center +
                            0.25*direction*size, 0.5*size );
    }
/* Migrate data to child-leaf nodes. If it is too much data they will recursi-*
 * vely grow up to (done by InsertData), which is not nice, but it shouldn't  *
 * happen too often, because we only added one element. For above to happen   *
 * all elements need to be in the same next-gen octant. If no data in this    *
 * leaf, then data.begin() != data.end(), therefore the while loop won't be   *
 * entered                                                                    */
    typename Datalist::iterator it = this->data.begin();
    while( it != this->data.end() ) {
        this->FindLeafContainingPos(it->pos)->InsertData( it->pos, it->object );
        this->data.erase( it++ );
    }
    assert( this->data.empty() );
#if DEBUG_OCTREE_NODE >= 10
    std::cerr << "All Grown Up!\n";
#endif
}

template<typename T_DTYPE, int T_DIM>
void Node<T_DTYPE,T_DIM>::InsertData( const VecD pos, T_DTYPE * const object )
{
    assert( this->IsLeaf() );
    Datapoint packed = { pos, object };
/* Push back no matter of the final size. It will be split up when splitting  *
 * the octant to 8 new octants anyway                                         */
    this->data.push_back( packed );
    if ( data.size() > this->maxdata )
        this->GrowUp();
}


template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::RemoveData( const VecD pos, T_DTYPE * const object )
{
    assert( this->IsLeaf() );
    typename Node::Datalist::iterator it = this->data.begin();
    while ( it != this->data.end() ) {
        if ( it->object == object )
            break;
        ++it;
    }
    if ( it == this->data.end() )
        return false;
    this->data.erase(it);
    /* Don't try this (at home) with root leaf */
    if (this->parent != NULL)
        this->parent->Rejuvenate();
    return true;
}


template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::HasOnlyLeaves( void ) const {
    if ( this->IsLeaf() )
        return false;
    else
        for ( int i=0; i < this->nchildren; ++i )
            if ( this->children[i] != NULL )
                if ( not this->children[i]->IsLeaf() )
                    return false;
    return true;
}

template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::Rejuvenate( void )
{
    assert( this->IsLeaf() == false );
/* Now check if we can collapse the parent node after removing that element.  *
 * It won't be possible, if the parent node has one or more non-leaf child    *
 * nodes, because that child node will already too much elements for the      *
 * parent node to store. This functionality could be enforced to be called by *
 * the in order to save time                                                  */
    if ( this->HasOnlyLeaves() ) {
        int sumelements = 0;
        for ( int i=0; i < this->nchildren; ++i )
            if ( children[i] != NULL ) {
                #if DEBUG_OCTREE_NODE >= 9
                std::cerr << "Child[" << i << "] = " << children[i]->data.size()
                          << std::endl;
                #endif
                sumelements += children[i]->data.size();
            }
        #if DEBUG_OCTREE_NODE >= 9
        std::cerr << "Rejuvenating Parent at " << center << " sized "
                  << size << " if " << sumelements << " =< "
                  << maxdata << std::endl;
        #endif
        if ( sumelements <= maxdata ) {
/* If we can collapse, then do so by copying all data into this node.         */
            assert( this->data.empty() );
            for ( int i=0; i < nchildren; ++i )
                if (  children[i] != this
                  and children[i] != NULL ) {
                    this->data.splice( this->data.end(),
                                       this->children[i]->data );
                    assert( this->children[i]->data.empty() );
                    delete children[i];
                    children[i] = NULL;
                }
            return true;
        }
    }
    return false;
}

template<typename T_DTYPE, int T_DIM>
const typename Node<T_DTYPE,T_DIM>::Node * Node<T_DTYPE,T_DIM>::getChildPtr( const int i ) const
{
    if ( i < this->nchildren )
        return children[i];
    else
        return NULL;
}


template<typename T_DTYPE, int T_DIM>
const typename Node<T_DTYPE,T_DIM>::Datapoint * Node<T_DTYPE,T_DIM>::getDataPtr( const int i ) const {
    typename Node::Datalist::const_iterator it = data.begin();
    if ( it == data.end() )
        return NULL;
    for (int j=0; j<i; ++j) {
        ++it;
        if ( it == data.end() )
            return NULL;
    }
    return &(*it);
}

template<typename T_DTYPE, int T_DIM>
int Node<T_DTYPE,T_DIM>::getLevel( void ) const {
    Node * tmp = this->parent;
    int level  = 0;
    while (tmp != NULL) {
        tmp = tmp->parent;
        level++;
    }
    return level;
}


} // namespace Octree

