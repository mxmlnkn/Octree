#pragma once

#define DEBUG_OCTREE_NODE 9




namespace Octree {

/****************************** Copy Constructor ******************************/
template<int T_DIM>
Node<T_DIM>::Node(const Node & src) {
    assert(false);
}

/********** Copy Assigment (Problematic because of const members :( ) *********/
template<int T_DIM>
typename Node<T_DIM>::Node& Node<T_DIM>::operator=(const Node & src) {
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

/********************************* Destructor *********************************/
template<int T_DIM>
Node<T_DIM>::~Node( void ) {
    for (int i=0; i<nchildren; i++)
        if ( this->children[i] != NULL )
            delete this->children[i];
}

/******************************** Constructor *********************************/
template<int T_DIM>
Node<T_DIM>::Node( void )
 : parent( NULL ), center( VecD(0) ), size( VecD(0) ), data(), 
   nchildren( compileTime::pow(2,T_DIM) )
{
    for (int i=0; i<nchildren; i++)
        this->children[i] = NULL;
}

template<int T_DIM>
Node<T_DIM>::Node(Node * const p_parent, VecD const p_center, VecD const p_size)
: parent( p_parent ), center( p_center ), size( p_size ), data(), 
  nchildren( compileTime::pow(2,T_DIM) )
{
    for (int i=0; i<nchildren; i++)
        this->children[i] = NULL;
}

template<int T_DIM>
bool Node<T_DIM>::IsInside ( const VecD & pos ) const {
    return (pos >= (this->center - 0.5*this->size)) and
           (pos <  (this->center + 0.5*this->size));
}

/********************************** IsLeaf ************************************/
template<int T_DIM>
bool Node<T_DIM>::IsLeaf( void ) const {
/* We assume the correctness of this class, i.e. all children are NULL or all *
 * children are not NULL                                                      */
    return children[0] == NULL;
}

template<int T_DIM>
typename Node<T_DIM>::VecI Node<T_DIM>::ConvertNumberToDirection( const int number ) const {
    VecI direction;
    int tmp = number;
    for (int i=0; i<T_DIM; ++i) {
        direction[i] = ( (tmp & 1) == 0 ) ? +1 : -1;
        tmp = tmp >> 1;
    }
#if DEBUG_OCTREE_NODE >= 11
    std::cerr << "number: " << number << " => " << direction << std::endl;
#endif
    //assert( ConvertDirectionToNumber(direction) == number );
    return direction;
}

template<int T_DIM>
int Node<T_DIM>::ConvertDirectionToNumber( const VecI direction ) {
    int tmp = 0;
    for (int i=T_DIM-1; i>=0; --i) {
        if ( direction[i]!=-1 and direction[i]!=+1 )
            return -1;
        tmp = ( tmp << 1) | (1-direction[i])/2;
    }
#if DEBUG_OCTREE_NODE >= 11
    std::cerr << "direction: " << direction << " => " << tmp << std::endl;
#endif
    return tmp;
}

template<int T_DIM> 
Node<T_DIM> * Node<T_DIM>::FindLeafContainingPos( const VecD & pos ) {
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
    assert( tmp->IsLeaf() );
    return tmp;
}

template<int T_DIM>
void Node<T_DIM>::GrowUp( void ) {
    for (int i=0; i<nchildren; ++i) {
        VecD direction = VecD( ConvertNumberToDirection( i ) );
        this->children[i] = new class Node( this, this->center +
                            0.25*direction*size, 0.5*size );
    }
    assert( this->data.empty() );
}

template<int T_DIM>
bool Node<T_DIM>::HasOnlyLeaves( void ) const {
    if ( this->IsLeaf() )
        return false;
    else
        for ( int i=0; i < this->nchildren; ++i )
            if ( this->children[i] != NULL )
                if ( not this->children[i]->IsLeaf() )
                    return false;
    return true;
}

template<int T_DIM>
bool Node<T_DIM>::Rejuvenate( void )
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

template<int T_DIM>
const typename Node<T_DIM>::Node * Node<T_DIM>::getChildPtr( const int i ) const
{
    if ( i < this->nchildren )
        return children[i];
    else
        return NULL;
}

template<int T_DIM>
typename Node<T_DIM>::VecD Node<T_DIM>::getSize( void ) const {
    return this->size;
}

template<int T_DIM>
int Node<T_DIM>::getLevel( void ) const {
    Node * tmp = this->parent;
    int level  = 0;
    while (tmp != NULL) {
        tmp = tmp->parent;
        level++;
    }
    return level;
}

template<int T_DIM>
int Node<T_DIM>::countLeaves( void ) {
    int sum = 0;
    iterator it = this->begin();
    while( it != this->end() ) {
        if ( (it++)->IsLeaf() )
            sum +=1;
    }
    return sum;
}

template<int T_DIM>
Node<T_DIM> * Node<T_DIM>::getNeighbor( const VecI targetDir, const VecI periodic ) {
    /* We want to find the neighbor of same size if possible, so go down as   *
     * much levels as we did have to go up. Depending on the direction the    *
     * target neighbor node is in and the position the Current Node is in the *
     * current parent, we have to go down                                     *
     *    +-----------+-----+-----+                                           *
     *    |           |     |     |                                           *
     *    |           |     |     |                                           *
     *    |           |     |     |                                           *
     *    | (0.5,0.5) |-----+-----+                                           *
     *    |           |     |     |                                           *
     *    |           |     |     |                                           *
     *    |           |     |     |                                           *
     *    +--+--+--+--+-----+-----+                                           *
     *    |  |  |  | P|     |     |                                           *
     *    +--+--+--+--+  T  |     |                                           *
     *    |  |  |  |  |     |     |                                           *
     *    +--+--+--+--+-----+-----+                                           *
     *    |  |  |  |  |     |     |                                           *
     *    +--+--+--+--+     |     |                                           *
     *    |  |  |  |  |     |     |                                           *
     *    +--+--+--+--+-----+-----+                                           */
    Node * root = this;
    while ( root->parent != NULL )
        root = root->parent;
    VecD theoreticalPosition = this->center + this->size * VecD(targetDir);
    /* if we couldn't find a neighbor, because we are e.g. at the border of   *
     * the octree global area (and we are not periodic (yet) ), then return   *
     * NULL                                                                   */
    for ( int i=0; i<T_DIM; ++i ) {
        if ( not periodic[i] )
            continue;
        double & t = theoreticalPosition[i];
        if ( t >  0.5*root->size[i] )
            t -= root->size[i];
        if ( t < -0.5*root->size[i] )
            t += root->size[i];
    }
    if ( periodic.sum() == T_DIM )
        assert( root->IsInside( theoreticalPosition) );
    if ( not root->IsInside( theoreticalPosition) )
        return NULL;
    /* a) try to neighbor anticipating it to be of the same size -> use       *
     *    FindleafContaining. This will get smaller or larger node            */
    Node * neighbor = root->FindLeafContainingPos( theoreticalPosition );
    /* b) As it's not very useful and deterministic to get a smaller node as  *
     *    a neighbor, we have to go up until we have a node of the same size  */
    while ( neighbor->size < this->size )
        neighbor = neighbor->parent;
    return neighbor;
}


} // namespace Octree

