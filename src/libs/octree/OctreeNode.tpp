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

/************************* ConvertNumberToDirection ***************************/
template<int T_DIM>
inline typename Node<T_DIM>::VecI Node<T_DIM>::ConvertNumberToDirection( const int number ) const {
    VecI direction;
    int tmp = number;
    for (int i=0; i<T_DIM; ++i) {
        direction[i] = ( (tmp & 1) == 0 ) ? +1 : -1;
        tmp = tmp >> 1;
    }
    return direction;
}

template<int T_DIM>
inline int Node<T_DIM>::ConvertDirectionToNumber( const VecI direction ) {
    int tmp = 0;
    for (int i=T_DIM-1; i>=0; --i) {
        if ( direction[i]!=-1 and direction[i]!=+1 )
            return -1;
        tmp = ( tmp << 1) | (1-direction[i])/2;
    }
    return tmp;
}

/*************************** FindLeafContainingPos ****************************/
template<int T_DIM>
inline Node<T_DIM> * Node<T_DIM>::FindLeafContainingPos( const VecD & pos ) {
    if ( not this->IsInside(pos) )
        return NULL;
/* Prone to rounding errors :(, except if worldsize is e.g. 1 and center is   *
 * 0, because in that case all sizes and new centers will be in the form of   *
 * 1/2^N which can be represented exactly with floating points!               *
 *   => Use this internally and overlay it with user-sizes !                  */
    Node * tmp = this;
    while( ! tmp->IsLeaf() ) {
/* E.g. center is (0,0) then pos =(0.3,-0.1) will with Greater() result in    *
 * (1,0) which in turn will be converted to (1,-1).                           */
        VecD direction = 2*VecI( pos.GreaterOrEqualThan( tmp->center )) - 1;
        tmp = tmp->children[ ConvertDirectionToNumber( direction ) ];
    }
    assert( tmp->IsLeaf() );
    assert( tmp->IsInside(pos) );
    return tmp;
}

template<int T_DIM>
inline Node<T_DIM> const * Node<T_DIM>::FindLeafContainingPos( const VecD & pos ) const {
    return const_cast<Node<T_DIM> const *>( const_cast<Node<T_DIM>*>(this)->FindLeafContainingPos( pos ) );
}

/*********************************** begin ************************************/
template<int T_DIM>
inline typename Node<T_DIM>::iterator Node<T_DIM>::begin( int ordering ) {
    iterator it( this, ordering );
    return it.begin();
}

template<int T_DIM>
inline typename Node<T_DIM>::iterator Node<T_DIM>::end( void ) {
    iterator it;
    return it.end();
}

/******************************************************************************/

template<int T_DIM>
inline const typename Node<T_DIM>::Node * Node<T_DIM>::getChildPtr( const int i ) const
{
    if ( i < this->nchildren )
        return children[i];
    else
        return NULL;
}

template<int T_DIM>
inline typename Node<T_DIM>::VecD Node<T_DIM>::getSize( void ) const {
    return this->size;
}

template<int T_DIM>
inline bool Node<T_DIM>::IsInside ( const VecD & pos ) const {
    return (pos >= (this->center - 0.5*this->size)) and
           (pos <  (this->center + 0.5*this->size));
}

/********************************** IsLeaf ************************************/
template<int T_DIM>
inline bool Node<T_DIM>::IsLeaf( void ) const {
/* We assume the correctness of this class, i.e. all children are NULL or all *
 * children are not NULL                                                      */
    return children[0] == NULL;
}

template<int T_DIM>
inline int Node<T_DIM>::getLevel( void ) const {
    Node * tmp = this->parent;
    int level  = 0;
    while (tmp != NULL) {
        tmp = tmp->parent;
        level++;
    }
    return level;
}

/*********************************** GrowUp ***********************************/
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
                sumelements += (int) children[i]->data.size();
            }
        /* If we can collapse, then do so by copying all data into this node. */
        if ( sumelements <= maxdata ) {
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
inline bool Node<T_DIM>::HasOnlyLeaves( void ) const {
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
void Node<T_DIM>::DeleteChildren( void )
{
    for ( int i=0; i < nchildren; ++i )
        if ( children[i] != NULL ) {
            delete children[i];
            children[i] = NULL;
        }
}

/******************************** getMinLevel *********************************/
template<int T_DIM>
int Node<T_DIM>::getMinLevel( void ) {
    int minLevel = 255;
    for ( iterator it = this->begin(); it != it.end(); ++it )
        if ( it->IsLeaf() ) {
            int curLevel = it->getLevel();
            assert( curLevel <= 255 );
            if ( minLevel > curLevel )
                minLevel = curLevel;
        }
    return minLevel;
}

/******************************** getMaxLevel *********************************/
template<int T_DIM>
int Node<T_DIM>::getMaxLevel( void ) {
    int maxLevel = 0;
    for ( iterator it = this->begin(); it != it.end(); ++it )
        if ( it->IsLeaf() ) {
            int curLevel = it->getLevel();
            assert( curLevel >= 0 );
            if ( maxLevel < curLevel )
                maxLevel = curLevel;
        }
    return maxLevel;
}

template<int T_DIM>
inline int Node<T_DIM>::countLeaves( void ) {
    int sum = 0;
    iterator it( this, 0 );
    while( it != it.end() ) {
        if ( (it++)->IsLeaf() )
            sum +=1;
    }
    return sum;
}

/******************************** getNeighbor *********************************/
template<int T_DIM>
inline Node<T_DIM> * Node<T_DIM>::getNeighbor
( const VecI targetDir, const VecI periodic )
{
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

/******************************** getNeighbors ********************************/
template<int T_DIM>
typename std::list<Node<T_DIM>*> Node<T_DIM>::getNeighbors
( const VecI targetDir, const VecI periodic )
{
/* Neighbors may be not leaves (e.g. right neighbor of 3 is is parent of 1).  *
 * In that case iterate over it's children. 4 of 1 will with this +-----+-+-+ *
 * strategy iterate only one time, because  it2.begin() is a      |     |4|1| *
 * leaf. To check whether the node or its children are direct     |  5  +-+-+ *
 * neighbors, we test if getNeighbor in the opposite direction    |     |3|2| *
 * returns the current node (it). Because this only works from    +-----+-+-+ *
 * small to big, also test the other way around. E.g. left of 4 is 5, but     *
 * right of 5 will yield parent of 4                                          */
    std::list<Node<T_DIM>*> neighbors;
    VecI opdir = getOppositeDirection<T_DIM>( targetDir );
    Node* neighbor = this->getNeighbor(targetDir, periodic);
    if ( neighbor != NULL )
        for ( iterator it = neighbor->begin(); it != it.end(); ++it ) {
            Node* revneighbor = it->getNeighbor( opdir, periodic );
            assert( revneighbor != NULL );
            if ( revneighbor->IsInside( this->center ) )
                neighbors.push_back( &(*it) );
        }
    return neighbors;
}


} // namespace Octree

