/*
ToDo:
  - SVG Output
  - Coarsening
  - GetSize Function
  - later: parallelism
  - max depth
  - Find Error (Grow Up, Rejuvenate ausschließen ... )
  - Add CheckIntegrity Method, which checks for parents, child-pointers, ...
  - Test all methods singlehandedly
  + Bigger goal: parallelize ...
  - two kind of errors: one after 61 moves -> Octree.tpp
  - Weird boxes arrangement after timestep 2 -> Svg.tpp
*/

#pragma once

#include <cassert>
#include <iostream>
#include <cstdlib>   // malloc
#include <list>
#include <stack>
#include <ctime> // only for debug in operator<<
#include "Vector.h"

#define OCTREE_MAX_OBJECTS_PER_LEAF 3
/* Sets max level depth to 10, meaning at maximum 8^10 = 2^30 leaf nodes can  *
 * be created. Much more doesn't make sense anyway, because of memory and     *
 * because dx i.e. size becomes so small, that the centers of neighbors       *
 * decoded in  doubles won't differe from each other, meaning: center+size    *
 * =center. This happens at least at the sides (center ~ 1) (not in the       *
 * center at 0) for size == DBL_EPSILON = 2.22e-16. Size is halved at every   *
 * level, so size(level) = 1/2^level == 2.22e-16 => level = 52                */
#define OCTREE_MAX_DEPTH 10



namespace Octree {
template<typename T_DTYPE, int T_DIM>
class Octree;
}

/* This outputs a nice formatted version of the tree.                         */
template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out,
                          const Octree::Octree<T_DTYPE,T_DIM>& tree );

namespace Octree {

template<typename T_DTYPE, int T_DIM>
class Node {
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;
    typedef struct { VecD pos; T_DTYPE* object; } Datapoint;
    typedef std::list<Datapoint> Datalist;

/* Saving parent makes traversing easiers. Also a Node can get it's neighbor  *
 * through its parent, e.g. to check whether coarsening can be done           */
    class Node<T_DTYPE,T_DIM> * parent;

/* In order to avoid floating point errors center and size elements must be   *
 * in the form of 2 ^ ( ..., -1, 0, 1, ... ) ! Root Node is supposed have     *
 *     size = VecD(1), center = VecD(0)                                       */
    VecD center;
    VecD size;
    
private:
/* Pointer of these are NULL, if its a leaf node. The mapping between array   *
 * index and position is being done by ConvertNumberToDirection    -------    *
 * and vice versa. So 0 is the lower left back octant, while      | 1 | 3 |   *
 * 1 is the lower left front octant and so on. They are not        -------    *
 * stored in the order they should be traversed , therefore a     | 0 | 2 |   *
 *  a traversing mapping additionally                              -------    */
    const int nchildren = compileTime::pow(2,T_DIM);
    Node * children[ compileTime::pow(2,T_DIM) ];

/* Max data might not necessarily be uphold, if for example more than         *
 * OCTREE_MAX_OBJECTS_PER_LEAF are for some reason at the exact same point    *
 * the Octree would otherwise allocate infinite amount of smaller boxes. You  *
 * can use data.empty() to check whether this is a leaf or an interior node.  *
 * data holds only Pointers. User has to worry about allocating and freeing   *
 * himself                                                                    */
    const int maxdata = OCTREE_MAX_OBJECTS_PER_LEAF;
    Datalist data;
    
/* converts 001 to (-1,-1,1), 101 to (-1,1,-1) and so on. With this we can    *
 * map 0...7 to all eight neighbors in an Octree or respectively 0...3 for a  *
 * Quadtree or 0...2^N-1 for a tree in N dimensions                           */
    VecI ConvertNumberToDirection( const int  ) const;
    int  ConvertDirectionToNumber( const VecI ) const;
/* Traverses the Octree non-recursively deeper and deeper if necessary to     *
 * find the octant lead which encompasses a given position                    */
    Node * FindLeafContainingPos( const VecD & pos );
/* Give birth to eight new children with the correct coordinates and sizes.   *
 * It's private, so no reason to not give this a funny name :)                */
    void GrowUp(void);
/* Tries to collapse all children into the parent                             */
    bool Rejuvenate(void);
    bool HasOnlyLeaves(void) const;

public:
/* Forbid empty constructor, because we always have to create a node with a   *
 * center, size and parent !                                                  */
    Node(void);
/* We need no copy constructor, instead just move the pointers to these Nodes */
    Node(const Node &);
/* Free all allocated children if there are any                               */
    ~Node(void);
    Node& operator=(const Node &);
/* As all children pointer will be set to NULL a Node will be effectively     *
 * initialized as a Lead Node, but it can "grow up" and become a parent :)    *
 * if it has grown big enough, i.e. has exceeded OCTREE_MAX_OBJECTS_PER_LEAF  */
    Node(Node * const parent, VecD const cent, VecD const size);

/* Some Functions granting read access to private variables                   */
    const Node *     getChildPtr( const int i ) const;
    const Datapoint * getDataPtr( const int i ) const;
    
/* Test if a point lies inside the space described by this node               */
    bool IsInside ( const VecD & pos ) const;
    
/* This is for example needed to decide whether we have reached an end in our *
 * recursive search for some data                                             */
    bool IsLeaf( void ) const;
/* Insert data pointer at spatial position pos into the octtree. This doesn't *
 * check for double insertion, so beware. In order to avoid recursive calls   *
 * InsertData asserts that it is indeed a leaf. It should only be a           *
 * performance issue on PCs, because of pushing and popping all registers,    *
 * but I hate recursive calls. It crashed my Mindstorm NXT Robot              */
    void InsertData( const VecD pos, T_DTYPE * const object );
/* Remove Data from Octree. pos is needed for fast search for the object and  *
 * the data pointer is needed to identify the object with 100% accuracy for   *
 * the "rare" case two objects are at the same position. E.g. in a maxwell    *
 * solver to fit the initial conditions electrons and ions are initialized at *
 * the same location (!) in order to make initializing E more easy, because   *
 * it becomes zero everywhere. Returns false, if it couldn't find the datum.  */
    bool RemoveData( const VecD pos, T_DTYPE * const object );

/* <> would also suffice after operator<< instead of <T_DTYPE,T_DIM>, but one *
 * of both is needed or else the linker can't find the appropriate template   *
 * global function                                                            */
    friend std::ostream& operator<< <T_DTYPE,T_DIM>( std::ostream& out, const class Octree<T_DTYPE,T_DIM>& tree );
    friend class Octree<T_DTYPE,T_DIM>;
};


template<typename T_DTYPE, int T_DIM>
class Octree {
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;

public:
    VecD center, size;
    typedef class Node<T_DTYPE,T_DIM> Nodetype;
    typedef class Node<T_DTYPE,T_DIM> Node;
    Nodetype * root; // root can be leaf or child!

    Octree(void);
    Octree(const VecD center, const VecD size);
    Octree(const Octree &);
    ~Octree(void);
    Octree& operator=(const Octree &);

    VecD FindData( T_DTYPE * const object ) const;
    const Nodetype * GetNodePtr( const VecD center ) const;
    void InsertData( const VecD pos, T_DTYPE * const object );
    /* Returns false if datum not found */
    bool RemoveData( const VecD pos, T_DTYPE * const object );
    bool MoveData( const VecD pos, T_DTYPE * const object, const VecD newpos );

    friend std::ostream& operator<< <T_DTYPE,T_DIM>( std::ostream& out, const Octree<T_DTYPE,T_DIM>& tree );
};

} // namespace Octree


