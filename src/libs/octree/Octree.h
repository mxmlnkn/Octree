/*
ToDo:
  - GetSize Function / Get Max Levels
  - max depth
  - Add CheckIntegrity Method, which checks for:
      . parents
      . child-pointers
      . are all points really inside their boxes ?
  = Bigger goal: parallelize
*/

#pragma once

#include <cassert>
#include <iostream>
#include <cstdlib>   // malloc
#include <list>
#include <stack>
#include <ctime> // only for debug in operator<<
#include "math/TVector.h"

#ifndef OCTREE_MAX_OBJECTS_PER_LEAF
	#define OCTREE_MAX_OBJECTS_PER_LEAF 3
#endif
/* Sets max level depth to 10, meaning at maximum 8^10 = 2^30 leaf nodes can  *
 * be created. Much more doesn't make sense anyway, because of memory and     *
 * because dx i.e. size becomes so small, that the centers of neighbors       *
 * decoded in  doubles won't differe from each other, meaning: center+size    *
 * =center. This happens at least at the sides (center ~ 1) (not in the       *
 * center at 0) for size == DBL_EPSILON = 2.22e-16. Size is halved at every   *
 * level, so size(level) = 1/2^level == 2.22e-16 => level = 52                */
#define OCTREE_MAX_DEPTH 10


namespace Octree {

template<int T_DIM> class Octree;

}

/* This outputs a nice formatted version of the tree.                         */
template<int T_DIM>
std::ostream& operator<<( std::ostream& out,
                          Octree::Octree<T_DIM>& tree );

#include "OctreeNode.h"

namespace Octree {

namespace Ordering {
    const int Morton   = 0;
    const int GrayCode = 1;
    const int Hilbert  = 2;
}

template< int T_DIM>
class Octree {
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;

public:
    VecD center, size;
    typedef class Node<T_DIM> Node;
    typedef typename Node::iterator iterator;
    Node * root; // root can be leaf or child!

    Octree(void);
    /* This is a feature of Octree. It can convert the fixed internal sizes   *
     * and positions of OctreeNode to arbitrary ones. Internal sizes are      *
     * fixed to avoid floating point rounding errors                          */
    Octree(const VecD center, const VecD size);
    Octree(const Octree &);
    ~Octree(void);
    Octree& operator=(const Octree &);
    
    VecD FindData( void * const object ) const;
    /* Returns Pointer of Node with its center equaling exactly the given one */
    const Node * GetNodePtr( const VecD center ) const;
    void InsertData( const VecD pos, void * const object );
    /* Returns false if datum not found */
    bool RemoveData( const VecD pos, const int index );
    Node * FindLeafContainingPos( const VecD & pos );
    
    /* Returns level at which first leaf can be found */
    int getMinLevel( void );
    /* Returns maximum depth of tree */
    int getMaxLevel( void );

    /* Returns true, if it found not structural error. It checks for:         *
     *  - Datapoint (particle) being inside the Octree-cell its bein stored in*
     *  - all children pointing to the parent in which they are being stored  *
     *  - leaf nodes containing illegitimate children                         *
     *  - pointer to childs being NULL / missing                              */
    bool CheckIntegrity( void ) const;

    /* returns iterator with only root-node in todo stack */
    typename Node::iterator begin(int ordering = Ordering::Morton) const;
    /* returns iterator with empty stack */
    typename Node::iterator end(void) const;

    friend std::ostream& operator<< <T_DIM>( std::ostream& out, Octree<T_DIM>& tree );
    /* internal coords have a range of [-0.5,0.5] in every dimension */
    VecD toInternalCoords( const VecD pos ) const;
    VecD toGlobalCoords( const VecD pos ) const;
};



} // namespace Octree



#include "OctreeNode.tpp"
#include "OctreeIterator.tpp"
#include "Octree.tpp"
