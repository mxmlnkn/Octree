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
#include "Vector.h"

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

template<typename T_DTYPE, int T_DIM> class Octree;

}

/* This outputs a nice formatted version of the tree.                         */
template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out,
                          Octree::Octree<T_DTYPE,T_DIM>& tree );


#include "OctreeNode.h"


namespace Octree {



template<typename T_DTYPE, int T_DIM>
class Octree {
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;

public:
    VecD center, size;
    typedef class Node<T_DTYPE,T_DIM> Node;
    Node * root; // root can be leaf or child!

    Octree(void);
    Octree(const VecD center, const VecD size);
    Octree(const Octree &);
    ~Octree(void);
    Octree& operator=(const Octree &);

    VecD FindData( T_DTYPE * const object ) const;
    const Node * GetNodePtr( const VecD center ) const;
    void InsertData( const VecD pos, T_DTYPE * const object );
    /* Returns false if datum not found */
    bool RemoveData( const VecD pos, T_DTYPE * const object );
    bool MoveData( const VecD pos, T_DTYPE * const object, const VecD newpos );

    /*bool CheckIntegrityParents( void );
    bool CheckIntegrityPositions( void );
    bool CheckIntegrityChilds( void ); */
    bool CheckIntegrity( void );

    class iterator {
    private:
        typedef struct{ int ichild; Node* node; } tododata;
        std::stack<tododata> todo;
    public:
        iterator(void);
        iterator( Node * );
        ~iterator(void);
        iterator(const iterator &);
        void operator= (const iterator &);

        iterator& operator++(void);
        bool operator==(const iterator &);
        bool operator!=(const iterator &);
        Node& operator*(void) const;
        Node* operator->(void) const;
    };
    iterator begin(void) const;
    iterator end(void) const;

    friend std::ostream& operator<< <T_DTYPE,T_DIM>( std::ostream& out, Octree<T_DTYPE,T_DIM>& tree );
};



} // namespace Octree



#include "OctreeNode.tpp"
#include "OctreeIterator.tpp"
#include "Octree.tpp"
