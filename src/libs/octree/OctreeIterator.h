#pragma once

#include <ctime>
#include "OctreeNode.h"

namespace Octree{

template<int T_DIM>
class OctreeIterator {
friend class Node<T_DIM>;
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;
    typedef Node<T_DIM> Node;
    typedef struct{ int ichild; Node* node; int orientation; } ToDoData;

private:
    Node* root;
    std::stack<ToDoData> todo;
    std::list<Node*> done; /* for ordering method 3 (Rows) */
    int curpos;
/* Ordering 0: Morton (Depth first traversal), 1: GrayCode, 2: Hilbert    */
    int ordering;
    std::vector<Node*> traversallist;
/* [ordering method,see above][parentOrientiation][n-th child to traverse]*
 * => child index like used to access children-array of a Node and        *
 *  orientation for that child                                            *
 * ToDo: Allocate and Calculate this dynamically on very first OctreeIterator   *
 *       access. Also this only works for 2D !!!                          */
    static struct orderingTableStruct {
        struct {
            int ordering   [ compileTime::pow(2,T_DIM-1) ][ compileTime::pow(2,T_DIM) ];
            int orientation[ compileTime::pow(2,T_DIM-1) ][ compileTime::pow(2,T_DIM) ];
        } GrayCode;
        struct {
            int ordering   [24][ compileTime::pow(2,T_DIM) ];
            int orientation[24][ compileTime::pow(2,T_DIM) ];
        } Hilbert;
    } orderingTable;

public:
    inline OctreeIterator(void);
    OctreeIterator( Node *, int ordering = 0 );
    inline ~OctreeIterator(void);
    inline OctreeIterator(const OctreeIterator &);
    inline OctreeIterator & operator= (const OctreeIterator &);

    OctreeIterator & operator++(void);
    inline OctreeIterator operator++( int unused );
    inline bool operator==(const OctreeIterator &);
    inline bool operator!=(const OctreeIterator &);
    inline Node& operator*(void) const;
    inline Node* operator->(void) const;
    inline OctreeIterator begin(void);
    /* return iterator with empty stack among other things (void constructor) */
    inline OctreeIterator end(void) const;
};

} // namespace Octree
