#pragma once

#include <cassert>
#include <iostream>
#include <cstdlib>   // malloc
#include <vector>
#include <stack>
#include <list>
#include <algorithm> // find
#include "math/TVector.h"
#include "math/TBaseMatrix.h" // convertLinearToVectorIndex
#include "Directions.h" // getOppositeDirection
#include "teestream/TeeStream.h" // tout
#include "CompileTime.h"


namespace Octree {

template<int T_DIM>
class Node {
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;
    typedef std::vector<void*> Datalist;

/* Saving parent makes traversing easiers. Also a Node can get it's neighbor  *
 * through its parent, e.g. to check whether coarsening can be done           */
    class Node<T_DIM> * parent;

/* In order to avoid floating point errors center and size elements must be   *
 * in the form of 2 ^ ( ..., -1, 0, 1, ... ) ! Root Node is supposed have     *
 *     size = VecD(1), center = VecD(0)                                       */
    VecD center;
    VecD size;

/* Max data might not necessarily be uphold, if for example more than         *
 * OCTREE_MAX_OBJECTS_PER_LEAF are for some reason at the exact same point    *
 * the Octree would otherwise allocate infinite amount of smaller boxes. You  *
 * can use data.empty() to check whether this is a leaf or an interior node.  *
 * data holds only Pointers. User has to worry about allocating and freeing   *
 * himself                                                                    */
    const int maxdata = OCTREE_MAX_OBJECTS_PER_LEAF;
    Datalist data;

private:
/* Pointer of these are NULL, if its a leaf node. The mapping between array   *
 * index and position is being done by ConvertNumberToDirection    -------    *
 * and vice versa. So 0 is the lower left back octant, while      | 1 | 0 |   *
 * 1 is the lower left front octant and so on. They are not        -------    *
 * stored in the order they should be traversed , therefore a     | 3 | 2 |   *
 *  a traversing mapping additionally                              -------    */
    const int nchildren;
    Node * children[ compileTime::pow(2,T_DIM) ];


/* converts 001 to (+1,+1,-1), 101 to (-1,-1,+1) and so on. The binary        *
 * representation is essentially treated as an array of sign bits. With this  *
 * we can map 0...7 to all eight neighbors in an Octree or respectively 0...3 *
 * for a Quadtree or 0...2^N-1 for a tree in N dimensions                     */
    VecI ConvertNumberToDirection( const int number ) const;
/* If the direction contains an entry which is not 1 or -1 then it's not the  *
 * format we expect => return -1                                              */
    static int ConvertDirectionToNumber( const VecI );
/* Traverses the Octree non-recursively deeper and deeper if necessary to     *
 * find the octant lead which encompasses a given position                    */
    Node * FindLeafContainingPos( const VecD & pos );

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


    class iterator {
    private:
        typedef struct{ int ichild; Node* node; int orientation; } tododata;
        std::stack<tododata> todo;
        std::list<Node*> done; /* for ordering method 3 (Rows) */
        int curpos;
    /* Ordering 0: Morton (Depth first traversal), 1: GrayCode, 2: Hilbert    */
        int ordering;
    /* [ordering method,see above][parentOrientiation][n-th child to traverse]*
     * => child index like used to access children-array of a Node and        *
     *  orientation for that child                                            *
     * ToDo: Allocate and Calculate this dynamically on very first iterator   *
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
        iterator(void);
        iterator( Node *, int ordering = 0 );
        ~iterator(void);
        iterator(const iterator &);
        iterator & operator= (const iterator &);

        iterator& operator++(void);
        iterator operator++( int unused );
        bool operator==(const iterator &);
        bool operator!=(const iterator &);
        Node& operator*(void) const;
        Node* operator->(void) const;
        iterator end(void);
    };
/* returns iterator with only root-node in todo stack */
    iterator begin(int ordering = 0);
/* returns iterator with empty stack */
    iterator end(void);


/* Some Functions granting read access to private variables                   */
    const Node * getChildPtr( const int i ) const;
/* Because this Octree initially was intended to hold Particles with          *
 * attributes, it is possible to store MAX_OBJECTS_PER_LEAF void pointers     *
 * The parameter i specifies which one of those pointers we want.             *
 * This may become deprecated, as the Octree should be kept to its basic      *
 * functionality. The user can store a pointer to an array himself, if he     *
 * needs that                                                                 */
    void* getDataPtr( const int i ) const;
    VecD getSize( void ) const;

/* Test if a point lies inside the space described by this node               */
    bool IsInside ( const VecD & pos ) const;

/* Give birth to eight new children with the correct coordinates and sizes.   *
 * It's private, so no reason to not give this a funny name :)                */
    void GrowUp(void);
/* Tries to collapse all children into the parent                             */
    bool Rejuvenate(void);
    bool HasOnlyLeaves(void) const;
    void DeleteChildren(void);

/* This is for example needed to decide whether we have reached an end in our *
 * recursive search for some data                                             */
    bool IsLeaf( void ) const;

/* Returns the depth for this node. The root node returns 0. If the root leaf *
 * becomes a parent, then the children will return 1. This being calculated   *
 * with the parents member                                                    */
    int getLevel( void ) const;
/* Returns level at which first leaf can be found */
    int getMinLevel( void );
/* Returns maximum depth of tree */
    int getMaxLevel( void );
    
/* Returns total amount of leaves, including those contained by childnodes if *
 * existing                                                                   */
    int countLeaves( void );
/* returns neighboring cell in direction. Works with diagonal directions and  *
 * also with VecI(0), returning the leaf itself then.                         *
 * If the neighbor cell in that direction is larger and has no                *
 * childnodes, then that is returned. If this cell is larger than the neigh-  *
 * boring cells in that direction, then the cell of same size will be         *
 * returned. The format of direction is a movement vector, which will just be *
 * added to the position of this node. BEWARE!!! Returns NULL if not periodic *
 * and on border, meaning there is no neighbor                                */
    Node * getNeighbor( const VecI direction, const VecI periodic );
/* returns all neighbors of that cell as an iterator with a prefilled todo-   *
 * stack of leaf nodes. Because they are leaf nodes the iterator++ will just  *
 * pop one after another until it's empty. This just expands the neighbor     *
 * if it is not a leaf node, because it has smaller children than the source  */
    iterator getNeighbors( const VecI direction, const VecI periodic );


/* <> would also suffice after operator<< instead of <T_DIM>, but one *
 * of both is needed or else the linker can't find the appropriate template   *
 * global function                                                            */
    friend std::ostream& operator<< <T_DIM>( std::ostream& out, class Octree<T_DIM>& tree );
    friend class Octree<T_DIM>;
};



} // namespace Octree
