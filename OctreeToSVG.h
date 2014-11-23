#pragma once

#include <iostream>
#include <fstream>
#include <map>

#include "Octree.h"

namespace Octree {

template<typename T_DTYPE, int T_DIM>
class OctreeToSVG {
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;
    typedef class Octree<T_DTYPE,T_DIM> Octreetype;
    typedef class Node<T_DTYPE,T_DIM> Node;
    
private:
    std::ofstream out;
    Octreetype tree;
    const double borderx;
    const double bordery;
    const double height;
    const double width;
    VecD imagesize;
    VecD imageborder;
    
    int currentTime = 1;
    const int DUR = 2; //2s per update
    const int STROKE_WIDTH = 1;
    
    class OrCompareVecD {
        public:
        bool operator()( const VecD & a, const VecD & b ) const {
            bool result = false;
            for (int i=0; i<T_DIM; i++)
                result = result or (a[i] < b[i]);
            return result;
    } };
    
/* Boxes Drawn is being called with the centers of each square, as those are  *
 * unique. For each unique center a unique ID and a boolean, whether it is    *
 * currently visible (the last <set/> was "stroke:white" or not. If a box can *
 * be found in this map, than it must already have been drawn with <rect/>    */
    typedef struct { int id; bool visible; } Keyvalues;
    int NBoxesDrawn = 0;
    std::map<VecD,Keyvalues,OrCompareVecD> BoxesDrawn;
    typename std::map<VecD,Keyvalues,OrCompareVecD>::iterator BoxesDrawnIt;

public:
    OctreeToSVG(void) { assert(false); }
    OctreeToSVG(const Octreetype & tree, const std::string filename);
    ~OctreeToSVG(void) {
        out << "</svg>\n";
        out.close();
    }
    void PrintGrid(void);
    void PrintPositions(void);
    void AnimateUpdated( const Octreetype & newtree );
};


} // namespace Octree
