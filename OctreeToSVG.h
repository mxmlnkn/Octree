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
    std::map<size_t,int> BoxesDrawn;
    std::map<char,int>::iterator BoxesDrawnIt;

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
