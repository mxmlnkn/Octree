#pragma once

#include <iostream>
#include <fstream>
#include <map>

#include "Octree.h"

namespace Octree {

template<typename T_DTYPE, int T_DIM>
class OctreeToSvg {
public:
    typedef Vec<double,T_DIM> VecD;
    typedef Vec<int   ,T_DIM> VecI;
    const int dim = T_DIM;
    typedef class Octree<T_DTYPE,T_DIM> Octreetype;
    typedef class Node<T_DTYPE,T_DIM> Node;
    
//private:
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
    
    /* Thanks to https://stackoverflow.com/questions/16362231/ */
    struct StrictWeakOrderingVecD {
        bool operator()( const VecD & a, const VecD & b ) const {
            for (int i=0; i<T_DIM; i++) {
                if (a[i] < b[i]) return true;
                if (a[i] > b[i]) return false;
            }
            return false;
    } };
/* Boxes Drawn is being called with the centers of each square, as those are  *
 * unique. For each unique center a unique ID and a boolean, whether it is    *
 * currently visible (the last <set/> was "stroke:white" or not. If a box can *
 * be found in this map, than it must already have been drawn with <rect/>    */
    typedef struct { int id; bool visible; } Keyvalues;
    typedef std::map<VecD,Keyvalues,StrictWeakOrderingVecD> VecDMap;
    int NboxesDrawn = 0;
    VecDMap boxesDrawn;
    typename VecDMap::iterator boxesDrawnIt;

public:
    OctreeToSvg(void) { assert(false); }
    OctreeToSvg(const Octreetype & tree, const std::string filename);
    ~OctreeToSvg(void) {
        out << "</svg>\n";
        out.close();
    }
    void PrintGrid(void);
    void PrintPositions(void);
    void AnimateUpdated( const Octreetype & newtree );
    VecD convertToImageCoordinates( const VecD pos );
};


} // namespace Octree