#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <ctime>
#include <sstream>
#include <string>

#include "Octree.h"
#include "opengl/3D-Projection.h"
#include "math/Matrix.h"

namespace Octree {

template<int T_DIM>
class OctreeToSvg {
public:
    const int dim = T_DIM;
    typedef class Octree<T_DIM> Octreetype;
    typedef class Node<T_DIM> Node;

//private:
    std::ofstream out;

    Octreetype tree;
    const Octreetype * treesrc;

    Vec<int,2> imagesize;    // in px
    Vec<int,2> imageborder;  // in px
    CameraData Camera; // only needed if T_DIM == 3
    Matrix mvp;

    int currentTime = 1;
    const int DUR = 2; //2s per update
    const int STROKE_WIDTH = 1;

    /* Thanks to https://stackoverflow.com/questions/16362231/ */
    struct StrictWeakOrderingVecD {
        bool operator()( const Vec<double,T_DIM> & a, const Vec<double,T_DIM> & b ) const {
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
    typedef std::map<Vec<double,T_DIM>,Keyvalues,StrictWeakOrderingVecD> VecDMap;
    int NboxesDrawn = 0;
    VecDMap boxesDrawn;

public:
    OctreeToSvg(void) { assert(false); }
    OctreeToSvg( const Octreetype & tree, const std::string filename,
                 bool timestamp = true, int height = 600 );
    OctreeToSvg & operator=( const OctreeToSvg & src ) { assert(false); };
    OctreeToSvg( const OctreeToSvg & src ) { assert(false); };
    ~OctreeToSvg(void) {}
    void PrintGrid(void);
    void PrintPositions(void);
    void AnimateUpdated( const Octreetype & newtree );
    Vec<double,2> convertToImageCoordinates( const Vec<double,T_DIM> pos );
    void close(void) {
        out << "</svg>\n";
        out.close();
    }
};

} // namespace Octree

#include "OctreeToSvg.tpp"
