#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "Octree.h"
#include "opengl/3D-Projection.h"
#include "math/Matrix.h"
#include "CompileTime.h"  // M_PI

#ifndef DEBUG_OCTREE_SVG
    #define DEBUG_OCTREE_SVG 0
#endif


namespace Octree {

template<int T_DIM>
class OctreeToSvg {
public:
    const int dim = T_DIM;
    typedef class Octree<T_DIM> OctreeType;
    typedef class Node<T_DIM> Node;

//private:
    std::ofstream out;

    OctreeType tree;
    const OctreeType * treesrc;

    Vec<int,2> imagesize;    // in px
    Vec<int,2> imageborder;  // in px
    CameraData Camera; // only needed if T_DIM == 3
    Matrix mvp;

    double currentTime = 0;
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
    OctreeToSvg( const OctreeType & tree, const std::string filename,
                 bool timestamp = true, int height = 600 );
    OctreeToSvg & operator=( const OctreeToSvg & src ) { assert(false); };
    OctreeToSvg( const OctreeToSvg & src ) { assert(false); };
    ~OctreeToSvg(void) {}
    void PrintGrid(void);
    /* colorfunc will be called with value in [0,1) and returns VecI with rgb *
     * delay >= 1./64. => 8 frames per second at max achievable with SVG      */
    template<typename T_FUNCTOR>
    void PrintTraversal( int pordering, T_FUNCTOR colorfunc, double delay = 1./64. );
    void PrintTraversal( int pordering, double delay = 1./64. );
    void AnimateUpdated( const OctreeType & newtree, double p_t = -1 );
    Vec<double,2> convertToImageCoordinates( const Vec<double,T_DIM> pos );
    void close(void);
};

} // namespace Octree

#include "OctreeToSvg.tpp"
