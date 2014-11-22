/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe
 
*/

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc

namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime

#include "Vector.h"
#include "Octree.h"
#include "TeeStream.h"

#include "Vector.tpp"
#include "Octree.tpp"

#define SIMDIM 2
typedef Vec<double,SIMDIM> VecD;


int main( int argc, char **argv )
{
    tout.Open("");

    const int NValues = 20;
    VecD positions[NValues];
    int  data     [NValues] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    positions[ 0][0] =  0.680369; positions[ 0][1] = -0.132787;
    positions[ 1][0] = -0.585351; positions[ 1][1] = -0.054848;
    positions[ 2][0] = -0.158091; positions[ 2][1] = -0.867561;
    positions[ 3][0] = -0.477332; positions[ 3][1] =  0.798686;
    positions[ 4][0] =  0.752312; positions[ 4][1] =  0.770033;
    positions[ 5][0] =  0.211714; positions[ 5][1] = -0.661169;
    positions[ 6][0] =  0.935703; positions[ 6][1] = -0.695199;
    positions[ 7][0] =  0.597030; positions[ 7][1] =  0.767231;
    positions[ 8][0] = -0.853328; positions[ 8][1] =  0.706876;
    positions[ 9][0] =  0.178442; positions[ 9][1] =  0.947039;
    positions[10][0] =  0.054220; positions[10][1] =  0.745018;
    positions[11][0] = -0.009938; positions[11][1] = -0.720338;
    positions[12][0] =  0.494300; positions[12][1] = -0.701520;
    positions[13][0] = -0.723740; positions[13][1] =  0.268552;
    positions[14][0] =  0.709980; positions[14][1] =  0.930955;
    positions[15][0] =  0.487417; positions[15][1] =  0.908338;
    positions[16][0] =  0.124996; positions[16][1] = -0.501056;
    positions[17][0] = -0.959923; positions[17][1] = -0.662497;
    positions[18][0] =  0.015413; positions[18][1] =  0.891808;
    positions[19][0] = -0.589081; positions[19][1] =  0.753129;
    VecD size(2), center(0);
    Octree::Octree<int,SIMDIM> tree( center, size );
    for (int i=0; i<20; i++) {
        tree.InsertData( positions[i], &(data[i]) );
    }
    std::cout << tree;
}