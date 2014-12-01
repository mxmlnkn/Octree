/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe

*/

#ifndef M_PI
#   define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc, srand, rand, RAND_MAX
#include <random>   // normal_distribution

namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime


#define SIMDIM 2

#include "Vector.h"
#include "BaseMatrix.h"
#define YEE_CELL_TIMESTEPS_TO_SAVE 3
#include "YeeCell.h"
typedef BaseMatrix<class YeeCell,SIMDIM> CellMatrix;
#include "TeeStream.h"

#include "Vector.tpp"
typedef Vec<double,SIMDIM> VecD;
typedef Vec<int   ,SIMDIM> VecI;

#include <pngwriter.h>

#define N_CELLS_X 256 // should be power of two, because of octree !
#define N_CELLS_Y 128 // should be same as above, because octree isn't able to have a different amount of equal sized cells in two directions !!!
#define CELL_SIZE_X 1
#define CELL_SIZE_Y 1

double t_spawn_func( int t ) {
	int Nlambda  = 20;
	double sigma = Nlambda/2.;
	double t0    = Nlambda;
	return sin(1.*t/Nlambda);
	return exp(-pow(t-t0,2)/(2.*sigma*sigma));
	return 1./(sigma*sqrt(2.*M_PI))*exp(-pow(t-t0,2)/(2.*sigma*sigma));
	return ( t < Nlambda ? 1 : 0 );
}

int main( int argc, char **argv )
{
    tout.Open("out");

    /* create data buffer with cells which will be inserted into the octree   *
	 * and fill it with some initial data (zero or random)                    */
	VecI size(0); size[0]=N_CELLS_X; size[1]=N_CELLS_Y;
	CellMatrix data(size);

	for ( int timestep=0; timestep < 400; ++timestep ) {
        int tprev = (YEE_CELL_TIMESTEPS_TO_SAVE + timestep-1) % YEE_CELL_TIMESTEPS_TO_SAVE;
        int tcur  = (timestep  ) % YEE_CELL_TIMESTEPS_TO_SAVE;
        int tnext = (timestep+1) % YEE_CELL_TIMESTEPS_TO_SAVE;
		
		VecI pos(0); pos[0] = 0;
		//std::cout << tprev << "," << tcur << "," << tnext << std::endl;
		data[pos].E[tprev][0] = t_spawn_func( timestep-1 );
        data[pos].E[tcur ][0] = t_spawn_func( timestep   );
        data[pos].E[tnext][0] = t_spawn_func( timestep+1 );

		for ( int curCell = 0; curCell < N_CELLS_X; curCell++ ) {
			VecI xprev(0); xprev[0] = curCell-1;
			VecI xcur (0); xcur [0] = curCell  ;
			VecI xnext(0); xnext[0] = curCell+1;
			//std::cout << curCell << "," << xprev << "," << xcur << "," << xnext << std::endl;
			double utcxp = xprev[0] >= 0 ? data[xprev].E[tcur][0] : 0;
			double utcxn = xnext[0] < N_CELLS_X ? data[xnext].E[tcur][0] : 0;
			
			double S = 1;
			data[xcur].E[tnext][0] = S*S * ( utcxn - 2*data[xcur].E[tcur][0] +  utcxp ) + 2*data[xcur].E[tcur][0] - data[xcur].E[tprev][0];
		}
		
		/*for ( int ix=0; ix < N_CELLS_X; ++ix ) {
			VecI pos; pos[0]=ix; pos[1]=0;
			std::cout << data[pos].E[tnext][0];
		}
		std::cout << std::endl;*/
		
		if (timestep % 2 == 0) {
			static int framecounter = 0;
			char filename[100];
			sprintf( filename, "output/E_%05i.png", framecounter++ );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=0;
				image.plot( ix,iy, data[pos].E[tnext][0], -data[pos].E[tnext][0], 0.0);
			}
			image.close();
		}
	}

	/* Create png with as much pixels as cells with initially black background*/
	/*pngwriter image( N_CELLS_X,N_CELLS_Y,0,"E.png" );
    for ( int i=0; i < data.getSize().product(); ++i ) {
        VecI pos = data.getVectorIndex( i );
		image.plot( pos[0], pos[1], data[i].E(0)[0], data[i].E(0)[1], data[i].E(0)[2]);
	}
	image.close();*/
}