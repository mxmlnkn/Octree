/*

rm testOctree.exe; g++ testOctree.cpp -o testOctree.exe -Wall -std=c++0x; ./testOctree.exe

ToDo:
  - check courant criterium !!! (automatically)
  - Check for energy conservation automatically (problematic for initial wave)
  - draw double slit into png in gray and shift wave from red-green to only red
  - calculate and draw into the png the intensity after the obstacle
  - Introduce SI and internal Units!
  - above is for CELL_SIZE_X = 1.0
  - for CELL_SIZE_X = 5.0 the wave doesn't move to the right :S?
  - why no tprev needed like in 1D wave equation ... 
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


#define SIMDIM 3

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

#define N_CELLS_X 128 // should be power of two, because of octree !
#define N_CELLS_Y 128 // should be same as above, because octree isn't able to have a different amount of equal sized cells in two directions !!!
#define N_CELLS_Z 1
#define CELL_SIZE_X 2.0
#define CELL_SIZE_Y 2.0
#define CELL_SIZE_Z 2.0
#define Dt 1.0

double t_spawn_func( double t ) {
	double Nlambda  = 40;
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

    /* create data buffer with cells initially with zeros in all entries      */
	VecI size(0); size[0]=N_CELLS_X; size[1]=N_CELLS_Y; size[2]=N_CELLS_Z;
	CellMatrix data(size);
	for ( int x = N_CELLS_X/2; x < N_CELLS_X/2+10; x++ )
	for ( int y = 0; y < N_CELLS_Y; y++ ) {
		VecI pos(0); pos[0]=x; pos[1]=y;
		data[pos].rhoprime = 1e40;
		data[pos].sigma    = 1e40;
	}


	for ( int timestep=0; timestep < 400; ++timestep )
	{
			const int X = 0;
			const int Y = 1;
			const int Z = 2;

        //int tprev = (YEE_CELL_TIMESTEPS_TO_SAVE + timestep-1) % YEE_CELL_TIMESTEPS_TO_SAVE;
        int tcur  = (timestep  ) % YEE_CELL_TIMESTEPS_TO_SAVE;
        int tnext = (timestep+1) % YEE_CELL_TIMESTEPS_TO_SAVE;

		VecI pos(0); pos[X] = 0;
		for ( int y = 0; y < N_CELLS_Y; y++ ) {
			pos[Y] = y;
			//data[pos].E[tprev][Z] = t_spawn_func( timestep-1 );
			data[pos].E[tcur ][Z] = t_spawn_func( timestep     );
			data[pos].H[tcur ][Y] = t_spawn_func( timestep-0.5 );
			//data[pos].E[tnext][Z] = t_spawn_func( timestep+1 );
		}

		for ( int x = 0; x < N_CELLS_X; x++ )
		for ( int y = 0; y < N_CELLS_Y; y++ )
		for ( int z = 0; z < 1; z++ ) {
			VecI pos; pos[X]=x; pos[Y]=y; pos[Z]=z;

			VecI xprev=pos; xprev[X]--;
			VecI xnext=pos; xnext[X]++;

			VecI yprev=pos; yprev[Y]--;
			VecI ynext=pos; ynext[Y]++;

			VecI zprev=pos; zprev[Z]--;
			VecI znext=pos; znext[Z]++;
			
			/* Make it periodic in y-direction */
			if (yprev[Y] <  0        ) yprev[Y] += N_CELLS_Y;
			if (ynext[Y] >= N_CELLS_Y) ynext[Y] -= N_CELLS_Y;

			/* Update all H components */
			{
			double Eztcyp = y-1 >= 0        ? data[yprev].E[tcur][Z] : 0;
			double Eztcyn = y+1 < N_CELLS_Y ? data[ynext].E[tcur][Z] : 0;
			//double Eytczp = z-1 >= 0        ? data[zprev].E[tcur][Y] : 0;
			//double Eytczn = z+1 < N_CELLS_Z ? data[znext].E[tcur][Y] : 0;
			data[pos].H[tnext][X] = -2.*Dt / (2*data[pos].mu + data[pos].rhoprime*Dt ) * ( ( Eztcyn - Eztcyp ) / CELL_SIZE_Y /* - ( Eytczn - Eytczp ) / CELL_SIZE_Z */ - 0.5*data[pos].rhoprime * data[pos].H[tcur][X]);
			}{
			//double Extczp = z-1 >= 0        ? data[zprev].E[tcur][X] : 0;
			//double Extczn = z+1 < N_CELLS_Z ? data[znext].E[tcur][X] : 0;
			double Eztcxp = x-1 >= 0        ? data[xprev].E[tcur][Z] : 0;
			double Eztcxn = x+1 < N_CELLS_X ? data[xnext].E[tcur][Z] : 0;
			data[pos].H[tnext][Y] = -2.*Dt / (2*data[pos].mu + data[pos].rhoprime*Dt ) * ( /*( Extczn - Extczp ) / CELL_SIZE_Z - */ ( Eztcxn - Eztcxp ) / CELL_SIZE_X  - 0.5*data[pos].rhoprime * data[pos].H[tcur][Y] );
			}{
			double Eytcxp = x-1 >= 0        ? data[xprev].E[tcur][Y] : 0;
			double Eytcxn = x+1 < N_CELLS_X ? data[xnext].E[tcur][Y] : 0;
			double Extcyp = y-1 >= 0        ? data[yprev].E[tcur][X] : 0;
			double Extcyn = y+1 < N_CELLS_Y ? data[ynext].E[tcur][X] : 0;
			data[pos].H[tnext][Z] = -2.*Dt / (2*data[pos].mu + data[pos].rhoprime*Dt ) * ( ( Eytcxn - Eytcxp ) / CELL_SIZE_X - ( Extcyn - Extcyp ) / CELL_SIZE_Y - 0.5*data[pos].rhoprime * data[pos].H[tcur][Z] );
			}
			/* Now update all E components */
			{
			double Hztcyp = y-1 >= 0        ? data[yprev].H[tcur][Z] : 0;
			double Hztcyn = y+1 < N_CELLS_Y ? data[ynext].H[tcur][Z] : 0;
			//double Hytczp = z-1 >= 0        ? data[zprev].H[tcur][Y] : 0;
			//double Hytczn = z+1 < N_CELLS_Z ? data[znext].H[tcur][Y] : 0;
			data[pos].E[tnext][X] = 2.*Dt / (2*data[pos].epsilon + data[pos].sigma*Dt ) * ( ( Hztcyn - Hztcyp ) / CELL_SIZE_Y /* - ( Hytczn - Hytczp ) / CELL_SIZE_Z */ - 0.5*data[pos].sigma * data[pos].E[tcur] [X]);
			}{
			//double Hxtczp = z-1 >= 0        ? data[zprev].H[tcur][X] : 0;
			//double Hxtczn = z+1 < N_CELLS_Z ? data[znext].H[tcur][X] : 0;
			double Hztcxp = x-1 >= 0        ? data[xprev].H[tcur][Z] : 0;
			double Hztcxn = x+1 < N_CELLS_X ? data[xnext].H[tcur][Z] : 0;
			data[pos].E[tnext][Y] = 2.*Dt / (2*data[pos].epsilon + data[pos].sigma*Dt ) * ( /*( Hxtczn - Hxtczp ) / CELL_SIZE_Z - */ ( Hztcxn - Hztcxp ) / CELL_SIZE_X  - 0.5*data[pos].sigma * data[pos].E[tcur][Y] );
			}{
			double Hytcxp = x-1 >= 0        ? data[xprev].H[tcur][Y] : 0;
			double Hytcxn = x+1 < N_CELLS_X ? data[xnext].H[tcur][Y] : 0;
			double Hxtcyp = y-1 >= 0        ? data[yprev].H[tcur][X] : 0;
			double Hxtcyn = y+1 < N_CELLS_Y ? data[ynext].H[tcur][X] : 0;
			data[pos].E[tnext][Z] = 2.*Dt / (2*data[pos].epsilon + data[pos].sigma*Dt ) * ( ( Hytcxn - Hytcxp ) / CELL_SIZE_X - ( Hxtcyn - Hxtcyp ) / CELL_SIZE_Y - 0.5*data[pos].sigma * data[pos].E[tcur][Z] );
			}
		}

		/*for ( int ix=0; ix < N_CELLS_X; ++ix ) {
		    VecI pos; pos[0]=ix; pos[1]=0;
		    std::cout << data[pos].E[tnext][Z];
		}
		std::cout << std::endl;*/
		
		if (timestep % 2 == 0) {
			static int framecounter = 0;
			framecounter++;
			char filename[100];
			{
			sprintf( filename, "output/Ex_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix,iy, data[pos].E[tnext][X], -data[pos].E[tnext][X], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Ey_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix,iy, data[pos].E[tnext][Y], -data[pos].E[tnext][Y], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Ez_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix,iy, data[pos].E[tnext][Z], -data[pos].E[tnext][Z], 0.0);
			}
			image.close();
			}

			{
			sprintf( filename, "output/Hx_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix,iy, data[pos].H[tnext][X], -data[pos].H[tnext][X], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Hy_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix,iy, data[pos].H[tnext][Y], -data[pos].H[tnext][Y], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Hz_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix,iy, data[pos].H[tnext][Z], -data[pos].H[tnext][Z], 0.0);
			}
			image.close();
			}

			std::cout << "Image " << framecounter << std::endl;
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