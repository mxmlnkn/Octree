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
#define INF (1.0/0.0)


#include <iostream>
#include <cmath>    // sin
#include <cstdlib>  // malloc, srand, rand, RAND_MAX
#include <random>   // normal_distribution

#include <pngwriter.h>


namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime


/* General libraries */
#include "Vector.h"
#include "Vector.tpp"
#include "BaseMatrix.h"
#include "TeeStream.h"

/* Simulation Code Includes */
#include "Parameters.cpp"
typedef Vec<double,SIMDIM> VecD;
typedef Vec<int   ,SIMDIM> VecI;
#define YEE_CELL_TIMESTEPS_TO_SAVE 100
#include "YeeCell.h"
typedef BaseMatrix<class YeeCell,SIMDIM> CellMatrix;


#define N_CELLS_X NUMBER_OF_CELLS_X // should be power of two, because of octree !
#define N_CELLS_Y NUMBER_OF_CELLS_Y // should be same as above, because octree isn't able to have a different amount of equal sized cells in two directions !!!
#define N_CELLS_Z NUMBER_OF_CELLS_Z

double t_spawn_func( double t ) {
	double Nlambda  = 40;
	double sigma = Nlambda/2.;
	double t0    = Nlambda;
	return sin(2.*M_PI*t/Nlambda);
	return ( t < Nlambda ? 1.0 : 0.0 );
	return exp(-pow(t-t0,2)/(2.*sigma*sigma));
	return 1./(sigma*sqrt(2.*M_PI))*exp(-pow(t-t0,2)/(2.*sigma*sigma));
}

int main( int argc, char **argv )
{
    tout.Open("out");

    tout << "MUE0                 : " << MUE0                       << "\n";
    tout << "EPS0                 : " << EPS0                       << "\n";
    tout << "SPEED_OF_LIGHT       : " << SPEED_OF_LIGHT             << "\n";
    tout << "CELL_SIZE_SI         : " << CELL_SIZE_SI               << "\n";
    tout << "CELL_SIZE            : " << CELL_SIZE_SI / UNIT_LENGTH << "\n";
    tout << "PARTICLE_RADIUS      : " << PARTICLE_RADIUS            << "\n";
    tout << "NUMBER_OF_CELLS_X    : " << NUMBER_OF_CELLS_X          << "\n";
    tout << "NUMBER_OF_CELLS_Y    : " << NUMBER_OF_CELLS_Y          << "\n";
    tout << "NUMBER_OF_CELLS_Z    : " << NUMBER_OF_CELLS_Z          << "\n";
    tout << "DELTA_T_SI           : " << DELTA_T_SI                 << "\n";
    tout << "UNIT_ENERGY          : " << UNIT_ENERGY                << "\n";
    tout << "UNIT_MOMENTUM        : " << UNIT_MOMENTUM              << "\n";
    tout << "UNIT_ANGULAR_MOMENTUM: " << UNIT_ANGULAR_MOMENTUM      << "\n";
    tout << "ELECTRON_MASS        : " << ELECTRON_MASS              << "\n";
    tout << "\n";

    /* create data buffer with cells initially with zeros in all entries      */
	VecI size(0);
	for ( int i = 0; i < SIMDIM; i++ )
		size[i] = NUMBER_OF_CELLS[i];
	CellMatrix data(size);
	
	for ( int i = 0; i < data.getSize().product(); i++ ) {
		data[i].epsilon    = EPS0;
		data[i].mu         = MUE0;
		data[i].rhoprime = 0;
		data[i].sigma    = 0;
	}
	/*for ( int x = N_CELLS_X/2; x < N_CELLS_X/2+10; x++ )
	for ( int y = 0; y < N_CELLS_Y; y++ ) {
		VecI pos(0); pos[0]=x; pos[1]=y;
		data[pos].epsilon  = EPS0;
		data[pos].mu       = MUE0;
		data[pos].rhoprime = 1e40;
		data[pos].sigma    = 1e40;
	}*/


	const double imp0 = 377;
	for ( int timestep=0; timestep < 40; ++timestep )
	{
		const int X = 0;
		const int Y = 1;
		const int Z = 2;

		const int tcur  = 0;
		const int tnext = 0;
		
		VecI pos(0); pos[X] = 0;
		for ( int y = 0; y < N_CELLS_Y; y++ ) {
			pos[Y] = y;
			data[pos].E[0][Z] = t_spawn_func( timestep );
		}
		
		for ( int x = 0; x < N_CELLS_X; x++ )
		for ( int y = 0; y < N_CELLS_Y; y++ )
		for ( int z = 0; z < 1; z++ ) {
			VecI pos; pos[X]=x; pos[Y]=y; pos[Z]=z;

			VecI xprev=pos; xprev[X]--;
			VecI xnext=pos; xnext[X]++;
			
			/* Now update all E components */
			double Hytcxn = x   < N_CELLS_X ? data[pos  ].H[0][Y] : 0;
			double Hytcxp = x-1 >= 0        ? data[xprev].H[0][Y] : 0;
			data[pos].E[0][Z] = data[pos].E[0][Z] + (Hytcxn - Hytcxp) * imp0;
		}
		
		for ( int x = 0; x < N_CELLS_X; x++ )
		for ( int y = 0; y < N_CELLS_Y; y++ )
		for ( int z = 0; z < 1; z++ ) {
			VecI pos; pos[X]=x; pos[Y]=y; pos[Z]=z;

			VecI xprev=pos; xprev[X]--;
			VecI xnext=pos; xnext[X]++;

			/* Update all H components */
			double Eztcxn = x+1 < N_CELLS_X ? data[xnext].E[0][Z] : 0;
			double Eztcxp = x   >= 0        ? data[pos  ].E[0][Z] : 0;
			data[pos].H[0][Y] = data[pos].H[0][Y] + (Eztcxn - Eztcxp ) / imp0;
		}
		
		std::cout << "E_z after Timestep: " << timestep << std::endl;
		for ( int ix=0; ix < N_CELLS_X; ++ix ) {
		    VecI pos; pos[X]=ix; pos[Y]=0; pos[Z]=0;
		    std::cout << data[pos].E[0][Z] << " ";
		}
		std::cout << std::endl;
		std::cout << "H_y after Timestep: " << timestep << std::endl;
		for ( int ix=0; ix < N_CELLS_X; ++ix ) {
		    VecI pos; pos[X]=ix; pos[Y]=0; pos[Z]=0;
		    std::cout << data[pos].H[tnext][Y] << " ";
		}
		std::cout << std::endl;

		if (timestep % 2 == 0) {
			/* Beware! PNGWriter begins counting with 1 in image coordinates */
			static int framecounter = 0;
			framecounter++;
			char filename[100];
			{
			sprintf( filename, "output/Ex_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, 0.5+double(ix % 2), 0.5+double(ix % 2), 0.5+double(ix % 2));
			}
			image.close();
			}{
			sprintf( filename, "output/Ey_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].E[tnext][Y], -data[pos].E[tnext][Y], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Ez_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].E[tnext][Z], -data[pos].E[tnext][Z], 0.0);
			}
			image.close();
			}

			{
			sprintf( filename, "output/Hx_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].H[tnext][X], -data[pos].H[tnext][X], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Hy_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].H[tnext][Y], -data[pos].H[tnext][Y], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Hz_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos; pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].H[tnext][Z], -data[pos].H[tnext][Z], 0.0);
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