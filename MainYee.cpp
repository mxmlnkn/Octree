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

#define DEBUG_MAIN_YEE 99

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
#define YEE_CELL_TIMESTEPS_TO_SAVE 2
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
	for ( int x = N_CELLS_X/2; x < N_CELLS_X/2+10; x++ )
	for ( int y = 0; y < N_CELLS_Y; y++ ) {
		const int wy = 40;
		if ( y <=  N_CELLS_Y/2 - wy/2 or y >= N_CELLS_Y/2 + wy/2 ) {
			VecI pos(GUARDSIZE); pos[0]=x; pos[1]=y;
			data[pos].epsilon  = 2*EPS0;
			//data[pos].mu       = INF;
		}
	}


	for ( int timestep=0; timestep < 400; ++timestep )
	{
		const int X = 0;
		const int Y = 1;
		const int Z = 2;

		/* For YEE_CELL_TIMESTEPS_TO_SAVE == 2 this switches 0 and 1 as the   *
		 * saving time step                                                   */
        int tcur  = (timestep  ) % YEE_CELL_TIMESTEPS_TO_SAVE;
        int tnext = (timestep+1) % YEE_CELL_TIMESTEPS_TO_SAVE;
		/*tcur  = 0;
		tnext = 0; */

		/* Function Generator on left side creates sine wave */
		VecI pos(GUARDSIZE); pos[X] = GUARDSIZE; pos[Z] = GUARDSIZE;
		for ( int y = 0; y < N_CELLS_Y - 2*GUARDSIZE; y++ ) {
			pos[Y] = y + GUARDSIZE;
			data[pos].E[tcur][Z] = t_spawn_func( timestep );
		}

		/**********************************************************************
		 * Initialize y-Guard with copied values to make boundaries periodic  *
		 *             -------------                                          *
		 *            |G G G G G G G|    Ny = 7, NGuard = 2                   *
		 *            |G G G G G G G|    periodic y:                          *
		 *  y     -------------  G G|     - x = 0...Nx, z = 0...Nz            *
		 *     0 |G G G G G G G| G G|     - [x,3] -> [x,0]                    *
		 *  ^  1 |G G G G G G G| G G|     - [x,4] -> [x,1]                    *
		 *  |  2 |G G C C C G G| G G|     - [x,2] -> [x,5]                    *
		 *  |  3 |G G C C C G G| G G|     - [x,3] -> [x,6]                    *
		 *  |  4 |G G C C C G G|----                                          *
		 *  |  5 |G G G G G G G|      [x,Ny-1-2*NGuard+0] -> [x,0]            *
		 *  |  6 |G G G G G G G|      [x,Ny-1-2*NGuard+1] -> [x,1]            *
		 *  |     -------------       [x,NGuard+0]        -> [x,Ny-NGuard+0]  *
		 *  |     0 1 2 3 4 5 6       [x,NGuard+1]        -> [x,Ny-NGuard+1]  *
		 *  -----------------> x                                              *
		 **********************************************************************/
		for ( int x = 0; x < N_CELLS_X; x++ )
		for ( int z = 0; z < N_CELLS_Z; z++ )
		for ( int y = 0; y < GUARDSIZE; y++ ) {
			VecI posTo(0); posTo[X]=x; posTo[Z]=z;
			VecI posFrom(posTo);

			/* [x,Ny-1-2*NGuard+0] -> [x,0] */
			posFrom[Y] = N_CELLS_Y-1-2*GUARDSIZE+y;
			posTo[Y]   = y;
			data[posTo].E[tcur] = data[posFrom].E[tcur];

			#if DEBUG_MAIN_YEE >= 100
            if (timestep == 3)
                std::cout << "(" << "Ez.[" << posFrom << "]=" << data[posFrom].E[tcur][Z] << " -> "
                << "(" << "Ez.[" << posTo << "]=" << data[posTo].E[tcur][Z] << std::endl;
			#endif

			/* [x,NGuard+0] -> [x,Ny-NGuard+0] */
			posFrom[Y] = GUARDSIZE+y;
			posTo  [Y] = N_CELLS_Y-GUARDSIZE+y;
			data[posTo].E[tcur] = data[posFrom].E[tcur];

			#if DEBUG_MAIN_YEE >= 100
            if (timestep == 3)
                std::cout << "(" << "Ez.[" << posFrom << "]=" << data[posFrom].E[tcur][Z] << " -> "
                << "(" << "Ez.[" << posTo << "]=" << data[posTo].E[tcur][Z] << "\n\n";
			#endif
		}


		#if DEBUG_MAIN_YEE >= 100
		std::cout << "E_z[x=GUARDSIZE,y=1...6...122...127,z=GUARDSIZE] before Timestep and after updating periodic boundary guards: " << timestep << std::endl;
		for ( int ix=0; ix < 8; ++ix ) {
			for ( int iy=0; iy < 3; ++iy ) {
				VecI pos(GUARDSIZE); pos[X]+=ix; pos[Y]=iy;
				std::cout << data[pos].E[tcur][Z] << " ";
			}
			std::cout << "... ";
			for ( int iy=125; iy < N_CELLS_Y; ++iy ) {
				VecI pos(GUARDSIZE); pos[X]+=ix; pos[Y]=iy;
				std::cout << data[pos].E[tcur][Z] << " ";
			}
			std::cout << std::endl;
		}
		#endif

		for ( int x = 0; x < N_CELLS_X - 2*GUARDSIZE; x++ )
		for ( int y = 0; y < N_CELLS_Y - 2*GUARDSIZE; y++ )
		for ( int z = 0; z < N_CELLS_Z - 2*GUARDSIZE; z++ ) {
			VecI pos;
			pos[X]=x + GUARDSIZE;
			pos[Y]=y + GUARDSIZE;
			pos[Z]=z + GUARDSIZE;

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
			double nom = (1.0 + 0.5*data[pos].rhoprime*DELTA_T / data[pos].mu);
			double Da  = (1.0-0.5*data[pos].rhoprime*DELTA_T/data[pos].mu)/nom;
			double Dbx = DELTA_T / ( data[pos].mu * CELL_SIZE_X ) / nom;
			double Dby = DELTA_T / ( data[pos].mu * CELL_SIZE_Y ) / nom;
			double Dbz = DELTA_T / ( data[pos].mu * CELL_SIZE_Z ) / nom;

			#if DEBUG_MAIN_YEE >= 100
				if ( data[ynext].E[tcur] != data[pos].E[tcur] )
					std::cout << "!!!! E[" << ynext << "]=" << data[ynext].E[tcur] << " != E[" << pos << "]=" << data[pos].E[tcur] << "\n";
			#endif

			data[pos].H[tnext][X] = Da * data[pos].H[tcur][X] +
				 Dbz*( data[znext].E[tcur][Y]-data[pos].E[tcur][Y] )
				-Dby*( data[ynext].E[tcur][Z]-data[pos].E[tcur][Z] );
			data[pos].H[tnext][Y] = Da * data[pos].H[tcur][Y] +
				 Dbx*( data[xnext].E[tcur][Z]-data[pos].E[tcur][Z] )
				-Dbz*( data[znext].E[tcur][X]-data[pos].E[tcur][X] );
			data[pos].H[tnext][Z] = Da * data[pos].H[tcur][Z] +
				 Dby*( data[ynext].E[tcur][X]-data[pos].E[tcur][X] )
				-Dbx*( data[xnext].E[tcur][Y]-data[pos].E[tcur][Y] );
		}

		/* Reinitialize y-guard, after updating H everywhere! */
		for ( int x = 0; x < N_CELLS_X; x++ )
		for ( int z = 0; z < N_CELLS_Z; z++ )
		for ( int y = 0; y < GUARDSIZE; y++ ) {
			VecI pos(0); pos[X]=x; pos[Z]=z;
			pos[Y] = y;
			VecI mirrorpos(pos); mirrorpos[Y] = pos[Y] + (N_CELLS_Y-1-GUARDSIZE);
			data[pos].H[tnext] = data[mirrorpos].H[tnext];
			pos[Y] = N_CELLS_Y-1 - y;
			mirrorpos[Y] = pos[Y] - (N_CELLS_Y-GUARDSIZE);
			data[pos].H[tnext] = data[mirrorpos].H[tnext];
		}

		/* Now update all E components */
		for ( int x = 0; x < N_CELLS_X - 2*GUARDSIZE; x++ )
		for ( int y = 0; y < N_CELLS_Y - 2*GUARDSIZE; y++ )
		for ( int z = 0; z < N_CELLS_Z - 2*GUARDSIZE; z++ ) {
			VecI pos;
			pos[X]=x + GUARDSIZE;
			pos[Y]=y + GUARDSIZE;
			pos[Z]=z + GUARDSIZE;

			VecI xprev=pos; xprev[X]--;
			VecI xnext=pos; xnext[X]++;

			VecI yprev=pos; yprev[Y]--;
			VecI ynext=pos; ynext[Y]++;

			VecI zprev=pos; zprev[Z]--;
			VecI znext=pos; znext[Z]++;

			double nom = (1.0 + 0.5*data[pos].sigma*DELTA_T / data[pos].epsilon);
			double Ca = (1.0 - 0.5*data[pos].sigma*DELTA_T / data[pos].epsilon)
					  / nom;
			double Cbx = DELTA_T / ( data[pos].epsilon * CELL_SIZE_X ) / nom;
			double Cby = DELTA_T / ( data[pos].epsilon * CELL_SIZE_Y ) / nom;
			double Cbz = DELTA_T / ( data[pos].epsilon * CELL_SIZE_Z ) / nom;

			data[pos].E[tnext][X] = Ca * data[pos].E[tcur][X] +
				 Cby*( data[pos].H[tnext][Z]-data[yprev].H[tnext][Z] )
				-Cbz*( data[pos].H[tnext][Y]-data[zprev].H[tnext][Y] );
			data[pos].E[tnext][Y] = Ca * data[pos].E[tcur][Y] +
				 Cbz*( data[pos].H[tnext][X]-data[zprev].H[tnext][X] )
				-Cbx*( data[pos].H[tnext][Z]-data[xprev].H[tnext][Z] );
			data[pos].E[tnext][Z] = Ca * data[pos].E[tcur][Z] +
				 Cbx*( data[pos].H[tnext][Y]-data[xprev].H[tnext][Y] )
				-Cby*( data[pos].H[tnext][X]-data[yprev].H[tnext][X] );
		}

		/* Reinitialize y-guard, after updating H everywhere! */
		for ( int x = 0; x < N_CELLS_X; x++ )
		for ( int z = 0; z < N_CELLS_Z; z++ )
		for ( int y = 0; y < GUARDSIZE; y++ ) {
			VecI pos(0); pos[X]=x; pos[Z]=z;
			pos[Y] = y;
			VecI mirrorpos(pos); mirrorpos[Y] = pos[Y] + (N_CELLS_Y-1-GUARDSIZE);
			data[pos].E[tnext] = data[mirrorpos].E[tnext];
			pos[Y] = N_CELLS_Y-1 - y;
			mirrorpos[Y] = pos[Y] - (N_CELLS_Y-GUARDSIZE);
			data[pos].E[tnext] = data[mirrorpos].E[tnext];
		}

		#if DEBUG_MAIN_YEE>=90
		std::cout << "E_z after Timestep: " << timestep << std::endl;
		for ( int ix=0; ix < 8; ++ix ) {
		    VecI pos(GUARDSIZE); pos[X]=ix;
		    std::cout << data[pos].E[tnext][Z] << " ";
		}
		std::cout << std::endl;
		std::cout << "H_y after Timestep: " << timestep << std::endl;
		for ( int ix=0; ix < 8; ++ix ) {
		    VecI pos(GUARDSIZE); pos[X]=ix;
		    std::cout << data[pos].H[tnext][Y] << " ";
		}
		std::cout << std::endl;
		#endif

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
				VecI pos(GUARDSIZE); pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].E[tnext][X], -data[pos].E[tnext][X], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Ey_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos(GUARDSIZE); pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].E[tnext][Y], -data[pos].E[tnext][Y], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Ez_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos(GUARDSIZE); pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].E[tnext][Z], -data[pos].E[tnext][Z], 0.0);
				/* Output Blende */
				if (data[pos].epsilon > EPS0) {
					double gray = data[pos].epsilon / EPS0 / 10.0;
					image.plot_blend( ix+1,iy+1, 0.5, gray, gray, gray);
				}
			}
			image.close();
			}

			{
			sprintf( filename, "output/Hx_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos(GUARDSIZE); pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].H[tnext][X], -data[pos].H[tnext][X], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Hy_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos(GUARDSIZE); pos[0]=ix; pos[1]=iy;
				image.plot( ix+1,iy+1, data[pos].H[tnext][Y], -data[pos].H[tnext][Y], 0.0);
			}
			image.close();
			}{
			sprintf( filename, "output/Hz_%05i.png", framecounter );
			pngwriter image( N_CELLS_X,N_CELLS_Y, 0, filename );
			for ( int ix=0; ix < N_CELLS_X; ++ix )
			for ( int iy=0; iy < N_CELLS_Y; ++iy ) {
				VecI pos(GUARDSIZE); pos[0]=ix; pos[1]=iy;
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