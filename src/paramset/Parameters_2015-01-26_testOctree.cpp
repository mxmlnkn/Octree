#pragma once

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>   // FLT_EPSILON
#include "math/TVector.h"

#ifndef M_PI
#   define M_PI 3.14159265358979323846
#endif
#define INF (1.0/0.0)

//componentsConfig.param
const uint16_t SIMDIM                  = 2;
const bool     PERIODIC_FORCE          = false;     // if false CONSIDERATION_RATIO still has effect, except if it is higher than the largest possible distance in one cell
const double   CONSIDERATION_RATIO     = 100.;       // Force between particles will be calculated if distance between those two is smaller than CONSIDERATION_RATIO * CELL_SIZE. If 0 then force is not calculated periodically
      int      NUMBER_OF_STEPS         = 1e6;
const bool     FULL_N_SQUARED_FORCE    = false;    // if false, then force is calculated only with the particles in the same cell
const uint32_t RANDOM_SEED             = 43907340; // if 0, then time() is used
const bool     STOPPING                = false;    // if some momentum ist lost every time step. Can be used to determine state of lowest energy

//GasConfig.param
const double GAS_DENSITY_SI            = 1.e30;    // 1/m^3

//physicalConstants.param
const double ELEMENTARY_CHARGE_SI      = 1.602176565e-19;  // C
const double PROTON_MASS_SI            = 1.6726231e-27;    // kg
const double SPEED_OF_LIGHT_SI         = 2.99792458e8;     // m/s
const double MUE0_SI                   = M_PI * 4.e-7;     // N/A^2 = kg*m/C^2
const double EPS0_SI                   = 1.0/(MUE0_SI*SPEED_OF_LIGHT_SI*SPEED_OF_LIGHT_SI);    // C^2/J*m, 8.854187817e-12

//GridConfig.param
      int GUARDSIZE                    = 1;
      int NUMBER_OF_CELLS_X            = 50; // these is for the INITIAL_OCTREE_REFINEMENT ...
      int NUMBER_OF_CELLS_Y            = 50;  // ... real number of cells will be higher
      int NUMBER_OF_CELLS_Z            = 50;    // ... at critical places
      int NUMBER_OF_PARTICLES_PER_CELL = 26;   // NUM instead of NUMBER_OF in picongpu is also inconsistent, and there are other longer names, soo ...
const double CELL_SIZE_SI              = 1e-9; // m
const double CELL_SIZE_X_SI            = CELL_SIZE_SI;   // (!!!) picongpu naming with width, height and depth seems to be too random => could lead to mixups
const double CELL_SIZE_Y_SI            = CELL_SIZE_SI;
const double CELL_SIZE_Z_SI            = 1./0.;
const double DELTA_T_SI                = 0.6*CELL_SIZE_SI / SPEED_OF_LIGHT_SI;
      int BOUNDARY_CONDITION           = 1;              //0:periodic, 1:reflecting, 2:adhering
const int isPeriodic[3]                = {0,0,0};

//particleConfig.param
      int    SPECIES                   = 3;                    // 1:only electrons, 2:only ions, 3:both in equal parts (NUMBER_OF_PARTICLES_PER_CELL must be even!)
const double ELECTRON_TEMPERATURE_keV  = 0.200;                // keV !
const double ION_TEMPERATURE_keV       = 0.100;                // keV
const double ELECTRON_MASS_SI          = 9.109382e-31;         // kg
const double ELECTRON_CHARGE_SI        =-ELEMENTARY_CHARGE_SI; // C
const double ION_MASS_SI               = ELECTRON_MASS_SI;     // kg
const double ION_CHARGE_SI             = ELEMENTARY_CHARGE_SI; // C
const double PARTICLE_INIT_DRIFT_GAMMA = 1.0;                  // 1.0 means v = 0
const double PARTICLE_RADIUS_SI        = 0.01*CELL_SIZE_SI; //0.84087e-15;          // m (proton_radius used here) //normally in PIC: CELL_SIZE_SI / 2.0;, where CELL_SIZE ~ DEBYE_LENGTH
const uint16_t DEFAULT_PARTICLE_SHAPE  = 1;  //00:point-point, 01:ball-ball (CIC radialsymmetric equivalent), 99:sphere-sphere (all except 00 will use ION_RADIUS)

//output
int PNG_INTERVAL            = 100;
int PRINT_INTERVAL          = 100;
int PRINTF_INTERVAL         = 10;
int PRINTF_SIMDATA_INTERVAL = std::min( int(ceil( 1e-18 / DELTA_T_SI )), 1 );
// for dt = 1e-19 => Nprint = 10. for 1e-17 it prints every time step, so that we can see something

//LASER:
const double LAMBDA_SI                 = 20*CELL_SIZE_X_SI;

//COLORS:


//================================== Units ===================================//
const double UNITCONV_keV_to_Joule     = 1.e3 * ELEMENTARY_CHARGE_SI;
const double UNITCONV_Joule_to_keV     = 1.0 / UNITCONV_keV_to_Joule;

const double UNIT_LENGTH               = CELL_SIZE_SI;
const double UNIT_MASS                 = ELECTRON_MASS_SI;
const double UNIT_CHARGE               = ELEMENTARY_CHARGE_SI;
const double UNIT_TIME                 = DELTA_T_SI;

//derived units
const double UNIT_VELOCITY             = UNIT_LENGTH / UNIT_TIME;
const double UNIT_ENERGY               = UNIT_MASS * UNIT_LENGTH*UNIT_LENGTH / (UNIT_TIME*UNIT_TIME); // in general not eV !
const double UNIT_MOMENTUM             = UNIT_MASS * UNIT_VELOCITY;
const double UNIT_ANGULAR_MOMENTUM     = UNIT_LENGTH * UNIT_MASS * UNIT_VELOCITY;
const double UNIT_FORCE                = UNIT_MASS * UNIT_LENGTH / (UNIT_TIME*UNIT_TIME);

//physicalConstants.unitless   //In Picongpu the variable names change to Q_EL, M_EL and so on, which I find inconsistent -.- (!!!)
const double ELEMENTARY_CHARGE         = ELEMENTARY_CHARGE_SI / UNIT_CHARGE;  // C
const double SPEED_OF_LIGHT            = SPEED_OF_LIGHT_SI / UNIT_VELOCITY;  // m/s
const double ELECTRON_CHARGE           = ELECTRON_CHARGE_SI / UNIT_CHARGE;
const double ELECTRON_MASS             = ELECTRON_MASS_SI / UNIT_MASS;
const double MUE0                      = MUE0_SI / (UNIT_MASS * UNIT_LENGTH / (UNIT_CHARGE*UNIT_CHARGE));    //! magnetic constant must be double 3.92907e-39
const double EPS0                      = 1. / ( MUE0*SPEED_OF_LIGHT*SPEED_OF_LIGHT );  //! electric constant must be double 2.54513e+38

//particleConfig.unitless
const double ELECTRON_TEMPERATURE_SI   = ELECTRON_TEMPERATURE_keV * UNITCONV_keV_to_Joule; // J
const double ELECTRON_TEMPERATURE      = ELECTRON_TEMPERATURE_SI / UNIT_ENERGY;
const double ION_TEMPERATURE_SI        = ION_TEMPERATURE_keV * UNITCONV_keV_to_Joule;      // J
const double ION_TEMPERATURE           = ION_TEMPERATURE_SI / UNIT_ENERGY;
const double ION_CHARGE                = ION_CHARGE_SI / UNIT_CHARGE;
const double ION_MASS                  = ION_MASS_SI / UNIT_MASS;
const double PARTICLE_RADIUS           = PARTICLE_RADIUS_SI / UNIT_LENGTH;
      int NUMBER_OF_PARTICLES          = NUMBER_OF_PARTICLES_PER_CELL * NUMBER_OF_CELLS_X * NUMBER_OF_CELLS_Y * NUMBER_OF_CELLS_Z;
      int NUMBER_OF_EONS_PER_CELL      = (SPECIES == 1) ? NUMBER_OF_PARTICLES_PER_CELL : ( (SPECIES == 2) ? 0 : NUMBER_OF_PARTICLES / 2 );
      int NUMBER_OF_IONS_PER_CELL      = (SPECIES == 2) ? NUMBER_OF_PARTICLES_PER_CELL : ( (SPECIES == 1) ? 0 : NUMBER_OF_PARTICLES / 2 );

//GridConfig.unitless
const double DELTA_T                   = DELTA_T_SI / UNIT_TIME;
const double CELL_SIZE_X               = CELL_SIZE_X_SI / UNIT_LENGTH;
const double CELL_SIZE_Y               = CELL_SIZE_Y_SI / UNIT_LENGTH;
const double CELL_SIZE_Z               = CELL_SIZE_Z_SI / UNIT_LENGTH;
Vec<double,SIMDIM> CELL_SIZE           = Vec<double,SIMDIM>( CELL_SIZE_X, CELL_SIZE_Y, CELL_SIZE_Z );
const double CELL_SIZE_MIN             = std::min( CELL_SIZE_X, std::min( CELL_SIZE_Y, CELL_SIZE_Z ) );
Vec<int,SIMDIM> NUMBER_OF_CELLS        = Vec<int,SIMDIM>( NUMBER_OF_CELLS_X, NUMBER_OF_CELLS_Y, NUMBER_OF_CELLS_Z );
Vec<double,SIMDIM> SIM_SIZE            = Vec<double,SIMDIM>( NUMBER_OF_CELLS ) * CELL_SIZE;

//LASER:
const double LAMBDA                    = LAMBDA_SI / UNIT_LENGTH;

//SETUP:
int INITIAL_OCTREE_REFINEMENT          = 4; // yiels pow(INIT_..., [4,8]) octree cells minimum
int MAX_OCTREE_REFINEMENT              = 6; // yiels pow(INIT_..., [4,8]) octree cells minimum
int OCTREE_SETUP                       = 6;

int SIMULATION_SETUP                   = 4;
double ABSORBING_BORDER_THICKNESS      = 10*CELL_SIZE_X + FLT_EPSILON;
double ABSORBER_STRENGTH               = 2e4;
Vec<double,SIMDIM> M                   = 0;    // center of circle
double R                               = 0.3*SIM_SIZE.min(); // radius of circle
double MIRROR_WIDTH                    = SIM_SIZE[1] - 2*ABSORBING_BORDER_THICKNESS;
double MIRROR_LENGTH                   = SIM_SIZE[0]/8;
Vec<double,SIMDIM> MIRROR_CENTER       = Vec<double,SIMDIM>( 2*ABSORBING_BORDER_THICKNESS, 0.5*SIM_SIZE[1] );
double FOCAL_LENGTH                    = pow(MIRROR_WIDTH/2, 2) / (4*MIRROR_LENGTH);

int WAVE_SPAWN_SETUP                   = 3;
Vec<double,SIMDIM> SPAWN_AREA_SIZE     = Vec<double,SIMDIM>( 3*CELL_SIZE_X, 3*LAMBDA );
Vec<double,SIMDIM> SPAWN_POS           = ABSORBING_BORDER_THICKNESS + Vec<double,SIMDIM>( 0.0, M[1] + 0.6*R );
//Vec<double,SIMDIM> SPAWN_POS           =  MIRROR_CENTER + Vec<double,SIMDIM>( FOCAL_LENGTH, - 0.5 );
Vec<double,SIMDIM> WAVE_GUIDE_SIZE     = 0;//Vec<double,SIMDIM>( 7*LAMBDA, 15.0 );
