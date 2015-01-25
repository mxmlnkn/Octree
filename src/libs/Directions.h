#pragma once

#include "math/TVector.h"
#include <assert.h>

/**************************************************************************
 * 3 bits are not enough, they give only 8 possibilities. That wouldn't   *
 * include diagonal neighbors of all sorts (edges, corners in 3D)         *
 * Therefore we need tertiary math:                                       *
 *   0=>-1 from position, 1=>stay at pos, 2=>+1 from position along axis  *
 * So the axes in 3D become:                                              *
 *   axis 0: -1=left  , 0=stay, +1=right                                  *
 *   axis 1: -3=bottom, 0=stay, +3=top    (bottom+right <! stay)          *
 *   axis 2: -9=back  , 0=stay, +9=front                                  *
 * Difference between two states of axis n must be smaller than range     *
 * accessible by alle axes < n                                            *
 * We get the following structure:                                        *
 *                                                       back:            *
 *                                    stay:         --------------        *
 *                  front:        --------------   | 11 | 12 | 13 |       *
 *              --------------   |  2 |  3 |  4 |  |--------------|       *
 *       top:  | -7 | -6 | -5 |  |--------------|  |  8 |  9 | 10 |       *
 *             |--------------|  | -1 |  0 |  1 |  |--------------|       *
 *      stay:  |-10 | -9 | -8 |  |--------------|  |  5 |  6 |  7 |       *
 *             |--------------|  | -4 | -3 | -2 |   --------------        *
 *    bottom:  |-13 |-12 |-11 |   --------------                          *
 *              --------------                                            *
 *               ^    ^    ^                                              *
 *             left  stay right                                           *
 * As can be seen arbitrary diagonal neighbors can be easily derived by   *
 * adding instead von anding the values, like it would be the case for    *
 * a bitmask:                                                             *
 *   left + top = -1-+3 = 2   =>  correct                                 *
 * by shifting with +13 we can map this to positive numbers               *
 * These numbers although nicely to look at, proove problematic (tried to *
 * access array[-1]. It is better to have only positive values, which is  *
 * possible if we say: negative value -> 2*abs( negative value ). We get  *
 * the following structure:                                               *
 *   axis 0: 2=left  , 0=stay, +1=right                                   *
 *   axis 1: 6=bottom, 0=stay, +3=top    (bottom+right <! stay)           *
 *   axis 2: 18=back , 0=stay, +9=front                                   *
 *                                                       back:            *
 *                                    stay:         --------------        *
 *                  front:        --------------   | 23 | 21 | 22 |       *
 *              --------------   |  5 |  3 |  4 |  |--------------|       *
 *       top:  | 14 | 12 | 13 |  |--------------|  | 20 | 18 | 19 |       *
 *             |--------------|  |  2 |  0 |  1 |  |--------------|       *
 *      stay:  | 11 |  9 | 10 |  |--------------|  | 26 | 24 | 25 |       *
 *             |--------------|  |  8 |  6 |  7 |   --------------        *
 *    bottom:  | 17 | 15 | 16 |   --------------                          *
 *              --------------                                            *
 *               ^    ^    ^                                              *
 *             left  stay right                                           *
 **************************************************************************/

/* These are the same access-specifiers like in SimulationBox.h -> merge? */
const int X_AXIS = 0;
const int Y_AXIS = 1;
const int Z_AXIS = 2;

/* This basically is a 3x3 Matrix => use getLinearIndex from              *
 * BaseMatrix to calculate these directions on the fly. It should anyway  *
 * not be necessary for N Dimensions to use LEFT,RIGHT,...                */
const int LEFT   = 2;
const int RIGHT  = 1;
const int BOTTOM = 6;
const int TOP    = 3;
const int BACK   = 18;
const int FRONT  = 9;

/* These should be changed to using BaseMatrix::getVectorIndex.           *
 * Problem: Will have to check for (0,0,0,...) at every for loop !!!      *
 *   - as all the values would be positive we also should shift it by the *
 *     position of the center, to get directions relative to it           */
template<int T_DIM>
Vec<int,T_DIM> getDirectionVector( const int direction ) {
    Vec<int,T_DIM> vec(0);
    assert( direction >= 1 and direction <= pow(3,T_DIM)+0.0001 );
    /* This may not be needed for 2D, but it also doesn't make the        *
     * results wrong                                                      */
    int t_direction = direction;
    switch( t_direction % 3 ) {
        case 1: vec[X_AXIS] = +1; break;
        case 2: vec[X_AXIS] = -1; break;
    }
    t_direction /= 3;
    switch( t_direction % 3 ) {
        case 1: vec[Y_AXIS] = +1; break;
        case 2: vec[Y_AXIS] = -1; break;
    }
    t_direction /= 3;
    switch( t_direction % 3 ) {
        case 1: vec[Z_AXIS] = +1; break;
        case 2: vec[Z_AXIS] = -1; break;
    }
    return vec;
}

/* reverse of the above getDirectionVector function */
template<int T_DIM>
int getLinearDirection( const Vec<int,T_DIM> & v ) {
    int direction = 0;
    int prevrange = 1;
    for (int i=0; i<T_DIM; i++) {
        int value = v[i];
        if (value==-1) value = 2;
        direction += value * prevrange;
        prevrange *= 3;
    }
    assert( direction >=1 and direction <= pow(3,v.dim)+0.0001 );
    assert( getDirectionVector<T_DIM>( direction ) == v );
    return direction;
}

template<int T_DIM>
int getOppositeDirection( const int direction ) {
    return getLinearDirection( getDirectionVector<T_DIM>( direction )*(-1) );
}

/* returns the axis corresponding to a direction:                         *
 *   RIGHT,LEFT -> 0, TOP,BOTTOM -> 1, FRONT,BACK -> 2, ...               *
 * returns -1 if e.g. RIGHT+LEFT                                          */
template<int T_DIM>
int getAxis( const int direction) {
    Vec<int,T_DIM> v( getDirectionVector<T_DIM>( direction ) );
    int axis = -1;
    for (int i=0; i<v.dim; i++)
        if ( v[i] != 0 ) {
            if (axis == -1)
                axis = i;
            else
                /* 2nd axis found -> on diagonal axes return -1 */
                return -1;
        }
    return axis;
}
