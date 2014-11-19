#pragma once

namespace SimulationBox {

    /* This are predefined constant for addressing. E.g. if we want to        *
     * iterate only over the CORE or we want a value on the X_AXIS. These are *
     * are not sizes or widths!                                               */
    const int GUARD  = 1;
    const int BORDER = 2;
    const int CORE   = 4;

    const int X_AXIS = 0;
    const int Y_AXIS = 1;
    const int Z_AXIS = 2;

}