/* Some functors for PNG-Output */

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        color[0] = cell.E[0];
        color[1] = cell.H[0];
        color[2] =-cell.E[0];
        return color;
    }
} returnEandH;

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        color[0] = cell.E[0];
        color[1] = cell.E[0];
        return color;
    }
} returnEx;

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        color[0] = cell.E[1];
        color[1] =-cell.E[1];
        return color;
    }
} returnEy;

#define BLACKGROUND 0

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        #if BLACKGROUND == 1
            color[0] = cell.E[2];
            color[1] =-cell.E[2];
        #else
            color[0] = cell.E[2] > 0 ? 1.0                 : 1.0 + cell.E[2];
            color[1] = cell.E[2] < 0 ? 1.0                 : 1.0 - cell.E[2];
            color[2] = cell.E[2] > 0 ? 1.0 - 0.5*cell.E[2] : 1.0 + cell.E[2];
        #endif
        return color;
    }
} returnEz;

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        color[0] = cell.H[0];
        color[1] =-cell.H[0];
        return color;
    }
} returnHx;

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        color[0] = cell.H[1];
        color[1] =-cell.H[1];
        return color;
    }
} returnHy;

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        color[0] = cell.H[2];
        color[1] =-cell.H[2];
        return color;
    }
} returnHz;

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        return Vec<double,3>(cell.epsilon / EPS0 / 10.0 + cell.sigmaE / ABSORBER_STRENGTH);
    }
} returnn;
