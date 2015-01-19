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

struct {
    Vec<double,3> operator() ( YeeCell cell ) {
        Vec<double,3> color;
        color[0] = cell.E[2];
        color[1] =-cell.E[2];
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
        double gray = 0;
        if (cell.epsilon > EPS0)
            gray = cell.epsilon / EPS0 / 10.0;
        return Vec<double,3>(gray);
    }
} returnn;
