#pragma once

#include <list>

#ifndef M_PI
#   define M_PI 3.14159265358979323846
#endif
#define INF (1.0/0.0)

namespace compileTime {

/* Compile time power (also exact for integers) */
template<typename T>
inline constexpr T pow(const T base, unsigned const exponent) {
    return exponent == 0 ? 1 : base * pow<T>(base, exponent-1);
}

} // compileTime

/* Some simplifications when using std::list */

#define CONTAINS(list,value) (std::find(list.begin(), list.end(), value) != list.end() )

template<typename T_DTYPE>
std::ostream& operator<<( std::ostream& out, const std::list<T_DTYPE> ls ) {
    out << "{";
    for ( typename std::list<T_DTYPE>::const_iterator it = ls.begin();
    it != ls.end(); ++it )
        out << *it << ",";
    out << "}";
    return out;
}
