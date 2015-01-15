#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>

/* outstream similar to cout but outputs (also/only) into files (tout/terr) */
class TeeStream{
    public:
        int verbosity;
        std::ofstream filestream;
        TeeStream( void );
        TeeStream( int verbosity = 1 );
        TeeStream( const char * filename );
        void Open( std::string filename, int rank = 0 );
        ~TeeStream(void);
        /* for things like tout << "hi" << 3; */
        template <class T> TeeStream& operator<<(T val);
        /* for things like tout << std::flush */
        TeeStream& operator<< (std::ostream& (*pfun)(std::ostream&));
};

#include "TeeStream.cpp"

TeeStream tout(1);
TeeStream terr(0);
