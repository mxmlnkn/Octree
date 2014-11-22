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
        std::ofstream fileStream;
        TeeStream( void ) { }
        TeeStream( int verbosity = 1 ) { 
            this->verbosity = verbosity;
        }
        TeeStream( const char * filename ) {
            fileStream.open( filename, std::ofstream::out | std::ofstream::app );
        }
        void Open( std::string filename, int rank = 0 ) {
            /* Create Timestamp and rank strings for filenames */
            time_t t = time(0);   // get time now
            struct tm * now = localtime( &t );
            std::stringstream timestamp;
            timestamp << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
                      << now->tm_mday << "_" << now->tm_hour << "-"
                      << now->tm_min << "_";
            std::string timesr = timestamp.str();
            std::stringstream ranksr; ranksr << rank;

            fileStream.open( timesr + filename + std::string("_rank_") 
                             + ranksr.str() + std::string(".log"),
                             std::ofstream::out | std::ofstream::app );
        }
        ~TeeStream(void) {
            fileStream.close();
        }
        template <class T> TeeStream& operator<<(T & val) {
            fileStream << val;
            if ( this->verbosity == 1 )
                std::cout << val;
            fileStream.flush();
            return *this;
        }
        TeeStream& operator<< (std::ostream& (*pfun)(std::ostream&)) {
            pfun(fileStream);
            if ( this->verbosity == 1 )
                pfun(std::cout);
            return *this;
        }
};

TeeStream tout(1);
TeeStream terr(0);
