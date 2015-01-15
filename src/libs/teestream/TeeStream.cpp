TeeStream::TeeStream( void ) { }

TeeStream::TeeStream( int verbosity ) { 
    this->verbosity = verbosity;
}

TeeStream::TeeStream( const char * filename ) {
    filestream.open( filename, std::ofstream::out | std::ofstream::app );
}

void TeeStream::Open( std::string filename, int rank ) {
    /* Create Timestamp and rank strings for filenames */
    time_t t = time(0);   // get time now
    struct tm * now = localtime( &t );
    std::stringstream timestamp;
    timestamp << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
              << now->tm_mday << "_" << now->tm_hour << "-"
              << now->tm_min << "_";
    std::string timesr = timestamp.str();
    std::stringstream ranksr; ranksr << rank;

    filestream.open( timesr + filename + std::string("_rank_") 
                     + ranksr.str() + std::string(".log"),
                     std::ofstream::out | std::ofstream::app );
}

TeeStream::~TeeStream(void) {
    if ( filestream.is_open() )
        filestream.close();
}

template <class T> TeeStream& TeeStream::operator<<(T val) {
    if ( not filestream.is_open() ) {
        std::cerr << "You need to specify a filename with tout.open( "
            "filename, rank ), before writing to TeeStream object.\n";
        return *this;
    }
    filestream << val;
    if ( this->verbosity == 1 )
        std::cout << val;
    filestream.flush();
    return *this;
}

TeeStream& TeeStream::operator<< (std::ostream& (*pfun)(std::ostream&)) {
    if ( not filestream.is_open() ) {
        std::cerr << "You need to specify a filename with tout.open( "
            "filename, rank ), before writing to TeeStream object.\n";
        return *this;
    }
    pfun(filestream);
    if ( this->verbosity == 1 )
        pfun(std::cout);
    return *this;
}
