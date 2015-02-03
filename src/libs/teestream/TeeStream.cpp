TeeStream::TeeStream( int verbosityIn )
: verbosity(verbosityIn), filestream()
{}

void TeeStream::Open( std::string filename, int rank, std::string folder ) {
    std::stringstream fullname;
    fullname << folder << filename;
    if ( rank != -1 )
        fullname << "_rank_" << rank;
    fullname << ".log";
    filestream.open( fullname.str(), std::ofstream::out | std::ofstream::app );
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
