template<int T_DIM>
int ConvertVectorToLinearIndex( const Vec<int,T_DIM> & pos, const Vec<int,T_DIM> & size ) {
    assert( pos < size );
    assert( pos >= (Vec<int,T_DIM>)(0) );
    int index     = 0;
    int prevrange = 1;
    for (int i=T_DIM-1; i>=0; i--) {
        index     += pos[i] * prevrange;
        prevrange *= size[i];
    }
    assert( index < size.product() );
    return index;
}

template<int T_DIM>
Vec<int,T_DIM> ConvertLinearToVectorIndex( const int & linindex, const Vec<int,T_DIM> & size ) {
    assert( linindex < size.product() );
    Vec<int,T_DIM> index;
    int tmp = linindex;
    for (int i=T_DIM-1; i>=0; i--) {
        index[i]  = tmp % size[i];
        tmp       = tmp / size[i];
    }
    assert( tmp == 0 );
    assert( index <  size );
    assert( index >= (Vec<int,T_DIM>)(0) );
    return index;
}

/******************************** Constructor *********************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( void ) : size(0), data(NULL)
{}

/****************************** Copy Constructor ******************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( const BaseMatrix & m ) : size(m.size), data(NULL) {
    this->data = new T_DTYPE[size.product()];
    memcpy( this->data, m.data, sizeof(T_DTYPE) * m.size.product() );
}

/********************************* Destructor *********************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::~BaseMatrix( void ) {
    if ( this->data != NULL )
        delete[] this->data;
}

/***************************** Assignment Operator ****************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>& BaseMatrix<T_DTYPE,T_DIM>::operator= (const BaseMatrix & m)
{
    this->size = m.size;
    if ( this->data != NULL )
        delete[] this->data;
    this->data = new T_DTYPE[m.size.product()];

    memcpy( this->data, m.data, sizeof(T_DTYPE) * m.size.product() );
    return *this;
}

/***************************** Other Constructors *****************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( VecI size ) : size(size), data(NULL) {
    this->data = new T_DTYPE[size.product()];
}

/************************** Assignment Operators **************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>& BaseMatrix<T_DTYPE,T_DIM>::operator= (const T_DTYPE a) {
    for (int i=0; i < size.product(); i++)
        this->data[i] = a;
    return *this;
}

template<typename T_DTYPE, int T_DIM>
int BaseMatrix<T_DTYPE,T_DIM>::getLinearIndex( const VecI & pos ) const {
    return ConvertVectorToLinearIndex( pos, this->size );
}

template<typename T_DTYPE, int T_DIM>
Vec<int,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getVectorIndex( const int & linindex ) const {
    return ConvertLinearToVectorIndex( linindex, this->size );
}

/**************************** Access Operators ****************************/
template<typename T_DTYPE, int T_DIM>
T_DTYPE BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const int i ) const {
    assert( i < size.product() );
    return data[i];
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE & BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const int i ) {
    assert( i < size.product() );
    return data[i];
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const VecI pos ) const {
    assert( pos < size );
    return data[ getLinearIndex(pos) ];
}

template<typename T_DTYPE, int T_DIM>
T_DTYPE & BaseMatrix<T_DTYPE,T_DIM>::operator[] ( const VecI pos ) {
    assert( pos < size );
    return data[ getLinearIndex(pos) ];
}


/***************************** Other Functions ****************************/
template<typename T_DTYPE, int T_DIM>
Vec<int,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getSize( void ) const {
    return this->size;
}

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getPartialMatrix( const VecI & pos, const VecI & size ) const {
    assert( pos+size <= this->size );
    BaseMatrix<T_DTYPE,T_DIM> tmp( size );
    VecI ind(0);

    for ( int i=0; i<size.product(); i++ )
        tmp[i] = (*this)[ pos + tmp.getVectorIndex(i) ];

    return tmp;
}

template<typename T_DTYPE, int T_DIM>
void BaseMatrix<T_DTYPE,T_DIM>::insertMatrix( const VecI & pos, const BaseMatrix & m ) {
    assert( pos+m.getSize() <= this->size );
    for ( int i=0; i < m.getSize().product(); i++ ) {
        VecI index = pos + m.getVectorIndex(i);
        assert( index < this->size );
        (*this)[ index ] = m[i];
    }
}


/* Enables cout << Vec<int,2>(1); This alos works with fstream and therefore with cout */
template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out, const BaseMatrix<T_DTYPE,T_DIM>& m ) {
    typedef typename BaseMatrix<T_DTYPE,T_DIM>::VecI VecI; // xD ... typisch picongpu
    VecI size = m.getSize();
    out << "This " << size << "Matrix:" << std::endl;
    VecI ind(0);
    for ( ind[1]=0; ind[1]<size[1]; ind[1]++) {
        for ( ind[0]=0; ind[0]<size[0]; ind[0]++)
            out << m[ ind ] << " ";
        out << std::endl;
    }
    out << std::endl;
    return out;
}
