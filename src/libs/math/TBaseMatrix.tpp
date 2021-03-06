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

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( VecI sizeIn ) : size(sizeIn), data(NULL) {
    this->data = new T_DTYPE[this->size.product()];
}

/****************************** Copy Constructor ******************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::BaseMatrix( const BaseMatrix & m ) : size(m.size), data(NULL) {
    this->data = new T_DTYPE[m.size.product()];
    /* use this instead of memcpy, to call all constructors, if T_DTYPE is
     * e.g. a class! */
    for ( int i=0; i < size.product(); ++i )
        this->data[i] = m[i];
}

/********************************* Destructor *********************************/
template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM>::~BaseMatrix( void ) {
    if ( this->data != NULL )
        delete[] this->data;
}

/***************************** Assignment Operator ****************************/
template<typename T_DTYPE, int T_DIM>
template<typename T_ETYPE>
BaseMatrix<T_DTYPE,T_DIM> & BaseMatrix<T_DTYPE,T_DIM>::operator=
(const T_ETYPE a)
{
    for (int i=0; i < size.product(); i++)
        this->data[i] = a;
    /* it's possible, this is called recursively, if the elements of          *
     * BaseMatrix are again BaseMatrix, but the recursion will stop if        *
     * T_DTYPE is reached. e.g. BaseMatrix< BaseMatrix< double,3 >, 3 >       */
    return *this;
}

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM> & BaseMatrix<T_DTYPE,T_DIM>::operator=
(const BaseMatrix m)
{
    this->size = m.size;
    if ( this->data != NULL )
        delete[] this->data;
    this->data = new T_DTYPE[m.size.product()];

    for ( int i=0; i < size.product(); ++i )
        this->data[i] = m[i];
    return *this;
}

template<typename T_DTYPE, int T_DIM>
int BaseMatrix<T_DTYPE,T_DIM>::getLinearIndex( VecI pos ) const
{
    return ConvertVectorToLinearIndex( pos, this->size );
}

template<typename T_DTYPE, int T_DIM>
Vec<int,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getVectorIndex( int linindex ) const
{
    return ConvertLinearToVectorIndex( linindex, this->size );
}

/****************************** Access Operators ******************************/
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

/******************************* Other Functions ******************************/

template<typename T_DTYPE, int T_DIM>
typename BaseMatrix<T_DTYPE,T_DIM>::VecI
BaseMatrix<T_DTYPE,T_DIM>::getSize( void ) const
{
    return this->size;
}

template<typename T_DTYPE, int T_DIM>
void BaseMatrix<T_DTYPE,T_DIM>::setSize( VecI psize )
{
    BaseMatrix mnew(psize);

    /* copy old values into new matrix ... */
    for ( int i = 0; i < mnew.size.product(); ++i )
    {
        VecI vind = mnew.getVectorIndex(i);
        /* ... if that index already existed in the old matrix */
        if ( vind < this->size )
            /*better than memcpy, because it calls T_DTYPE assignment operator*/
            mnew[vind] = (*this)[vind];
    }

    (*this) = mnew;
}

template<typename T_DTYPE, int T_DIM>
BaseMatrix<T_DTYPE,T_DIM> BaseMatrix<T_DTYPE,T_DIM>::getPartialMatrix
( const VecI & pos, const VecI & sizePartialMatrix ) const
{
    assert( pos + sizePartialMatrix <= this->size );
    BaseMatrix<T_DTYPE,T_DIM> tmpPartialMatrix( sizePartialMatrix );
    VecI ind(0);

    for ( int i=0; i < sizePartialMatrix.product(); i++ )
        tmpPartialMatrix[i] = (*this)[ pos + tmpPartialMatrix.getVectorIndex(i) ];

    return tmpPartialMatrix;
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

/****************************** Resizing methods ******************************/

template<typename T_DTYPE, int T_DIM>
T_DTYPE InterpolateWithLagrange( const std::vector<Vec<double,T_DIM> > & x0,
    const std::vector<T_DTYPE> & y0, Vec<double,T_DIM> x )
{
    assert( x0.size() == y0.size() );
    int n = x0.size();

    T_DTYPE value(0);
    for ( int i=0; i<n; i++ ) {
        T_DTYPE summand(y0[i]);
        for ( int j=0; j<n; j++ ) {
            if ( i == j )
                continue;
            summand *= (x-x0[j]).norm() / (x0[i]-x0[j]).norm();
        }
        value += summand;
    }
    return value;
}

template<typename T_DTYPE, int T_DIM>
void BaseMatrix<T_DTYPE,T_DIM>::LagrangianResizeTo( BaseMatrix & target, double order ) const {
    assert(order > 0);

    for ( int i=0; i<target.size.product(); i++ ) {
        /* Find values we want to use to interpolate current target cell */
        VecD pos = ( VecD(0.5) + VecD(target.getVectorIndex(i)) ) /
                   VecD(target.size) * VecD(this->size);
        /* Cycle through all cells of source matrix ... maybe something better? */
        std::vector<Vec<double,T_DIM> > positions;
        std::vector<T_DTYPE> values;
        for ( int j=0; j<size.product(); j++ ) {
            VecD posSource = VecD(0.5) + VecD(getVectorIndex(j));
            if ( (posSource - pos).norm() <= order ) {
                positions.push_back(posSource);
                values.push_back( (*this)[j] );
            }
        }
        target[i] = InterpolateWithLagrange( positions, values, pos );
    }
}

template<typename T_DTYPE, int T_DIM>
void BaseMatrix<T_DTYPE,T_DIM>::NearestResizeTo( BaseMatrix & target ) const {
    assert( target.size.product() != 0);
    for ( int i=0; i<target.size.product(); i++ ) {
        /* Find values we want to use to interpolate current target cell */
        VecD pos = ( VecD(0.5) + VecD(target.getVectorIndex(i)) ) /
                   VecD(target.size) * VecD(this->size);
        /* Cycle through all cells of source matrix ... maybe something better? */
        int lastInd = -1;
        double lastDist = 1e5;
        for ( int j=0; j<size.product(); j++ ) {
            VecD posSource = VecD(0.5) + VecD(getVectorIndex(j));
            double curDist = (posSource - pos).norm();
            if ( curDist < lastDist ) {
                lastInd = j;
                lastDist = curDist;
            }
        }
        assert( lastInd != -1 );
        target[i] = (*this)[lastInd];
    }
}

/* Enables cout << BaseMatrix<int,2>(1); This also works with fstream and therefore with cout */
template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out, const BaseMatrix<T_DTYPE,T_DIM>& m ) {
    typedef typename BaseMatrix<T_DTYPE,T_DIM>::VecI VecI;
    VecI size = m.getSize();
    VecI index(0);
    bool isvector = size[0] == 1 or size[1] == 1 ? 1 : 0;

    for ( int i = 0; i < T_DIM - isvector; ++i )
        out << "{";

    for ( int i = 0; i < size.product(); ++i )
    {
        out << m[index];
        int closedBrackets = 0 - isvector;

        /* increment vector index */
        index[T_DIM-1] += 1;
        for ( int j=T_DIM-1; j>=0; j-- )
        {
            if ( index[j] >= size[j] )
            {
                if ( j > 0 ) {
                    index[j] = 0;
                    index[j-1] += 1;
                } else
                    assert( i == size.product()-1 );
                closedBrackets++;
            }
        }

        for ( int j = 0; j < closedBrackets; ++j )
            out << "}";

        /* if not last element, then add ',' and open brackets for next element */
        if ( i != size.product()-1 )
        {
            out << ",";
            for ( int j = 0; j < closedBrackets; ++j )
                out << "{";
        }
    }
    return out;
}
