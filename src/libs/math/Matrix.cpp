#pragma once


/******************************** Constructors ********************************/

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( void )
 : BaseMatrix<T_DTYPE,2>()
{ };

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( int m, int n )
 : BaseMatrix<T_DTYPE,2>( VecI(m,n) )
{ }

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( VecI psize )
 : BaseMatrix<T_DTYPE,2>( VecI( psize[0], psize[1] ) )
{ }

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::MathMatrix( const MathMatrix & m )
 : BaseMatrix<T_DTYPE,2>( m )
{ }

template<typename T_DTYPE>
MathMatrix<T_DTYPE>::~MathMatrix( void )
/* no explicit call for parent class destructor needed */
{ }

/**************************** Assignment operators ****************************/
// broadcast value to all elements
template<typename T_DTYPE>
template<typename T_ETYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator= ( T_ETYPE a )
{
    BaseMatrix<T_DTYPE,2>::operator=(a);
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator= (const MathMatrix m)
{
    /* Can't call copy assignment from base-class, because it expects a
     * BaseMatrix argument instead of MathMatrix. We would have to convert it
     * first */
    BaseMatrix<T_DTYPE,2>::operator=( (BaseMatrix<T_DTYPE,2>) m );
    return *this;
}

/************************ Binary Arithmetic Operators *************************/

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator+=(const MathMatrix &mat)
{
    assert( this->size == mat.size );
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] += mat[i];
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator-=(const MathMatrix &mat)
{
    assert( this->size == mat.size );
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] -= mat[i];
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator*=(const MathMatrix &mat)
{
    *this = (*this) * mat;
    return *this;
}

template<typename T_DTYPE>
template<typename T_ETYPE>
inline MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator/=(T_ETYPE a)
{
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] /= a;
    return *this;
}

template<typename T_DTYPE>
template<typename T_ETYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::operator*=(T_ETYPE a)
{
    for ( int i = 0; i < this->size.product(); ++i )
        (*this)[i] *= a;
    return *this;
}

/**************************** Arithmetic Operators ****************************/

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator+(const MathMatrix &mat) const
{
    MathMatrix<T_DTYPE> res( *this );
    res += mat;
    return res;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator-(const MathMatrix &mat) const
{
    MathMatrix<T_DTYPE> res( *this );
    res += mat;
    return res;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator*(const MathMatrix &mat) const
{
    assert( this->size[1] == mat.size[0] );
    int m = this->size[0];
    int n = mat.size[1];
    MathMatrix res(m,n);
    /* copy example element from mat to all elements in res, thereby applying *
     * dimensions, if T_DTYPE is a matrix and not double                      */
    res = mat(0,0)*0; /* implicit broadcast */

    for ( int j = 0; j < res.size[1]; j++ ) {       // go down every column
        for ( int i = 0; i < res.size[0]; i++ ) {   // go down every line
            T_DTYPE sum = res(i,j)*0; // copy res(i,j), to to copy the dimension, this also works, if T_DTYPE is a double!
            for ( int k = 0; k < this->size[1]; k++ ) { // dot product
	           sum += (*this)(i,k) * mat(k,j);
            }
            res(i,j) = sum;
        }
    }
    return res;
}

template<typename T_DTYPE>
template<typename T_ETYPE>
inline MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator/(T_ETYPE divisor) const
{
    return (*this) * ( 1.0 / divisor );
}

template<typename T_DTYPE>
template<typename T_ETYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::operator*(T_ETYPE a) const
{
    MathMatrix<T_DTYPE> res( this->size );
    for ( int i = 0; i < this->size.product(); ++i )
        res[i] = (*this)[i] * a;
    return res;
}

template<typename T_DTYPE>
bool MathMatrix<T_DTYPE>::operator==(const MathMatrix &mat) const
{
    if ( this->size != mat.size ) {
        return false;
    } else {
        for ( int i = 0; i < this->size.product(); ++i )
            if ( (*this)[i] != mat[i] )
                return false;
    }
    return true;
}

template<typename T_DTYPE>
inline bool MathMatrix<T_DTYPE>::operator!=(const MathMatrix &mat) const
{
    return not ( (*this) == mat );
}

/******************************** Test methods ********************************/

template<typename T_DTYPE>
inline bool MathMatrix<T_DTYPE>::isSquare( void ) const
{
    return this->size[0] == this->size[1];
}

template<typename T_DTYPE>
inline bool MathMatrix<T_DTYPE>::isVector( void ) const
{
    return ( this->size[0] == 1 || this->size[1] == 1 );
}

/****************************** Modifying methods *****************************/

/* k=0 is main diagonal, k>0 is above main diagonal, k<0 is below */
template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::setDiagonal( T_DTYPE a, int k )
{
    if ( k >= 0 ) {
        assert( k < this->size[0] );
        for ( int i = 0; i < this->size[1] - k; ++i )
            (*this)( i,i+k ) = a;
    } else {
        assert( -k < this->size[1] );
        for ( int i = 0; i < this->size[0] + k; ++i )
            (*this)( i-k,i ) = a;
    }
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::setDiagonal( T_DTYPE a[], int k )
{
    if ( k >= 0 ) {
        assert( k < this->size[0] );
        for ( int i = 0; i < this->size[1] - k; ++i )
            (*this)( i,i+k ) = a[i];
    } else {
        assert( -k < this->size[1] );
        for ( int i = 0; i < this->size[0] + k; ++i )
            (*this)( i-k,i ) = a[i];
    }
    return *this;
}

template<typename T_DTYPE>
template<typename T_ETYPE>
inline MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::setAll( T_DTYPE a ) {
    (*this) = a;
    return *this;
}

template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::setSize( int m, int n )
{
    BaseMatrix<T_DTYPE,2>::setSize( VecI(m,n) );
    return *this;
}

/******************************* Reading methods*******************************/

template<typename T_DTYPE>
T_DTYPE & MathMatrix<T_DTYPE>::operator()( int i, int j )
{
    return (*this)[ VecI(i,j) ];
}

template<typename T_DTYPE>
T_DTYPE MathMatrix<T_DTYPE>::operator()( int i, int j ) const
{
    return (*this)[ VecI(i,j) ];
}

template<typename T_DTYPE>
inline int MathMatrix<T_DTYPE>::getVectorDim(void) const {
    if ( this->size[1] == 1 )
        return this->size[0];
    else if ( this->size[0] == 1 )
        return this->size[1];
    else {
        assert(false);
        return -1;
    }
}

/* -1 errorcode and clearly not the dimension of Matrix */
template<typename T_DTYPE>
inline int MathMatrix<T_DTYPE>::getSquareDim(void) const {
    if ( isSquare() )
        return this->size[0];
    else {
        assert(false);
        return -1;
    }
}

/************************** Simple Matrix Arithmetics *************************/

template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::transpose( void ) const
{
    /* could try to do this by setting a flag... but that could be error prone*/
    MathMatrix trans( this->size[0], this->size[1] );
    for ( int i=0; i < this->size[0]; i++ )
        for (int j = 0; j < this->size[1]; j++ )
            trans(j,i) = (*this)(i,j);
    return trans;
}

template<typename T_DTYPE>
T_DTYPE MathMatrix<T_DTYPE>::trace( void ) const
{
    if ( not isSquare() ) {
        assert(false);
        return -1;
    }

    /*copy from this(0,0), in order to be correct if T_DTYPE is a larger class*/
    T_DTYPE sum = (*this)(0,0)*0;
    for ( int i = 0; i < this->size[0]; ++i )
        sum += (*this)(i,i);
    return sum;
}

/* Deletes row-th Row counting from 1 and all n after that */
template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::delRow( int irow, int nrows )
{
    int m = this->size[0];
    int n = this->size[1];

    assert( irow > 0 and irow < m );
    assert( irow + nrows <= m );
    MathMatrix res( m - nrows, n );

    /* we need a special i-counter, because when deleting rows this counter   *
     * should not be incremented                                              */
    int ires = 0;

    /* copy all rows ... */
    for ( int i = 0; i < m; ++i ) {
        /* ... except those to be deleted */
        if ( i < irow or i >= irow + nrows ) {
            /* copy all elements column-wise in row i */
            for ( int j = 0; j < n; ++j )
                res(ires,j) = (*this)(i,j);
            ires++;
        }
    }

    *this = res;
    return *this;
}

/* deletes ncols beginning with column icol */
template<typename T_DTYPE>
MathMatrix<T_DTYPE> & MathMatrix<T_DTYPE>::delCol(int icol, int ncols)
{
    int m = this->size[0];
    int n = this->size[1];

    assert( icol > 0 and icol <n );
    assert( icol + ncols <= n );
    MathMatrix res( m, n - ncols );

    /* we need a special i-counter, because when deleting rows this counter   *
     * should not be incremented                                              */
    int jres = 0;

    /* copy all columns ... */
    for ( int j = 0; j < n; ++j ) {
        /* ... except those to be deleted */
        if ( j < icol or j >= icol + ncols ) {
            /* copy all elements row-wise in column j */
            for ( int i = 0; i < m; ++i )
                res(i,jres) = (*this)(i,j);
            jres++;
        }
    }

    *this = res;
    return *this;
}

/*template<typename T_DTYPE>
MathMatrix<T_DTYPE> MathMatrix<T_DTYPE>::augment(const MathMatrix &mat) const {
    if (m != mat.m)
        return *this;

    MathMatrix merg(m,n+mat.n);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            merg(i,j) = data[i][j];
    for (int i=0; i<m; i++)
        for (int j=0; j<mat.n; j++)
            merg(i,n+j) = mat(i,j);

    return merg;
}*/

/********************** More advanced Matrix arithmetics **********************/
#if 1==0
/* Load from File? */
Matrix::Matrix(char* filename) : m(0), n(0), data(NULL) {
    std::ifstream file(filename);
    char s[10];
    int spaces=0, eols=0, line=0;
    do {
        file.getline(s, 10);
        for (int i=0; i<10; i++)
            if (s[i]==' ')
                spaces++;
        line++;
    } while(!file.eof() && file.fail() );

    if (file.fail()) {
        #if DEBUG >= 1
            cout << "Failure when determining number of columns\n";
        #endif
        return;
    }

    eols++;
    char * row = (char*) malloc(10*line*sizeof(char));
    while (!file.eof()) {
        file.getline(row,10*line);
        if (file.fail()) {                                  //TODO: if this happens, then the char-array is too short (but too lazy to fix this)
            #if DEBUG >= 1
                cout << "Error determining number of Rows in File (non standard conform file?) \"" << filename << "\"\n";
            #endif
            return;
        }
        eols++;
    };
    #if DEBUG >= 2
        cout << "Determined from \"" << filename << "\": Rows: " << eols << " Columns: " << spaces+1 << "\n";
    #endif
    file.close();

    Setup(eols,spaces+1);
    Load(filename);
    free(row);
    return;
}

double Matrix::Minor(int row, int col) const {                  //i counting from 1 (mathematical)!
    if ( !IsSquare() || (row*col<=0) || (row>m) || (col>m) )    //Det and Minor only possible on square matrices
        return 0;                                               //not the wisest error code :S
    Matrix mat = *this;
    mat.DelRow(row);
    mat.DelCol(col);
    return mat.Det();
}

double Matrix::Det(void) const {
    if (!IsSquare())                                            //Det only possible with square matrices
        return 0;
    if (m==1)                                                   //Determinant of 1x1 Matrix = one and only element
        return data[0][0];

    //recursive call Det->Minor->Det->... (Laplace expansion along the first column)
    double sum=0;
    const int j=1;                                              //first column
    for (int i=1; i<=m; i++) {
        if ( (i+j)%2 == 0 ) {
            sum += data[i-1][j-1] * Minor(i,j);
        } else {
            sum -= data[i-1][j-1] * Minor(i,j);
        }
    }
    return sum;
}

Matrix Matrix::Invert(void) const {
    if (!IsSquare())                                                //only possible with square matrices
        return *this;
#define JACOBI 1
#if CRAMER == 1
    return ( ((*this).Adjugate()) * (1/(*this).Det()) );
#elif JACOBI == 1
    Matrix id(m,m);
    id.SetDiagonal(1);
    Matrix inv = (this->Augment(id)).RowEchelon();

    //inv is now in row echelon form -> make the upper right corner zero
    for (int i=inv.m-1; i>=0; i--) {                        //Begin at m(-1 because array counts from 0)
        //Subtract i-th line from j-th line above
        for (int j=i-1; j>=0; j--) {
            #if DEBUG >=2
                cout << "inv(" << j << "," << i << ")(" << inv(i,j) << ") = factor\n";
            #endif
            //Subtract i-th line from this j-th line with factor inv(...,i) (because with the last line we make the rightmost row zero)
            double factor = inv(j,i);
            for (int k=0; k<inv.n; k++) {                   //begin with k, because the lower left corner is zero so that the subtrahend equals 0 in this area
                inv(j,k) -= factor * inv(i,k);
            }
        }
    }

    inv.DelCol(1,m);                                        //Delete augmented part
    return inv;
#endif
}

Matrix Matrix::Adjugate(void) const {
    if (m != n)
        return *this;
    Matrix adj(m,n);
    for (int i=1; i<=m; i++) {
        for (int j=1; j<=n; j++) {
            #if DEBUG >=2
                cout << "a[" << i << "," << j << "]: ";
            #endif
            if ((i+j)%2 ==0) {
                adj(i-1,j-1) = Minor(j,i);
            } else {
                adj(i-1,j-1) = -Minor(j,i);
            }
        }
    }
    return adj;
}

double Matrix::Norm(void) const {
    if (!IsVector()) {
        #if DEBUG == 1
            cout << "Can't calculate Euclidic Norm of Non-Vector\n";
        #endif
        return -1;                                              //well chosen errorcode, because Norm is non-negative
    }
    double sum=0;
    int d = GetVectorDim();
    for (int i=0; i<d; i++) {
        sum += (*this)[i] * (*this)[i];
    }

    return sqrt(sum);
}

Matrix Matrix::Abs(void) const {
    Matrix res(m,n);
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            res(i,j) = fabs( (*this)(i,j) );
        }
    }
    return res;
}

int Matrix::Rank(void) const {
    Matrix mat = RowEchelon();
    int rank = (m<n) ? (m) : (n);
    for (int i=rank-1; i>=0; i--)
        for (int j=0; j<n; j++) {
            #if DEBUG >= 2
                cout << "(" << i << "," << j << ")=" << abs(mat(i,j)) << "\n";
            #endif
            if ( fabs(mat(i,j)) != 0.0 )
                return i+1;
        }
    return 0;
}

Matrix Matrix::RowEchelon (void) const {             //TODO
    Matrix echelon = (*this);
    for (int i=0; i<echelon.m; i++) {
        //Set Diagonal value to 1 by dividing complete row with inv(i,i)
        double divisor = (fabs(echelon(i,i)) == 0.0) ? (1) : (echelon(i,i));       //abs(divisor), because 0 != -0

        for (int j=0; j<echelon.n; j++) {
            #if DEBUG >=2
                cout << "inv(" << i << "," << j << ")(" << echelon(i,j) << ") /= " << echelon(i,i) << "\n";
            #endif
            echelon(i,j) /= divisor;
        }

        //Subtract Multiple of i-th line
        for (int j=i+1; j<echelon.m; j++) {
            double factor = echelon(j,i);
            for (int k=0; k<echelon.n; k++) {
                echelon(j,k) -= factor * echelon(i,k);
            }
        }
    }
    return echelon;
}



int Matrix::Load(char* filename) {                      //returns not 0 if there was an error
    //misses dummy matrix if file reading doesn't work in the middle of overwriting the old matrix
    std::ifstream mf(filename);                              //mf = matrix file
    if (!mf) {
        #if DEBUG == 1
            cout << "Couldn't open " << filename << "\n";
        #endif
        return 1;
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            mf >> data[i][j];
            if (mf.bad())
                break;
        }
    }
    if (mf.bad()) {
        #if DEBUG == 1
            cout << "Reading error in " << filename << "\n";
        #endif
        return 1;
    }

    mf.close();
    return 0;
}

void Matrix::Save(char* filename) {
    std::ofstream file(filename);
    for (int i=0; i < m; i++) {
        for (int j=0; j < n; j++) {
            file << data[i][j];
            if (j < n-1)
                file << " ";
        }
        if (i < m-1)
            file << "\n";
    }
    file.close();
    return;
}
#endif

template<typename T_DTYPE, typename T_ETYPE>
MathMatrix<T_DTYPE> operator*( const T_ETYPE a, const MathMatrix<T_DTYPE> & rhs )
{
    return rhs * a;
}
