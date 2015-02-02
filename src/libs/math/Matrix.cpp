/* delete matrix from heap */
void Matrix::Clear(void) {
    for (int i=0; i<m; i++) {
        delete[] data[i];
    }
    delete[] data;

    #if DEBUG >= 2
        cout << "(" << m << "," << n << ")-Matrix data cleared\n";
    #endif
    return;
}

void Matrix::Setup(int pm, int pn) {
    this->m = pm; this->n = pn;
    data = new double* [m];
    for (int i=0; i<m; i++) {
        data[i] = new double[n];
        for (int j=0; j<n; j++)
            data[i][j]=0;
    }
    return;
}

/* Create new 2d-array with (m,n) being the dimension */
Matrix::Matrix(int pm, int pn) : m(0), n(0), data(NULL) {
    Setup(pm,pn);
    return;
}

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

template<int T_DIM> Matrix::Matrix(const Vec<double,T_DIM> v) : m(0), n(0), data(NULL)
{
    Setup(T_DIM,1);
    for (int i=0; i < m*n; i++)
        (*this)[i] = v[i];
}

Matrix::Matrix(const Matrix &mat) : m(0), n(0), data(NULL) {
    Setup(mat.m,mat.n);
    *this = mat;
    #if DEBUG >= 2
        cout << "Because only the matrix pointer was COPIED(constructor) prior and therefore the 2darray was deleted twice!\n";
    #endif
    return;
}

Matrix::Matrix(double x) : m(0), n(0), data(NULL) {
    Setup(1,1);
    (*this)(0,0)=x;
    return;
}

Matrix::~Matrix() {                                          //delete matrix from heap
    Clear();
    #if DEBUG >= 2
        cout << "(" << m << "," << n << ")-Matrix deleted!\n";
    #endif
    return;
}

template<int T_DIM>
Matrix::operator Vec<double,T_DIM> () const
{
    Vec<double,T_DIM> v;
    for (int i=0; i < T_DIM; i++)
        if ( i < m*n )
            v[i] = (*this)[i];
        else
            v[i] = 0;
    return v;
}

Matrix& Matrix::operator=(const Matrix &mat) {
    #if DEBUG >= 2
        cout << "Assigment of following in process:\n";
        MDIM dim = mat.GetDim();
        cout << this->m << "," << this->n << "\n";
        mat.Print();
    #endif
    if (this == &mat)                                       //Test for Self Assignment
        return *this;

    this->Clear();
    Setup(mat.m,mat.n);

    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            data[i][j] = mat(i,j);
    return *this;
}

Matrix& Matrix::operator=(const double val) {
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            data[i][j] = val;
    return *this;
}

template<int T_DIM>
Matrix& Matrix::operator=(const Vec<double,T_DIM> val)
{
    for (int i=0; i < std::min(T_DIM,m*n); i++)
        (*this)[i] = val[i];
    return *this;
}

double& Matrix::operator()(int i, int j) const{
    assert( i<m and j<n );
    return data[i][j];
}

double& Matrix::operator[](int x) const {
    if ( (m!=1 && x>=m) || (n!=1 && x>=n) || !IsVector() ) {
        #if DEBUG == 1
            cout << "Element reference [] out of bound or not Vector\n";
        #endif
        return data[0][0];
    }

    if (m==1)
        return data[0][x];
    else
        return data[x][0];
}

Matrix Matrix::operator+(const Matrix &mat) {
    if (m != mat.m || n != mat.n)
        return *this;

    Matrix sum(m,n);
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            sum(i,j) = mat(i,j) + (*this)(i,j);
        }
    }

    #if DEBUG >= 2
        MDIM dim = sum.GetDim();
        cout << "Sum-dim: " << dim.m << "," << dim.n << "\n";
    #endif
    return sum;
}

Matrix Matrix::operator-(const Matrix &mat) {          //Subtract
    if (m != mat.m || n != mat.n)
        return *this;

    Matrix diff(m,n);
    #if DEBUG >= 2
        cout << this->m << "," << this->n << "\n";
    #endif
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            diff(i,j) = (*this)(i,j) - mat(i,j);
        }
    }
    return diff;
}

Matrix Matrix::operator*(const Matrix &mat) {   //this * mat
    assert( n == mat.m );

    Matrix product(m,mat.n);
    for (int i=0; i<m; i++) {                   //go down every line
        for (int j=0; j<mat.n; j++) {           //go down every column
            double sum=0;
            for (int k=0; k<n; k++) {           //go down every column in *this and every line in mat
	           sum += (*this)(i,k) * mat(k,j);
            }
            product(i,j) = sum;
        }
    }
    return product;
}

Matrix Matrix::operator*(const double a) const
{
    Matrix product(m,n);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            product(i,j) = (*this)(i,j) * a;
    return product;
}

inline Matrix Matrix::operator/(double divisor) {                 //scalar Division
    return (*this)*(1/divisor);
}

Matrix Matrix::operator*(double factor) {                         //scalar Division
    Matrix product(m,n);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            product(i,j) = factor * (*this)(i,j);
    return product;
}

bool Matrix::operator==(const Matrix &mat) {
    if ( m!=mat.m || n!=mat.n )
        return false;
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            if ( (*this)(i,j) != mat(i,j) )
                return false;
    return true;
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

Matrix Matrix::Transpose(void) const {
    Matrix trans(n,m);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            trans(j,i) = (*this)(i,j);
    return trans;
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

double Matrix::Trace(void) const {
    if (!IsSquare())
        return -1;

    double sum;
    for (int i=0; i<m; i++)
        sum += (*this)(i,i);
    return sum;
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

inline bool Matrix::IsSquare(void) const{
    return m==n;
}

inline bool Matrix::IsVector(void) const {
    return (m==1 || n==1);
}

int Matrix::GetVectorDim(void) const {
    if (!IsVector())
        return -1;
    else if (m == 1)
        return n;
    else
        return m;
}

int Matrix::GetSquareDim(void) const {
    if (IsSquare())
        return m;
    else
        return -1;                                  //-1 errorcode and clearly not the dimension of Matrix
}




void Matrix::SetDiagonal(double val) {
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            if (i==j)
                (*this)(i,i)=val;
            else
                (*this)(i,j)=0;
    return;
}

void Matrix::SetDiagonal(double val[]) {
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            if (i==j)
                (*this)(i,i)=val[i];
            else
                (*this)(i,j)=0;
    return;
}

void Matrix::SetAll(double val) {
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            (*this)(i,j)=val;
    return;
}



void Matrix::DelRow(int row, int amount) {                          //Deletes row-th Row counting from 1 and all n after that
    Matrix mat(m-amount,n);
    int k1=0;
    const int i=row-1;

    for (int k0=0; k0<m; k0++) {
        if ( (k0 < i)  || (k0 >= i+amount) ) {
            for (int j=0; j<n; j++) {
                mat(k1,j) = data[k0][j];            //WHY THE FUCK is gcc warning me about this line when it works for example in operator+ -.-?! Well as long as it works: [Warning] passing `double' for converting 1 of `double& Matrix::operator()(int, int) const'
            }
            k1++;
        }
    }
    *this = mat;
    return;
}

void Matrix::DelCol(int column, int amount) {
    Matrix mat(m,n-amount);
    int l1=0;
    const int j=column-1;

    for (int i=0; i<m; i++) {
        for (int l0=0; l0<n; l0++) {
            if ( (l0 < j) || (l0 >= j+amount) ) {
                mat(i,l1) = data[i][l0];            //WHY THE FUCK is gcc warning me about this line when it works for example in operator+ -.-?! Well as long as it works: [Warning] passing `double' for converting 1 of `double& Matrix::operator()(int, int) const'
                l1++;
            }
        }
        l1=0;
    }
    *this = mat;
    return;
}

Matrix Matrix::Augment(const Matrix &mat) const {
    if (m != mat.m)
        return *this;

    Matrix merg(m,n+mat.n);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            merg(i,j) = data[i][j];
    for (int i=0; i<m; i++)
        for (int j=0; j<mat.n; j++)
            merg(i,n+j) = mat(i,j);

    return merg;
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

MDIM Matrix::GetDim(void) const {
    MDIM dim;
    dim.m=m;
    dim.n=n;
    return dim;
}

void Matrix::Print(void) const {
    #if DEBUG >= 2
        std::cout << "Printing: " << m << "x" << n << "\n";
    #endif
    for (int i=0; i < m; i++) {
        std::cout << "[";
        for (int j=0; j < n; j++) {
            std::cout << std::setw(4) << data[i][j] << " ";
        }
        std::cout << "]\n";
    }
    return;
}


template<typename T>
Matrix operator*( const T scalar, const Matrix & rhs )
{
    return rhs * scalar;
}
