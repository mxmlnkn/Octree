#pragma once

template<typename T_COEFF>
std::vector<std::string> Polynomial<T_COEFF>::splitByDelimiter
( std::string src, char delimiter ) const
{
    std::vector<std::string> res;

    /* search for delimiter */
    if ( delimiter == -1 )
        for ( auto it = src.begin(); it != src.end(); ++it )
            if ( not isalnum( *it ) )
                delimiter = *it;


    /* split src by delimiter */
    size_t pos1 = 0;
    while (true)
    {
        src = src.substr( pos1 );
        pos1 = src.find( delimiter );

        if ( pos1 != src.npos ) {
            std::string svar = src.substr( 0, pos1 );
            if ( svar.length() != 0 )
                res.push_back( svar );
            ++pos1;
        }
        else
        {
            pos1 = src.length() - 1;
            break;
        }
    }
    if ( src.length() != 0 )
        res.push_back( src );

    return res;
}

/******************************** Constructors ********************************/

template<typename T_COEFF>
Polynomial<T_COEFF>::Polynomial
( std::vector<std::string> pvarnames, std::string spol, int pn0 )
 : varnames(pvarnames), np(MAX_POWERS), n0(pn0)
{
    this->checkVariableNames();
    this->fromString( spol ); /* automatically sets size */
}

/* takes the first non alphanumerical character as a delimiter, else space is
 * the standard delimiter to convert a string of variable names to a list */
template<typename T_COEFF>
Polynomial<T_COEFF>::Polynomial
( std::string pvarnames, std::string spol, int pn0 )
 : varnames(), np(MAX_POWERS), n0(pn0)
{
    this->varnames = splitByDelimiter( pvarnames );
    this->checkVariableNames();
    this->fromString( spol ); /* automatically sets size */
}

template<typename T_COEFF>
Polynomial<T_COEFF>::Polynomial( const T_COEFF & src, int pn0 )
 : data(src), varnames(), np(data.getSize()[0]), n0(pn0)
{
    assert( data.isVector() );
    int levels = this->countLevels();
    #if POLYNOMIAL_DEBUG >= 100
        std::cerr << "Counted " << levels << " Levels when constructing from nested Matrix\n";
    #endif
    assert( levels < 'z'-'a' );
    char varname[2] = "a";
    for ( int i = 0; i < levels; ++i ) {
        this->varnames.push_back( std::string(varname) );
        ++varname[0];
    }
}

template<typename T_COEFF>
Polynomial<T_COEFF>::Polynomial
( std::vector<std::string> pvarnames, const T_COEFF & src, int pn0 )
 : data(src), varnames(pvarnames), np( data.getSize()[pn0] ), n0(pn0)
{
    assert( data.isVector() );
}

template<typename T_COEFF>
int Polynomial<T_COEFF>::checkVariableNames( void )
{
    if ( varnames.size() == 0 ) {
        assert( "No varnames given!" && false );
        return -1;
    }

    int wrongvarnames = 0;
    /* Only consisting of alphanum and beginning with letter */
    for ( auto it = varnames.begin(); it != varnames.end(); ++it )
    {
        /* Check if varname is maybe even empty */
        if ( it->length() == 0 ) {
            varnames.erase( it-- );
            continue;
        }
        /* Check first character */
        if ( not isalpha( (*it)[0] ) ) {
            ++wrongvarnames;
            continue;
        }
        /* Check rest */
        for ( auto its = it->begin(); its != it->end(); ++its )
            if ( not isalnum( *its ) ) {
                ++wrongvarnames;
                break;
            }
    }
    if ( wrongvarnames != 0 ) {
        std::cerr << "Variable Names:\n";
        for ( auto it = varnames.begin(); it != varnames.end(); ++it )
            std::cerr << *it << "\n";
        assert( wrongvarnames == 0);
    }
    return wrongvarnames;
}

/******************** Wrapper for recursive template calls ********************/
/* Counts recursion MathMatrix */
template<typename T_COEFF>
inline int Polynomial<T_COEFF>::countLevels( void ) const
{
    return countLevels( this->data );
}

/* Print the polynomial representation of given v, e.g. {0,1,2} -> x+2x**2 */
template<typename T_COEFF>
inline std::string Polynomial<T_COEFF>::toString( void ) const
{
    std::string result = toString( this->data );
    if ( result.length() == 0 )
        return std::string("0");
    else
        return result;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::fromOneVarString( std::string sr )
{
    fromOneVarString( sr, this->data );
    this->np = this->data.getSize()[0];
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::fromString( std::string sr )
{
    fromString( sr, this->data );
    this->np = this->data.getSize()[0];
    return *this;
}

/* convert vector to multiplication operator i.e. matrix */
template<typename T_COEFF>
inline T_COEFF Polynomial<T_COEFF>::toMultiplyer( void ) const
{
    return toMultiplyer( this->data, this->n0 );
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::setSize( int pnp )
{
    this->setSize( this->data, pnp );
    this->np = pnp;
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::setZero( int pn0 )
{
    std::string exportdata = this->toString();
    this->n0 = pn0;
    this->fromString( exportdata );
    return *this;
}

template<typename T_COEFF>
int Polynomial<T_COEFF>::getSize( void ) const
{
    return this->np;
}

template<typename T_COEFF>
int Polynomial<T_COEFF>::getZero( void ) const
{
    return this->n0;
}

template<typename T_COEFF>
Polynomial<T_COEFF> & Polynomial<T_COEFF>::reorderVariableNames
( std::string to )
{
    std::vector<std::string> newnames = this->splitByDelimiter( to );
    /* check if variables are only reordered, not renamed! */
    assert( varnames.size() == newnames.size() );
    for ( auto it = newnames.begin(); it != newnames.end(); ++it )
        assert( find( varnames.begin(), varnames.end(), *it ) != varnames.end() );

    std::string output = this->toString();
    this->varnames.clear();
    this->varnames = newnames;
    this->fromString(output);
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator+=( double b )
{
    (*this) += ::toString(b);
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator-=( double b )
{
    (*this) -= ::toString(b);
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator*=( double b )
{
    (*this) *= ::toString(b);
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator+=( std::string b )
{
    Polynomial<T_COEFF> pol = *this;
    pol.fromString(b);
    this->data += pol.data;
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator-=( std::string b )
{
    Polynomial<T_COEFF> pol = *this;
    pol.fromString(b);
    this->data -= pol.data;
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator*=( std::string b )
{
    Polynomial<T_COEFF> pol = *this;
    pol.fromString(b);
    (*this) *= pol;
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator+=( Polynomial<T_COEFF> b )
{
    this->data += b.data;
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator-=( Polynomial<T_COEFF> b )
{
    this->data -= b.data;
    return *this;
}

template<typename T_COEFF>
inline Polynomial<T_COEFF> & Polynomial<T_COEFF>::operator*=( Polynomial<T_COEFF> b )
{
    /* TODO!!! change np if necessary */
    #if POLYNOMIAL_DEBUG >= 99
        std::cerr << "multiply:\n " << this->data << "\nwith\n " << b.data << "\n";
        std::cerr << "toMultiplyer returned:\n " << toMultiplyer(this->data, this->n0) << "\n";
    #endif
    this->data = toMultiplyer(this->data, this->n0) * b.data;
    #if POLYNOMIAL_DEBUG >= 99
        std::cerr << "multiplied:\n " << this->data << "\n";
    #endif
    return *this;
}

template<typename T_COEFF>
template<typename T_ETYPE>
inline Polynomial<T_COEFF> Polynomial<T_COEFF>::operator+ ( T_ETYPE b ) const
{
    Polynomial<T_COEFF> tmp = *this;
    tmp += b;
    return tmp;
}

template<typename T_COEFF>
template<typename T_ETYPE>
inline Polynomial<T_COEFF> Polynomial<T_COEFF>::operator- ( T_ETYPE b ) const
{
    Polynomial<T_COEFF> tmp = *this;
    tmp -= b;
    return tmp;
}

template<typename T_COEFF>
template<typename T_ETYPE>
inline Polynomial<T_COEFF> Polynomial<T_COEFF>::operator* ( T_ETYPE b ) const
{
    Polynomial<T_COEFF> tmp = *this;
    tmp *= b;
    return tmp;
}

template<typename T_COEFF>
inline bool Polynomial<T_COEFF>::operator==( Polynomial<T_COEFF> b )
{
    return ( this->data == b.data ) and ( this->varnames == b.varnames );
}

template<typename T_COEFF>
void Polynomial<T_COEFF>::operator=( std::string b )
{
    this->fromString( b );
}

template<typename T_COEFF>
void Polynomial<T_COEFF>::operator=( const T_COEFF & src )
{
    int levelsold = this->countLevels();
    this->data = src;
    int levelsnew = this->countLevels();

    assert( data.isVector() );
    this->np = data.getSize()[0];
    /* leave n0 where it was */

    #if POLYNOMIAL_DEBUG >= 99
        std::cerr << "Counted " << levelsnew << " Levels when assigning from nested Matrix\n";
    #endif

    if ( levelsold != levelsnew )
    {
        this->varnames.clear();
        assert( levelsnew < 'z'-'a' );
        char varname = 'a';
        for ( int i = 0; i < levelsnew; ++i )
            this->varnames.push_back( varname );
    }
}

template<typename T_COEFF>
Polynomial<T_COEFF>::operator std::string(void)
{
    return this->toString();
}

/********************************** Integrate *********************************/
template<typename T_COEFF>
auto Polynomial<T_COEFF>::integrate
( std::string variable, std::string from, std::string to ) const
 -> Polynomial<decltype(this->data[0])>
{
    typedef Polynomial<decltype(this->data[0])> ReturnPol;

    /* create integrator by creating new polynom consisting of all variables lower than integration variable -> set up diagonal integration matrix using 1/k*toMultiplyer(fromString("1")) to fill it. Now broadcast that integration matrix into the polynom, which should recurisvely broadcast it the correct level */

    /* Determine which variable we want to integrate and reorder it to the    *
     * front, so that setting up the integration Matrix will be easier        */
    std::vector<std::string> newvarorder = this->varnames;
    auto it = find( newvarorder.begin(), newvarorder.end(), variable );
    assert( it != newvarorder.end() );
    newvarorder.erase( it );
    std::vector<std::string> resultvars = newvarorder;
    newvarorder.insert( newvarorder.begin(), variable );
    Polynomial<T_COEFF> integrand( newvarorder, this->toString() );
    std::cerr << "Created integrand with new variable order\n";

    /* Check if vector is long enough to contain extra power */
    auto maxpower = integrand.data[ integrand.data.getSize()[0] - 1 ];
    if ( maxpower != maxpower*0 ) {
        #if POLYNOMIAL_DEBUG >= 98
            std::cerr << "While integrating Vector needed to be elongated\n";
        #endif
        integrand.setSize( 2*np );
    }
    int newnp = integrand.data.getSize()[0];

    /* Set up minor diagonal integration Matrix, n0 doesn't influence         *
     * diagonal position, only coefficients inside the diagonal               */
    T_COEFF MI(newnp,newnp);
    setSize( MI, newnp, newnp ); // set size recursively
    ReturnPol unit( resultvars, "1" );
    for ( int i = 1; i < newnp; ++i )
        if ( i == n0-1 )
            MI(i,i-1) = unit.toMultiplyer() * INF;
        else
            MI(i,i-1) = unit.toMultiplyer() * ( 1.0 / double( i-n0 ) );
    #if POLYNOMIAL_DEBUG >= 98
        std::cerr << "Set up integration Matrix\n";
    #endif

    /* Create vectors for limits they contain all powers in their vectors     *
     * So if limit a is 'x-1', then vector will be ( 1, x-1, (x-1)**2, ... )  */
    ReturnPol a0( resultvars, from );
    ReturnPol b0( resultvars, to );
    std::vector<ReturnPol> a;
    std::vector<ReturnPol> b;
    a.push_back( unit );
    b.push_back( unit );
    for ( int i = 1; i < newnp; ++i ) {
        #if POLYNOMIAL_DEBUG >= 98
            std::cerr << "Calculate power " << i << "\n";
        #endif
        a.push_back( a0 * a[i-1] );
        b.push_back( b0 * b[i-1] );
    }

    /* Integrate and apply limits */
    #if POLYNOMIAL_DEBUG >= 98
        std::cerr << "Integrate integrand ( size: " << integrand.data.getSize() << " with MI ( size: " << MI.getSize() << " )\n";
        std::cerr << "MI: " << MI << "\n";
        std::cerr << "Integrand: " << integrand.data << "\n";
    #endif
    T_COEFF indefinite = MI * integrand.data;
    #if POLYNOMIAL_DEBUG >= 98
        std::cerr << "Applying limits\n";
    #endif
    assert( (size_t) indefinite.getSize()[0] == a.size() );
    ReturnPol result( resultvars, "0" );
    for ( int i = 0; i < indefinite.getSize()[0]; ++i )
        result += ( b[i] - a[i] ) * ReturnPol( indefinite[i], integrand.getZero() );

    return result;
}

/************************** Recursive template calls **************************/

template<typename T_COEFF>
template<typename T_ELEM>
MathMatrix<MathMatrix<T_ELEM>> & Polynomial<T_COEFF>::setSize
( MathMatrix<MathMatrix<T_ELEM>> & pol, int pnx, int pny )
{
    pol.setSize( pnx, pny );
    for ( int i = 0; i < pol.getSize().product(); ++i )
        setSize( pol[i], pnx, pny );
    /* MathMatrix::setSize sets newly added elements to zero automatically ! */
    return pol;
}

template<typename T_COEFF>
template<typename T_ELEM>
MathMatrix<T_ELEM> & Polynomial<T_COEFF>::setSize
( MathMatrix<T_ELEM> & pol, int pnx, int pny )
{
    pol.setSize( pnx, pny );
    return pol;
}

template<typename T_COEFF>
template<typename T_ELEM>
int Polynomial<T_COEFF>::countLevels
( MathMatrix<T_ELEM> v, int level )
{
    return countLevels( v[0], level+1 );
}

/* Last specialized call for recursion toPolynomial goes here */
template<typename T_COEFF>
template<typename T_ELEM>
int Polynomial<T_COEFF>::countLevels
( T_ELEM v, int level )
{
    return level;
}

/* Last specialized call for recursion toString goes here */
template<typename T_COEFF>
template<typename T_ELEM>
std::string Polynomial<T_COEFF>::toString
( const T_ELEM & v, int level ) const
{
    /* todo: first output ( after '(' ) doesn't have to have a plus sign */
    std::stringstream out;
    out << std::showpos << v << std::noshowpos;
    return out.str();
}

/* Print the polynomial representation of given v, e.g. {0,1,2} -> x+2x**2 */
template<typename T_COEFF>
template<typename T_ELEM>
std::string Polynomial<T_COEFF>::toString
( const MathMatrix<T_ELEM> & v, int level ) const
{
    std::stringstream out;

    /* Count summands to decide whether parentheses are necessary */
    int nsummands = 0;
    for ( int i = 0; i < v.size[0]; ++i )
        if ( v[i] != v[i]*0 )
            ++nsummands;

    bool firstprinted = false;
    for ( int i = 0; i < v.size[0]; ++i )
    {
        /* don't print summands whose coefficients are 0 */
        if ( v[i] == v[i]*0 )
            continue;

        std::string tmp = toString( v[i], level+1 );

        /* count coefficients by counting + and - */
        long int ncoeffs = 1;
        ncoeffs += std::count( tmp.begin(), tmp.end(), '+' );
        ncoeffs += std::count( tmp.begin(), tmp.end(), '-' );
        if ( tmp[0] == '+' or tmp[0] == '-' )
            --ncoeffs;

        /* only print a sign, if needed i.e. if tmp does not begin with one   *
         * and if this is not the very first string to be printed             */
        if ( tmp[0] != '+' and tmp[0] != '-' and firstprinted )
            out << "+";
        firstprinted = true;

        /* if power > 0 and return of toString == "1", then don't print it*/
        int power = i-n0;
        if ( ncoeffs > 1 and not ( nsummands == 1 or power == 0 ) ) {
            out << "(" << tmp << ")";
        } else {
            if ( tmp.compare( "+1" ) == 0 or tmp.compare( "1" ) == 0 ) {
                if ( power == 0 )
                    out << "1";
            } else if ( tmp.compare( "-1" ) == 0 ) {
                if ( power != 0 )
                    out << "-";
                else /* power == 0 */
                    out << "-1";
            } else
                out << tmp;
        }

        /* print variable name: x**(n) if not x**0 */
        if ( power != 0 )
        {
            /* print negative powers as fractions */
            if ( i-n0 < 0 )
                out << "/";
            out << varnames[level];
            /* output power only if necessary */
            if ( i-n0 < 0 )
                out << "**" << -(i-n0);
            else if ( i-n0 != 1 )
                out << "**" << i-n0;
        }
    }

    return out.str();
}

/* Last specialized call for recursion toMultiplyer goes here, T_ELEM is e.g. int */
template<typename T_COEFF>
template<typename T_ELEM>
MathMatrix<T_ELEM> Polynomial<T_COEFF>::toMultiplyer
( const MathMatrix<T_ELEM> & v, int n0 )
{
    int np = v.getVectorDim();
    MathMatrix<T_ELEM> mul(np,np);
    mul = 0;
    for ( int i = 0; i < v.getSize()[0]; ++i )
        mul.setDiagonal( v[i], n0-i );
    return mul;
}

/* convert vector to multiplication operator i.e. matrix */
template<typename T_COEFF>
template<typename T_ELEM>
MathMatrix<MathMatrix<T_ELEM>> Polynomial<T_COEFF>::toMultiplyer
( const MathMatrix<MathMatrix<T_ELEM>> & v, int n0 )
{
    assert( v.isVector() );
    int np = v.getVectorDim();

    /* set up result matrix with correct dimension recursively and 0-elements */
    MathMatrix<MathMatrix<T_ELEM>> mul(np,np);
    MathMatrix<T_ELEM> zero = toMultiplyer( v[0]*0, n0 );
    for ( int i = 0; i < mul.getSize().product(); ++i )
       mul[i] = zero;

    for ( int i = 0; i < v.getSize()[0]; ++i ) {
        /* recursive convert inner vectors to multipliers */
        /* for division set upper diagonal, for multipl. set lower diagonal */
        mul.setDiagonal( toMultiplyer( v[i], n0 ), n0-i );
    }
    return mul;
}

/* converts single coefficients to recursive data format: e.g. x**2 or 3, but *
 * not 3+x                                                                    */

template<typename T_COEFF>
template<typename T_ELEM>
MathMatrix<T_ELEM> & Polynomial<T_COEFF>::fromOneVarString
( std::string str, MathMatrix<T_ELEM> & result, int level ) const
{
    result.setSize(np,1) = 0;

    /* test string for correctness. may not contain ()  */
    if ( str.find('(') != str.npos or str.find(')') != str.npos or
         str.find(' ') != str.npos )
        assert(false);

    /* search for powers: x**n */
    std::string searchstr = varnames[level] + std::string("**");
    size_t foundpos = str.find( searchstr );
    assert( foundpos == 0 or foundpos == str.npos );

    /* if found power x**n */
    if ( foundpos == 0 )
    {
        std::string spower = str.substr( foundpos+varnames[level].length()+2 );
        int power = atoi( spower.c_str() );
        if ( n0 + power >= np )
            result.setSize( n0+power+1, 1 );
        result[n0 + power] = 1;
    }
    /* if found varname without power */
    else if ( str == varnames[level] )
    {
        if ( n0 + 1 >= np )
            result.setSize( n0+1+1, 1 );
        result[n0+1] = (T_ELEM) 1;
    }
    /* if string is number */
    else if ( foundpos == str.npos )
    {
        int dotsfound = 0;
        for ( size_t i = 0; i < str.length(); ++i )
            switch ( str[i] ) {
                case '.':
                    ++dotsfound;
                    break;
                case '+':
                case '-':
                    assert( i==0 );
                    break;
                case '0' ... '9':
                    break;
                default:
                    std::cout << "i:" << i;
                    assert( "maybe invalid varname in string" == false );
            }
        assert( dotsfound <= 1 );
        if ( dotsfound == 0 )
            result[n0] = (T_ELEM) atoi( str.c_str() );
        else
            result[n0] = (T_ELEM) atof( str.c_str() );
    }
    else if ( str.length() == 0 ) {
        /* interpret empty string as 0 */
    } else
        assert(false);

    return result;
}

template<typename T_COEFF>
template<typename T_ELEM>
MathMatrix<MathMatrix<T_ELEM>> & Polynomial<T_COEFF>::fromOneVarString
( std::string str, MathMatrix<MathMatrix<T_ELEM>> & result, int level ) const
{
    MathMatrix<T_ELEM> tmp(np,1);
    setSize( tmp   ,np ) = 0;
    setSize( result,np ) = 0;

    /* test string for correctness. may not contain +-()  */
    if ( str.find('(') != str.npos or str.find(')') != str.npos or
         str.find(' ') != str.npos )
        assert(false);

    /* search for powers: x**n */
    std::string searchstr = varnames[level] + std::string("**");
    size_t foundpos = str.find( searchstr );
    assert( foundpos == 0 or foundpos == str.npos );

    /* if found power x**n */
    int power = 0;
    if ( foundpos == 0 )
    {
        std::string spower = str.substr( foundpos+varnames[level].length()+2 );
        power = atoi( spower.c_str() );
        if ( n0 + power >= np ) {
            setSize( result, n0+power+1 );
            setSize( tmp   , n0+power+1 );
        }
        result[n0 + power] = fromOneVarString( "1", tmp, level+1 );
        if ( result[n0 + power].getSize()[0] > np ) {
            setSize( result, result[n0 + power].getSize()[0] );
            setSize( tmp   , result[n0 + power].getSize()[0] );
        }
    }
    /* if found varname without power */
    else if ( str == varnames[level] )
    {
        power = 1;
        if ( n0 + 1 >= np ) {
            setSize( result, n0+1+1 );
            setSize( tmp   , n0+1+1 );
        }
        result[n0+1] = fromOneVarString( "1", tmp, level+1 );
        if ( result[n0 + 1].getSize()[0] > np ) {
            setSize( result, result[n0 + 1].getSize()[0] );
            setSize( tmp   , result[n0 + 1].getSize()[0] );
        }
    }
    /* if string is number, convert it */
    else if ( foundpos == str.npos )
    {
        power = 0;
        result[n0] = fromOneVarString( str, tmp, level+1 );
        if ( result[n0].getSize()[0] > np ) {
            setSize( result, result[n0].getSize()[0] );
            setSize( tmp   , result[n0].getSize()[0] );
        }
    }
    else
        assert(false);

    return result;
}

template<typename T_COEFF>
template<typename T_ELEM>
MathMatrix<T_ELEM> & Polynomial<T_COEFF>::fromString
( std::string str, MathMatrix<T_ELEM> & result ) const
{
    #if POLYNOMIAL_DEBUG >= 100
        std::cout << "Input String   : " << str << "\n";
        std::cout << "Variable Names : ";
        for ( int i = 0; i < varnames.size(); ++i ) {
            std::cout << varnames[i] << ",";
        }
        std::cout << "\n";
    #endif

    result = fromOneVarString( "0", result );

    /* strip whitespaces */
    {size_t i = 0;
    while ( i < str.length() ) {
        if ( str[i] == ' ' )
            str.erase(i,1);
        else
            ++i;
    }}
    #if POLYNOMIAL_DEBUG >= 100
        std::cout << "Strip Whitespaces: " << str << "\n";
    #endif

    /* strip unnecessary parentheses, e.g. '(a+b)' or '(a)+b' reduce to 'a+b' */
    {assert(    std::count( str.begin(), str.end(), '(' )
            == std::count( str.begin(), str.end(), ')' ) );
    std::stack<size_t> ppos;
    bool pmsignfound = false;
    for ( size_t i = 0; i < str.length(); ++i ) {
        switch ( str[i] ) {
            case '(':
                ppos.push(i);
                break;
            case ')':
                if ( not pmsignfound ) {
                    str.erase( i,1 );
                    str.erase( ppos.top(), 1 );
                    ppos.pop();
                    i -= 2;
                }
                pmsignfound = false;
                break;
            case '+':
            case '-':
                pmsignfound = true;
                break;
        }
    }
    while ( str[0] == '(' and str[ str.length()-1 ] == ')' ) {
        str.erase( str.length()-1, 1 );
        str.erase( 0, 1 );
    }}
    #if POLYNOMIAL_DEBUG >= 100
        std::cout << "Strip unnecessary parentheses: " << str << "\n";
    #endif

    /* analyze highest level only, meaning: 2a*(...)+b-(...) call recursively
     * for (...) */

    MathMatrix<T_ELEM> summand = result*0;
    std::stack<size_t> ppos;
    bool summandinitialized = true;

    /* this is a trick to make this for loop apply the last summand */
    str.push_back('+');

    for ( size_t i = 0; i < str.length(); ++i )
    {
        #if POLYNOMIAL_DEBUG >= 100
            std::cout << "Evaluate " << str.substr(i) << ":\n" << std::flush;
        #endif

        bool tmpchanged = false;
        MathMatrix<T_ELEM> tmp = result*0;

        switch ( str[i] ) {
            /* Add summand to result and reinitialize it */
            case '+':
                #if POLYNOMIAL_DEBUG >= 100
                    std::cerr << "result += summand\n";
                    std::cerr << " result : " << result  << "\n";
                    std::cerr << " summand: " << summand << "\n";
                #endif
                result += summand;
                fromOneVarString( "+0", summand );
                summandinitialized = true;
                break;
            case '-':
                result += summand;
                fromOneVarString( "-1", summand );
                #if POLYNOMIAL_DEBUG >= 100
                    std::cout << "fromOneVarString('-1'): " << summand << "\n";
                #endif
                summandinitialized = false;
                break;

            /* Recursively evaluate things in parentheses */
            case '(':
            {
                /* there is a matching ')', because we checked correctness of string above */
                size_t pos0 = i; /* points to first '(' */
                size_t pos1 = i; /* will point to ')' of same level */
                int parentheses = 1;
                while ( parentheses != 0 )
                    if ( str[++pos1] == ')' )
                        --parentheses;
                /* e.g. "ab(cd)" -> pos0 = 2; pos1 = 5 */
                std::string enclosed = str.substr( pos0+1, pos1-pos0-1 );
                #if POLYNOMIAL_DEBUG >= 100
                    std::cout << "Recursively evaluate: " << enclosed << "\n";
                #endif
                tmp = fromString( enclosed, tmp );
                if ( tmp.getSize()[0] > np ) {
                    setSize( summand, tmp.getSize()[0] );
                    setSize( result , tmp.getSize()[0] );
                }
                #if POLYNOMIAL_DEBUG >= 100
                    std::cout << "Recursively evaluated: " << enclosed << " to " << toString(tmp) << "\n";
                #endif
                i = pos1;
                tmpchanged = true;
                break;
            }

            /* Variables may only begin with a letter */
            case 'A' ... 'Z':
            case 'a' ... 'z':
            {
                bool validvarname = false;
                size_t pos0 = i;
                size_t pos1 = i;
                for ( size_t j = 0; j < varnames.size(); ++j ) {
                    if ( str.substr( pos0 ).find( varnames[j] ) == 0 )
                    {
                        pos1 = i + varnames[j].length();
                        /* check if found thing really is that variable. E.g. *
                         * if "x1" found in "x13", then this isn't correct,   *
                         * because x13 could exist, if not user has to        *
                         * specify this with x1*3                             */
                        if ( ( str[pos1] >= 'A' and str[pos1] <= 'Z' ) or
                             ( str[pos1] >= 'a' and str[pos1] <= 'z' ) or
                               str[pos1] == '*'  or str[pos1] == '+' or
                               str[pos1] == '-'  or str[pos1] == '(' or
                               str[pos1] == ')' )
                        {
                            validvarname = true;
                        } else
                            continue;

                        if ( str[pos1] == '*' and str[pos1+1] == '*' ) {
                            pos1 += 2;
                            while ( str[pos1] >= '0' and str[pos1] <= '9' )
                                ++pos1;
                        }
                        std::string sp = str.substr( pos0, pos1-pos0 );
                        #if POLYNOMIAL_DEBUG >= 100
                            std::cout << "Convert Variable: " << sp << "\n";
                        #endif
                        tmp = fromOneVarString( sp, tmp );
                        if ( tmp.getSize()[0] > np ) {
                            setSize( summand, tmp.getSize()[0] );
                            setSize( result , tmp.getSize()[0] );
                        }
                        #if POLYNOMIAL_DEBUG >= 100
                            std::cout << " to: " << tmp << "\n" << std::flush;
                        #endif
                        tmpchanged = true;
                        i = pos1-1;
                    }
                }
                assert( validvarname == true );
                break;
            }

            /* Evaluate numbers appearing */
            case '0' ... '9':
            {
                tmpchanged = true;
                size_t j = i;
                int dotsfound = 0;
                while ( ( str[j] >= '0' and str[j] <= '9' ) or  str[j] == '.' ) {
                    ++j;
                    if ( str[j] == '.' )
                        ++dotsfound;
                }
                std::string number = str.substr( i, j-i );
                assert( dotsfound <= 1 );
                tmp = fromOneVarString( number, tmp );
                if ( tmp.getSize()[0] > np ) {
                    setSize( summand, tmp.getSize()[0] );
                    setSize( result , tmp.getSize()[0] );
                }
                tmpchanged = true;
                break;
            }
        }

        if ( tmpchanged ) {
            if ( summandinitialized ) {
                #if POLYNOMIAL_DEBUG >= 100
                    std::cout << "summand = tmp\n";
                #endif
                summand = tmp;
            } else {
                #if POLYNOMIAL_DEBUG >= 100
                    std::cout << "summand = summand * tmp\n";
                    auto tmpmul = toMultiplyer( summand );
                    std::cout << "Multiplyer dims: " << tmpmul.getSize() << "\n";
                    #if POLYNOMIAL_DEBUG >= 110
                        std::cout << "raw: " << tmpmul << "\n";
                    #endif
                #endif
                summand = toMultiplyer( summand, n0 ) * tmp;
            }
            summandinitialized = false;
        }
    }

    return result;
}

/*template<typename T_COEFF, typename T_ETYPE>
Polynomial<T_COEFF> operator*( const T_ETYPE a, const Polynomial<T_COEFF> & rhs )
{
    return rhs * a;
}*/

template<typename T_COEFF>
std::ostream& operator<<
( std::ostream& out, const Polynomial<T_COEFF> & p )
{
   out << p.toString();
   return out;
}
