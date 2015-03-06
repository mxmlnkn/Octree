#pragma once

/******************** Wrapper for recursive template calls ********************/
template<typename T_POLYNOM>
Polynomial<T_POLYNOM>::Polynomial( std::vector<std::string> pvarnames, std::string spol )
 : varnames(pvarnames), np(5), n0(0)
{
    this->checkVariableNames();
    this->fromString( spol );
}

/* takes the first non alphanumerical character as a delimiter, else space is
 * the standard delimiter to convert a string of variable names to a list */
template<typename T_POLYNOM>
Polynomial<T_POLYNOM>::Polynomial( std::string pvarnames, std::string spol )
 : varnames(), np(5), n0(0)
{
    /* search for delimiter */
    char delimiter = ' ';
    for ( auto it = pvarnames.begin(); it != pvarnames.end(); ++it )
        if ( not isalnum( *it ) )
            delimiter = *it;

    #if POLYNOMIAL_DEBUG >= 100
        std::cerr << "Found delimtier: " << delimiter << "\n";
    #endif

    /* split pvarnames by delimiter */
    int pos1 = 0;
    while (true)
    {
        pvarnames = pvarnames.substr( pos1 );
        pos1 = pvarnames.find( delimiter );

        if ( pos1 != pvarnames.npos ) {
            std::string svar = pvarnames.substr( 0, pos1 );
            if ( svar.length() != 0)
                this->varnames.push_back( svar );
            ++pos1;
        }
        else
        {
            pos1 = pvarnames.length() - 1;
            break;
        }
    }
    if ( pvarnames.length() != 0 )
        this->varnames.push_back( pvarnames );

    this->checkVariableNames();
    this->fromString( spol );
}

template<typename T_POLYNOM>
int Polynomial<T_POLYNOM>::checkVariableNames( void )
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

/* Counts recursion MathMatrix */
template<typename T_POLYNOM>
int Polynomial<T_POLYNOM>::countLevels( void ) const
{
    return countLevels( this->data );
}

/* Print the polynomial representation of given v, e.g. {0,1,2} -> x+2x**2 */
template<typename T_POLYNOM>
std::string Polynomial<T_POLYNOM>::toString( void ) const
{
    std::string result = toString( this->data );
    if ( result.length() == 0 )
        return std::string("0");
    else
        return result;
}

template<typename T_POLYNOM>
Polynomial<T_POLYNOM> & Polynomial<T_POLYNOM>::fromOneVarString( std::string sr )
{
    fromOneVarString( sr, this->data );
    this->np = this->data.getSize()[0];
    return *this;
}

template<typename T_POLYNOM>
Polynomial<T_POLYNOM> & Polynomial<T_POLYNOM>::fromString( std::string sr )
{
    fromString( sr, this->data );
    this->np = this->data.getSize()[0];
    return *this;
}

/* convert vector to multiplication operator i.e. matrix */
template<typename T_POLYNOM>
Polynomial<T_POLYNOM> Polynomial<T_POLYNOM>::toMultiplyer( void ) const
{
    return toMultiplyer( *this );
}

template<typename T_POLYNOM>
Polynomial<T_POLYNOM> & Polynomial<T_POLYNOM>::setSize( int pnp )
{
    this->setSize( this->data, pnp );
    this->np = pnp;
}

/************************** Recursive template calls **************************/

template<typename T_POLYNOM>
template<typename T_ELEM>
MathMatrix<MathMatrix<T_ELEM>> & Polynomial<T_POLYNOM>::setSize
( MathMatrix<MathMatrix<T_ELEM>> & pol, int pnp ) const
{
    pol.setSize( pnp, 1 );
    for ( int i = 0; i < pol.getSize().product(); ++i )
        setSize( pol[i], pnp );
    /* MathMatrix::setSize sets newly added elements to zero automatically ! */
    return pol;
}

template<typename T_POLYNOM>
template<typename T_ELEM>
MathMatrix<T_ELEM> & Polynomial<T_POLYNOM>::setSize
( MathMatrix<T_ELEM> & pol, int pnp ) const
{
    pol.setSize( pnp, 1 );
    return pol;
}


template<typename T_POLYNOM>
template<typename T_ELEM>
int Polynomial<T_POLYNOM>::countLevels
( MathMatrix<T_ELEM> v, int level ) const
{
    return countLevels( v[0], level+1 );
}

/* Last specialized call for recursion toPolynomial goes here */
template<typename T_POLYNOM>
template<typename T_ELEM>
int Polynomial<T_POLYNOM>::countLevels
( T_ELEM v, int level ) const
{
    return level;
}

/* Last specialized call for recursion toString goes here */
template<typename T_POLYNOM>
template<typename T_ELEM>
std::string Polynomial<T_POLYNOM>::toString
( const T_ELEM & v, int level ) const
{
    /* todo: first output ( after '(' ) doesn't have to have a plus sign */
    std::stringstream out;
    out << std::showpos << v << std::noshowpos;
    return out.str();
}

/* Print the polynomial representation of given v, e.g. {0,1,2} -> x+2x**2 */
template<typename T_POLYNOM>
template<typename T_ELEM>
std::string Polynomial<T_POLYNOM>::toString
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
        int ncoeffs = 1;
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
template<typename T_POLYNOM>
template<typename T_ELEM>
MathMatrix<T_ELEM> Polynomial<T_POLYNOM>::toMultiplyer
( const MathMatrix<T_ELEM> & v ) const
{
    int np = v.getVectorDim();
    MathMatrix<T_ELEM> mul(np,np);
    mul = 0;
    for ( int i = 0; i < v.getSize()[0]; ++i )
        mul.setDiagonal( v[i], n0-i );
    return mul;
}

/* convert vector to multiplication operator i.e. matrix */
template<typename T_POLYNOM>
template<typename T_ELEM>
MathMatrix<MathMatrix<T_ELEM>> Polynomial<T_POLYNOM>::toMultiplyer
( const MathMatrix<MathMatrix<T_ELEM>> & v ) const
{
    assert( v.isVector() );
    int np = v.getVectorDim();

    /* set up result matrix with correct dimension recursively and 0-elements */
    MathMatrix<MathMatrix<T_ELEM>> mul(np,np);
    MathMatrix<T_ELEM> zero = toMultiplyer( v[0]*0 );
    for ( int i = 0; i < mul.getSize().product(); ++i )
       mul[i] = zero;

    for ( int i = 0; i < v.getSize()[0]; ++i ) {
        /* recursive convert inner vectors to multipliers */
        /* for division set upper diagonal, for multipl. set lower diagonal */
        mul.setDiagonal( toMultiplyer( v[i] ), n0-i );
    }
    return mul;
}

/* converts single coefficients to recursive data format: e.g. x**2 or 3, but *
 * not 3+x                                                                    */

template<typename T_POLYNOM>
template<typename T_ELEM>
MathMatrix<T_ELEM> & Polynomial<T_POLYNOM>::fromOneVarString
( std::string str, MathMatrix<T_ELEM> & result, int level ) const
{
    result.setSize(np,1) = 0;

    /* test string for correctness. may not contain ()  */
    if ( str.find('(') != str.npos or str.find(')') != str.npos or
         str.find(' ') != str.npos )
        assert(false);

    /* search for powers: x**n */
    std::string searchstr = varnames[level] + std::string("**");
    int foundpos          = str.find( searchstr );
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
        for ( int i = 0; i < str.length(); ++i )
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

template<typename T_POLYNOM>
template<typename T_ELEM>
MathMatrix<MathMatrix<T_ELEM>> & Polynomial<T_POLYNOM>::fromOneVarString
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
    int foundpos          = str.find( searchstr );
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

template<typename T_POLYNOM>
template<typename T_ELEM>
MathMatrix<T_ELEM> & Polynomial<T_POLYNOM>::fromString
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
    {int i = 0;
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
    std::stack<int> ppos;
    bool pmsignfound = false;
    for ( int i = 0; i < str.length(); ++i ) {
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
    std::stack<int> ppos;
    bool summandinitialized = true;

    /* this is a trick to make this for loop apply the last summand */
    str.push_back('+');

    for ( int i = 0; i < str.length(); ++i )
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
                int pos0 = i; /* points to first '(' */
                int pos1 = i; /* will point to ')' of same level */
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
                int pos0 = i;
                int pos1 = i;
                for ( int j = 0; j < varnames.size(); ++j ) {
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
                int j = i;
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
                summand = toMultiplyer( summand ) * tmp;
            }
            summandinitialized = false;
        }
    }

    return result;
}

template<typename T_POLYNOM>
std::ostream& operator<<
( std::ostream& out, const Polynomial<T_POLYNOM> & p )
{
   out << p.toString();
   return out;
}
