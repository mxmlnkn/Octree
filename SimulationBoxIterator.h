
    class Iterator {
    public:
        /* Only public for convenience, shouldn't be tempered with !!! */
        bool core;
        bool border;
        bool guard;
        
        VecI icell;
        VecI ncells; // with Guard !!!
    private:

        bool inArea( VecI pos ) {
            bool inCore   = ( pos > VecI(2*T_GUARDSIZE-1) and
                              pos - VecI(2*T_GUARDSIZE) < ncells - VecI(4*T_GUARDSIZE) );
            bool inBorder = ( pos > VecI(T_GUARDSIZE-1) and
                              pos < ncells - VecI(T_GUARDSIZE) and not inCore );
            bool inGuard  = ( pos >= VecI(0) and pos < ncells
                              and not inBorder and not inCore );

            #if DEBUG_SIMBOX >= 2
                cout << endl << "[InArea: (" ;
                for ( int i=0; i<T_DIMENSION-1; i++ )
                    cout << pos[i] << ",";
                cout << pos[T_DIMENSION-1] << ") ";
                cout << "is in ";
                if (inCore  ) cout << "CORE ";
                if (inBorder) cout << "BORDER ";
                if (inGuard ) cout << "GUARD ";
                if (!inCore and !inBorder and !inGuard) cout << "None ";
                cout << ", localcells=(" ;
                for ( int i=0; i<T_DIMENSION-1; i++ )
                    cout << ncells[i] << ",";
                cout << ncells[T_DIMENSION-1] << ") ";
                cout << " GuardSize=" << guardsize;
                cout << "]" << endl;
            #endif
            
            bool allor = false;
            if ( guard )
                allor = allor or inGuard;
            if ( border )
                allor = allor or inBorder;
            if ( core )
                allor = allor or inCore;

            return allor;
        }

    public:
        Iterator( const int area, const VecI ncells ) {
            assert( area != 0 ); // would be trivial iterator
            assert( area  < 8 );
            assert( area != CORE+GUARD ); // more difficult to implement,
                                          // and makes no sense
            this->core   = area & CORE;
            this->border = area & BORDER;
            this->guard  = area & GUARD;
            icell = 0;
            this->ncells = ncells;
        }

        ~Iterator( void ) {};

        Iterator( const Iterator & src ) {
            this->core   = src.core  ;
            this->border = src.border;
            this->guard  = src.guard ;
            this->icell  = src.icell ;
            this->ncells = src.ncells;
        }

        Iterator begin( void ) const {
            Iterator tmp( *this );

            int upperleftcorner;
            if (guard)
                upperleftcorner = 0;
            else if (border)
                upperleftcorner = T_GUARDSIZE; // not T_GUARDSIZE-1 !
            else if (core)
                upperleftcorner = 2*T_GUARDSIZE;
            // Above is equivalent to cryptic: upperleftcorner = area / 2;

            tmp.icell = upperleftcorner;
            #if DEBUG_SIMBOX >= 2
                cout << endl << "[Iterator::begin()] (" ;
                for ( int i=0; i<T_DIMENSION-1; i++ )
                    cout << tmp.icell[i] << ",";
                cout << tmp.icell[T_DIMENSION-1] << ") " << endl;
            #endif
            return tmp;
        }

        /* This basically is (ncells-1)-begin() (vector algebra) */
        Iterator end( void ) const {
            Iterator tmp( *this );
            tmp.icell  = VecI(-1);
            return tmp;
        }

        Iterator& operator<<( const int index[T_DIMENSION] ) {
            VecI ind = VecI(index);
            assert( inArea(ind) );
            *this = this->begin();
            this->icell += ind;
            return *this;
        }

        Iterator & operator++( void ) { // only prefix, because it's faster!
            for ( int k=T_DIMENSION-1; k>=0; k-- ) {
                icell[k]++;
                #if DEBUG_SIMBOX >= 2
                    cout << "Current Position: (" ;
                    for ( int i=0; i<T_DIMENSION-1; i++ )
                        cout << icell[i] << ",";
                    cout << icell[T_DIMENSION-1] << ") ";
                    cout << "is" << ((inArea(icell)) ? " " : " not ") << "in ";
                    if (core ) cout << "CORE";
                    if (border) cout << " + BORDER";
                    if (guard ) cout << " + GUARD";
                    cout << endl;
                #endif
                if ( inArea(icell) )
                    break;
                else
                    icell[k] = this->begin().icell[k];
            }
            // If after this everything is at begin(), then overflow occured
            if ( this->begin().icell == this->icell )
                this->icell = VecI(-1);

            return *this;
        }

        bool operator==( const Iterator & it ) {
            bool alland = true;
            alland = alland & ( this->core   == it.core  );
            alland = alland & ( this->border == it.border );
            alland = alland & ( this->guard  == it.guard  );

            alland = alland & ( this->icell  == it.icell  );
            alland = alland & ( this->ncells == it.ncells );

            return alland;
        }

        bool operator!=( const Iterator & it ) {
            return !( (*this) == it );
        }

        void operator=( const Iterator & src ) {
            this->core   = src.core  ;
            this->border = src.border;
            this->guard  = src.guard ;
            this->icell  = src.icell ;
            this->ncells = src.ncells;
        }

        /**********************************************************************
         * If we have 2 slates of 3x4 length and we wanted the index          *
         * [i,j,k] = [2,3,2] , then the matrix lies in the memory like:       *
         *   [oooo|oooo|oooo] [oooo|oooo|oxoo]                                *
         * This means to address that element (above marked x) we need to     *
         * calculate:                                                         *
         *   (i=2)*[ (nj=3) * (nk=4) ] + (j=3)*[ (nk=4) ] + (k=2)*[ 1 ]       *
         * That argument in [] will be named 'prevrange'                      *
         **********************************************************************/
        int getCellIndex( void ) const {
            int index     = 0;
            int prevrange = 1;
            for (int i=T_DIMENSION-1; i>=0; i--) {
                index     += icell[i] * prevrange;
                prevrange *= ncells[i];
            }
            return index;
        }
    };
