#pragma once

#define DEBUG_OCTREE 9



namespace Octree {

/* Some methods we want do forbid from use!                                   */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::Node( void ) {
    assert(false);
}
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::Node(const Node &) {
    assert(false);
}
template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::Node& Node<T_DTYPE,T_DIM>::operator=(const Node &) {
    assert(false);
}

/* Destructor */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::~Node( void ) {
    for (int i=0; i<dim; i++)
        if ( this->children[i] != NULL )
            delete this->children[i];
}

/* Constructor */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::Node(Node * const parent, VecD const cent, VecD const size)
: parent( parent ), center( cent ), size( size ) {
    for (int i=0; i<dim; i++)
        this->children[i] = NULL;
}

template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::IsLeaf( void ) const {
/* We assume the correctness of this class, i.e. all children are NULL or all *
 * children are not NULL                                                      */
    return children[0] == NULL;
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::VecI Node<T_DTYPE,T_DIM>::ConvertNumberToDirection
( const int number ) const {
    VecI direction;
    int tmp = number;
    for (int i=T_DIM-1; i>=0; --i) {
        direction[i] = ( (tmp & 1) == 0 ) ? -1 : 1;
        tmp = tmp >> 1;
    }
#if DEBUG_OCTREE >= 11
    std::cerr << "number: " << number << " => " << direction << std::endl;
#endif
    //assert( ConvertDirectionToNumber(direction) == number );
    return direction;
}

template<typename T_DTYPE, int T_DIM>
int Node<T_DTYPE,T_DIM>::ConvertDirectionToNumber( const VecI direction ) const {
    int tmp = 0;
    for (int i=0; i<T_DIM; ++i) {
        tmp = ( tmp << 1) | (1+direction[i])/2;
    }
#if DEBUG_OCTREE >= 11
    std::cerr << "direction: " << direction << " => " << tmp << std::endl;
#endif
    return tmp;
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::Node& Node<T_DTYPE,T_DIM>::FindLeafContainingPos
( const VecD & pos ) {
/* Prone to rounding errors :(, except if worldsize is e.g. 1 and center is   *
 * 0, because in that case all sizes and new centers will be in the form of   *
 * 1/2^N which can be represented exactly with floating points!               *
 *   => Use this internally and overlay it with user-sizes !                  */   bool insideNode = (pos >= (this->center - 0.5*this->size)) and
                      (pos <  (this->center + 0.5*this->size));
    if (!insideNode) {
        std::cout << "pos: " << pos << " center: " << this->center << " size: "
                  << this->size << std::endl;
        assert( insideNode );
    }
    Node * tmp = this;
#if DEBUG_OCTREE >= 10
        std::cerr << "Try to find correct leaf for " << pos << std::endl;
#endif
    while( ! tmp->IsLeaf() ) {
/* E.g. center is (0,0) then pos =(0.3,-0.1) will with Greater() result in    *
 * (1,0) which in turn will be converted to (1,-1).                           */
        VecD direction = 2*VecI( pos.GreaterThan( tmp->center )) - 1;
        tmp = tmp->children[ ConvertDirectionToNumber( direction ) ];
#if DEBUG_OCTREE >= 10
        std::cerr << "Comparison with center " << tmp->center << " resulted in "
                  << pos.GreaterThan( tmp->center ) << " -> "
                  << 2*VecI( pos.GreaterThan( tmp->center )) - 1 << std::endl;
#endif
    }
    return *tmp;
}

template<typename T_DTYPE, int T_DIM>
void Node<T_DTYPE,T_DIM>::GrowUp( void ) {
#if DEBUG_OCTREE >= 10
    std::cerr << "Leaf at " << this->center << " growing Up!\n";
#endif
    for (int i=0; i<compileTime::pow(2,T_DIM); ++i) {
        VecD direction = VecD( ConvertNumberToDirection( i ) );
        this->children[i] = new class Node( this, this->center +
                            0.25*direction*size, 0.5*size );
    }
/* Migrate data to child-leaf nodes. If it is too much data they will recursi-*
 * vely grow up to, which is not nice, but it shouldn't happen too often,     *
 * because we only added one element. For above to happen alle elements       *
 * need to be in the same next-gen octant                                     */
    typename Datalist::iterator it = data.begin();
    while( it != data.end() ) {
        FindLeafContainingPos( it->pos ).InsertData( it->pos, it->object );
        this->data.erase( it++ );
    }
    assert( this->data.empty() );
#if DEBUG_OCTREE >= 10
    std::cerr << "All Grown Up!\n";
#endif
}

template<typename T_DTYPE, int T_DIM>
void Node<T_DTYPE,T_DIM>::InsertData( const VecD pos, T_DTYPE * const object ) {
    assert( this->IsLeaf() );
    Datapoint packed = { pos, object };
/* Push back no matter of the final size. It will be split up when splitting  *
 * the octant to 8 new octants anyway                                         */
    this->data.push_back( packed );
    if ( data.size() > this->maxdata )
        this->GrowUp();
}


template<typename T_DTYPE, int T_DIM>
void Node<T_DTYPE,T_DIM>::RemoveData( const VecD pos, T_DTYPE * const object ) {

}

template<typename T_DTYPE, int T_DIM>
void Node<T_DTYPE,T_DIM>::MoveData( const VecD pos, T_DTYPE * const object,
const VecD newpos ) {

}

/************************ Wrapper / Node manager class ************************/


template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::Octree( void ) {
    assert(false);
}

template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::Octree( const Octree & ) {
    assert(false);
}

template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::~Octree( void ) {
    delete this->root;
}

template<typename T_DTYPE, int T_DIM>
typename Octree<T_DTYPE,T_DIM>::Octree& Octree<T_DTYPE,T_DIM>::operator=(const Octree &)
{
    assert(false);
}

template<typename T_DTYPE, int T_DIM>
Octree<T_DTYPE,T_DIM>::Octree( const VecD center, const VecD size )
: center( center ), size( size ) {
/* root Node has no parent, therefore it's parent is set to NULL !            */
    this->root = new class Node<T_DTYPE,T_DIM>( NULL, VecD(0), VecD(1) );
}

template<typename T_DTYPE, int T_DIM>
void Octree<T_DTYPE,T_DIM>::InsertData( const VecD pos, T_DTYPE * const object ) {
    this->root->FindLeafContainingPos( pos / this->size ).
                InsertData( pos / this->size, object );
}

template<typename T_DTYPE, int T_DIM>
void Octree<T_DTYPE,T_DIM>::PrintToSVG( const std::string filename ) {
    /* Create Timestamp and rank strings for filenames */
    time_t t = time(0);   // get time now
    struct tm * now = localtime( &t );
    std::stringstream timestamp;
    timestamp << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
              << now->tm_mday << "_" << now->tm_hour << "-"
              << now->tm_min  << "_";
    std::string timesr = timestamp.str();
    std::ofstream out;
    out.open( timesr + filename + std::string(".svg"),
                     std::ofstream::out | std::ofstream::app );
    double borderx = 20;
    double bordery = 20;
    double height  = 800;
    double width   = height / this->size[1] * this->size[0];
    VecD imagesize  ; imagesize  [0] = width  ; imagesize  [1] = height ;
    VecD imageborder; imageborder[0] = borderx; imageborder[1] = bordery;
    out << "<?xml version=\"1.0\" standalone=\"no\"?>              \n"
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"       \n"
        << "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
        << "<svg                                                   \n"
        << "   xmlns   = \"http://www.w3.org/2000/svg\">           \n"
        << "   width   = \"" << 2*borderx+width  << "px\"          \n"
        << "   height  = \"" << 2*bordery+height << "px\"          \n"
        << "   version = \"1.1\"                                   \n"
        << "  >                                                    \n";
    /* Make Canvas Black */
    out << "  <rect                                     \n"
        << "     x      = \"" << 0                << "\"\n"
        << "     y      = \"" << 0                << "\"\n"
        << "     width  = \"" << 2*borderx+width  << "\"\n"
        << "     height = \"" << 2*bordery+height << "\"\n"
        << "     fill   = \"black\"                     \n"
        << "     stroke = \"none\"                      \n"
        << "  />                                        \n";

    typedef class Node<T_DTYPE,T_DIM> Node;
    typedef struct{ int ichild; Node* node; } tododata;
    std::stack<tododata> todo;
    tododata tmp = { 0, this->root };
    todo.push( tmp );
    while( !todo.empty() ) {
        Node & current = *(todo.top().node);
        if ( todo.top().ichild == 0 ) {
            /* Draw white stroked rectangle */
            /* We can just scale everything up, because Octree saves all data *
             * scaled down to a box at 0 with size 1 ! :)                     */
            VecD upperleft = VecD(0.5) + current.center;
            upperleft[1] = 1 - upperleft[1]; // flip along y-axis
            upperleft -= 0.5*current.size;
            upperleft *= imagesize;   // scale image up
            upperleft += imageborder;
            VecD rect_size = current.size * imagesize;
            out << "  <rect                                         \n"
                << "     x      = \"" << upperleft[0] << "px" << "\"\n"
                << "     y      = \"" << upperleft[1] << "px" << "\"\n"
                << "     width  = \"" << rect_size[0] << "px" << "\"\n"
                << "     height = \"" << rect_size[1] << "px" << "\"\n"
                << "     fill   = \"none\"                          \n"
                << "     stroke = \"white\"                         \n"
                << "     stroke-width = \"3px\"                     \n"
                << "  />                                            \n";
        }
        if ( current.IsLeaf() ) {
            typename Node::Datalist::iterator  it = current.data.begin();
            while ( it != current.data.end() ) {
                /* Draw object into the world! */
                /* out << it->pos * imagesize << *(it->object) */
                VecD center( VecD(0.5) + it->pos );
                center[1] = 1 - center[1];
                center *= imagesize;   // scale image up
                center += imageborder;
                out << "  <circle                                    \n"
                    << "     cx     = \"" << center[0] << "px" << "\"\n"
                    << "     cy     = \"" << center[1] << "px" << "\"\n"
                    << "     r      = \"3px\"                        \n"
                    << "     fill   = \"red\"                        \n"
                    << "     stroke = \"none\"                       \n"
                    << "  />                                         \n";
                ++it;
            }
            todo.pop();
        } else {
            if ( todo.top().ichild < current.nchildren ) {
                tmp.node   = current.children[ todo.top().ichild++ ];
                tmp.ichild = 0;
                todo.push(tmp);
            } else {
                todo.pop();
            }
        }
    }
    out << "</svg>\n";
    out.close();
}


} // namespace Octree




/******************************** Debug output ********************************/

template<typename T_DTYPE, int T_DIM>
std::ostream& operator<<( std::ostream& out, const Octree::Octree<T_DTYPE,T_DIM>& tree ) {
    typedef class Octree::Node<T_DTYPE,T_DIM> Node;
    typedef struct{ int ichild; Node* node; } tododata;
    std::stack<tododata> todo;
    tododata tmp = { 0, tree.root };
    todo.push( tmp );
    while( !todo.empty() ) {
        Node & current = *(todo.top().node);
        if ( current.IsLeaf() ) {
            for ( int i=0; i<todo.size()-1; ++i )
                out << "  ";
            out << "Leaf at "    << tree.size * current.center
                << " with size " << tree.size * current.size << std::endl;
            typename Node::Datalist::iterator  it = current.data.begin();
            while ( it != current.data.end() ) {
                for ( int i=0; i<todo.size(); ++i )
                    out << "  ";
                out << it->pos << " -> " << *(it->object) << std::endl;
                ++it;
            }
            todo.pop();
        } else {
            if ( todo.top().ichild == 0 ) {
                for ( int i=0; i<todo.size()-1; ++i )
                    out << "  ";
                out << "Parent at "  << tree.size * current.center
                    << " with size " << tree.size * current.size << std::endl;
            }
            if ( todo.top().ichild < current.nchildren ) {
                tmp.node   = current.children[ todo.top().ichild++ ];
                tmp.ichild = 0;
                todo.push(tmp);
            } else {
                todo.pop();
            }
        }
    }
    return out;
}

