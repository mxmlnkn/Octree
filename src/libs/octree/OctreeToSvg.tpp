#pragma once

#ifndef DEBUG_OCTREE_SVG
    #define DEBUG_OCTREE_SVG 0
#endif

namespace Octree {

/******************************** Constructor *********************************/
template<int T_DIM>
OctreeToSvg<T_DIM>::OctreeToSvg
( const Octreetype & ptree, const std::string filename, bool timestamp, int height )
: out(), tree( ptree ), treesrc( &ptree ),
  imagesize( int( double(height) / ptree.size[1] * ptree.size[0] ), height ),
  imageborder( 20,20 ),
  Camera( /* Eye */ Vec<double,3>(0) + Vec<double,3>(0.5, 0.5, 5.0),
          /* looks at */ ptree.center )/* camera is in internal octree units */,
  mvp(4,4), boxesDrawn()
{
	Matrix projection_matrix = CalcProjection( M_PI/2.0, double(imagesize[0]) /
        double(imagesize[1]), 0, 1.5 );
    #if DEBUG_OCTREE_SVG == 0
        tout << "view: " << CalcView(Camera) << "\n";
        tout << "projection: " << projection_matrix << "\n";
        tout << "model matrix: " << CalcModelMatrix( Vec<double,3>(0), 1, 1, 0 ) << "\n";
        tout << "mvp: " << mvp << "\n\n";
    #endif
    mvp = CalcModelMatrix( Vec<double,3>(0), 1, 1, 0 ) * projection_matrix * CalcView(Camera);

    /* Create Timestamp and rank strings for filenames */
    time_t t = time(0);   // get time now
    struct tm * now = localtime( &t );
    std::stringstream fname;
    if (timestamp)
        fname << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
              << now->tm_mday << "_" << now->tm_hour << "-" << now->tm_min
              << "_";
    fname << filename << ".svg";
    this->out.open( fname.str(), std::ofstream::out );

    out << "<?xml version=\"1.0\" standalone=\"no\"?>"              "\n"
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\""       "\n"
        << "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">""\n"
        << "<svg"                                                   "\n"
        << " xmlns:xlink = \"http://www.w3.org/1999/xlink\""        "\n"
        << " xmlns   = \"http://www.w3.org/2000/svg\""              "\n"
        << " width   = \"" << 2*imageborder[0]+imagesize[0] <<    "\"\n"
        << " height  = \"" << 2*imageborder[1]+imagesize[1] <<    "\"\n"
        << " version = \"1.1\""                                     "\n"
        << ">"                                                      "\n"
        << "<defs>"                                                 "\n"
        << "  <style type=\"text/css\">"                            "\n"
        << "   rect {"                                              "\n"
        << "    fill : none;"                                       "\n"
        << "    stroke : white;"                                    "\n"
        << "    stroke-width : " << STROKE_WIDTH <<                 "\n"
        << "   }"                                                   "\n"
        << "   circle {"                                            "\n"
        << "    fill : red;"                                        "\n"
        << "    stroke : none"                                      "\n"
        << "   }"                                                   "\n"
        << "  </style>"                                             "\n"
        << "</defs>"                                                "\n";
    /* Make Canvas Black */
    out << "<rect"                                            "\n"
        << " x     =\"" << 0                             << "\"\n"
        << " y     =\"" << 0                             << "\"\n"
        << " width =\"" << 2*imageborder[0]+imagesize[0] << "\"\n"
        << " height=\"" << 2*imageborder[1]+imagesize[1] << "\"\n"
        << " style =\"fill:black;stroke:none\""               "\n"
        << "/>"                                               "\n";
    out << std::flush;
}


template<>
void OctreeToSvg<2>::PrintGrid(void) {
    /* Update internal copy of the tree, before printing */
    this->tree = *(this->treesrc);

    struct tododata{ int ichild; const Node* node; };
    std::stack<tododata> todo;
    tododata tmp = { /* ichild */ 0, /* node */ this->tree.root };
    todo.push( tmp );
    while( !todo.empty() ) {
        const Node * currentnode = todo.top().node;
        if ( todo.top().ichild == 0 ) {
            Keyvalues toBeStored = { /*ID*/ NboxesDrawn++, /*visible*/ true };
            this->boxesDrawn[ currentnode->center ] = toBeStored;
            /* Draw white stroked rectangle */
            /* We can just scale everything up, because Octree saves all data *
             * scaled down to a box at 0 with size 1 ! :)                     */
            Vec<double,2> upperleft = Vec<double,2>(0.5) + currentnode->center;
            upperleft[1] = 1 - upperleft[1]; // flip along y-axis
            upperleft -= 0.5*currentnode->size;
            upperleft *= imagesize;   // scale image up
            upperleft += imageborder;
            Vec<double,2> rect_size = currentnode->size * Vec<double,2>(imagesize);
            out << "<rect"                            "\n"
                << " id    =\"" << toBeStored.id << "\"\n"
                << " x     =\"" << upperleft[0] <<  "\"\n"
                << " y     =\"" << upperleft[1] <<  "\"\n"
                << " width =\"" << rect_size[0] <<  "\"\n"
                << " height=\"" << rect_size[1] <<  "\"\n"
                << "/>"                               "\n";
        }
        if ( currentnode->IsLeaf() )
            todo.pop();
        else {
            if ( currentnode->getChildPtr( todo.top().ichild ) != NULL ) {
                tmp.node   = currentnode->getChildPtr( todo.top().ichild++ );
                tmp.ichild = 0;
                todo.push( tmp );
            } else {
                todo.pop();
            }
        }
    }
    out << std::flush;
}

template<>
void OctreeToSvg<3>::PrintGrid(void) {
    /* Update internal copy of the tree, before printing */
    this->tree = *(this->treesrc);

    struct tododata{ int ichild; const Node* node; };
    std::stack<tododata> todo;
    tododata tmp = { /* ichild */ 0, /* node */ this->tree.root };
    todo.push( tmp );
    while( !todo.empty() ) {
        const Node * currentnode = todo.top().node;
        if ( todo.top().ichild == 0 ) {
            Keyvalues toBeStored = { /*ID*/ NboxesDrawn++, /*visible*/ true };
            this->boxesDrawn[ currentnode->center ] = toBeStored;
            /* Draw white stroked rectangle */
            /* We can just scale everything up, because Octree saves all data *
             * scaled down to a box at 0 with size 1 ! :)                     */
            Vec<double,2> upperleft = Vec<double,2>(0.5) + Vec<double,2>( currentnode->center[0], currentnode->center[1] );
            upperleft[1] = 1 - upperleft[1]; // flip along y-axis
            upperleft -= 0.5*Vec<double,2>( currentnode->size[0], currentnode->size[1] );
            upperleft *= imagesize;   // scale image up
            upperleft += imageborder;
            Vec<double,2> rect_size = Vec<double,2>( currentnode->size[0], currentnode->size[1]) * Vec<double,2>(imagesize);
            out << "<rect"                            "\n"
                << " id    =\"" << toBeStored.id << "\"\n"
                << " x     =\"" << upperleft[0] <<  "\"\n"
                << " y     =\"" << upperleft[1] <<  "\"\n"
                << " width =\"" << rect_size[0] <<  "\"\n"
                << " height=\"" << rect_size[1] <<  "\"\n"
                << "/>"                               "\n";
        }
        if ( currentnode->IsLeaf() )
            todo.pop();
        else {
            if ( currentnode->getChildPtr( todo.top().ichild ) != NULL ) {
                tmp.node   = currentnode->getChildPtr( todo.top().ichild++ );
                tmp.ichild = 0;
                todo.push( tmp );
            } else {
                todo.pop();
            }
        }
    }
    out << std::flush;
}

template<int T_DIM>
void OctreeToSvg<T_DIM>::AnimateUpdated( const Octreetype & newtree )
{
    /*std::cerr << "-boxesDrawn-\n";
    typename std::map<VecD,Keyvalues>::iterator it = boxesDrawn.begin();
    while( it != boxesDrawn.end() ) {
        std::cerr << boxesDrawnIt->first << " => ("
                  << boxesDrawnIt->second.id << ","
                  << boxesDrawnIt->second.visible << ")\n";
        ++it;
    }*/

    typedef struct{ int ichild; const Node* node; } tododata;
    std::stack<tododata> todo;
    tododata tmp = { 0, newtree.root };
    todo.push( tmp );
    while( !todo.empty() ) {
        const Node * currentnode = todo.top().node;

        if ( todo.top().ichild == 0 )
        {
            Vec<double,2> boxcenter = currentnode->center;
            typename VecDMap::iterator boxesDrawnIt = boxesDrawn.find( boxcenter );
            #if DEBUG_OCTREE_SVG >= 10
            std::cerr << boxcenter << " already drawn? "
                      << !(boxesDrawnIt == boxesDrawn.end() )
                      << " boxesDrawnSize: " << boxesDrawn.size() << std::endl;
            if ( boxesDrawnIt != boxesDrawn.end() )
                std::cerr << " Found: " << boxesDrawnIt->first << " => ("
                          << boxesDrawnIt->second.id << ","
                          << boxesDrawnIt->second.visible << ")\n";
            #endif
            /* A new octant was created. We need to draw that first */
            if ( boxesDrawnIt == boxesDrawn.end() ) {
/* By setting visibility to false we ensure that in the If-case thereafter    *
 * the stroke will be set to white per <set/>. We do not set it right now to  *
 * white, because then it would be visible from the start !                   */
                Keyvalues toBeStored = { NboxesDrawn++, /*visible*/ false };
                this->boxesDrawn[ currentnode->center ] = toBeStored;
                boxesDrawnIt = boxesDrawn.find( boxcenter );

                #if DEBUG_OCTREE_SVG >= 10
                std::cerr << "Box around " << boxcenter << " was created !!! "
                          << "ID : " << toBeStored.id << "("
                          << boxesDrawnIt->second.id << ")\n";
                #endif

                Vec<double,2> upperleft = Vec<double,2>(0.5) + currentnode->center;
                upperleft[1] = 1 - upperleft[1]; // flip along y-axis
                upperleft -= 0.5*currentnode->size;
                upperleft *= imagesize;   // scale image up
                upperleft += imageborder;
                Vec<double,2> rect_size = currentnode->size * imagesize;

                out << "<rect"                            "\n"
                    << " id    =\"" << toBeStored.id << "\"\n"
                    << " x     =\"" << upperleft[0] <<  "\"\n"
                    << " y     =\"" << upperleft[1] <<  "\"\n"
                    << " width =\"" << rect_size[0] <<  "\"\n"
                    << " height=\"" << rect_size[1] <<  "\"\n"
                    << " style =\"stroke:none\""          "\n"
                    << "/>"                               "\n";
            }
            /* Some Box which was drawn, but is not visible, reappeard */
            if (   boxesDrawnIt != boxesDrawn.end() and
                 ! boxesDrawnIt->second.visible )
            {
                #if DEBUG_OCTREE_SVG >= 10
                std::cerr << "Box around " << boxcenter << " reappeared !!! "
                          << "ID : " << boxesDrawnIt->second.id << "\n";
                #endif

                out << "<set"                                            "\n"
                    << " xlink:href=\"#" << boxesDrawnIt->second.id << "\"\n"
                    << " attributeName=\"stroke\""                       "\n"
                    << " begin=\"" << DUR*(currentTime+0.5) <<        "s\"\n"
                    << " to   =\"white\""                                "\n"
                    << "/>"                                              "\n";
                boxesDrawnIt->second.visible = true;
            }
        }

/* If the Current Node is not a leaf, then increment child-index and push the *
 * next child to be processed                                                 */
        else
        {
            if ( currentnode->getChildPtr( todo.top().ichild ) != NULL ) {
                tmp.node   = currentnode->getChildPtr( todo.top().ichild++ );
                tmp.ichild = 0;
                todo.push(tmp);
            } else {
                todo.pop();
            }
        }
    }

/* Now iterate trough all the Boxes which were already drawn at least once    *
 * check if they are visible, if so, then check if they still exist in        *
 * newtree. If not they will be set to not visible in the map and with <set/> */

    for ( typename VecDMap::iterator boxesDrawnIt = boxesDrawn.begin();
          boxesDrawnIt != boxesDrawn.end(); ++boxesDrawnIt )
    {
        Vec<double,2> boxcenter = boxesDrawnIt->first;
        if ( boxesDrawnIt->second.visible and
             newtree.GetNodePtr( boxcenter ) == NULL )
        {
            #if DEBUG_OCTREE_SVG >= 10
            std::cerr << "Box around " << boxcenter << " vanished!!! "
                      << "ID: " << boxesDrawnIt->second.id << "\n";
            #endif
            out << "<set"                                            "\n"
                << " xlink:href=\"#" << boxesDrawnIt->second.id << "\"\n"
                << " attributeName=\"stroke\""                       "\n"
                << " begin=\"" << DUR*(currentTime+0.5) << "s\""     "\n"
                << " to   =\"none\""                                 "\n"
                << "/>"                                              "\n";
            boxesDrawnIt->second.visible = false;
        }
    }

    currentTime += 1;
    this->tree = newtree;
    out << std::flush;
}

template<int T_DIM>
Vec<double,2> OctreeToSvg<T_DIM>::convertToImageCoordinates
( Vec<double,T_DIM> pos )
{
    pos   += 0.5; // shift internal coordinate from [-0.5,0.5) to [0,1) in every(!) axis

    Matrix pmvp(2,3);
    pmvp(0,0) = 1; pmvp(0,1) = 0; pmvp(0,2) = -1.0/sqrt(2);// pmvp(0,3) = 0;
    pmvp(1,0) = 0; pmvp(1,1) = 1; pmvp(1,2) = -1.0/sqrt(2);// pmvp(1,3) = 0;

    /*tout << "pos=" << pos << " -> newpos="
         <<  pmvp * Matrix(pos) << ")";*/

    Vec<double,2> newpos;
    if (T_DIM == 3) {
        Vec<double,2> tmppos = static_cast<Vec<double,2>>( pmvp*Matrix(pos) );
        /*for (int i=0; i<4; i++)
            tmppos /= tmppos[3];*/
        /*newpos[0] = tmppos[0] - tmppos[2] / sqrt(2);
        newpos[1] = tmppos[1] - tmppos[2] / sqrt(2);*/
        /* z in [0,1) -> x,y in [-1/sqrt(2),1) => need to scale and shift picture coords */
        newpos = (tmppos + 1.0/sqrt(2)) / (1 + 1.0/sqrt(2));
        assert( newpos[0] >= 0.0 and newpos[0] <= 1.0 );
        assert( newpos[1] >= 0.0 and newpos[1] <= 1.0 );
    } else if (T_DIM == 2)
        newpos = Vec<double,2>( pos[0], pos[1] );
    else
        assert( T_DIM == 2 or T_DIM == 3 );

    newpos[1] = 1 - newpos[1];  // flip along y-axis
    newpos *= imagesize;        // scale image up
    newpos += imageborder;
    //tout << " -> " << newpos << "\n";
    return newpos;
}

} // namespace Octree


