#pragma once


#define DEBUG_OCTREE_SVG 0


namespace Octree {




template<typename T_DTYPE, int T_DIM>
OctreeToSvg<T_DTYPE,T_DIM>::OctreeToSvg (
  const Octree<T_DTYPE,T_DIM> & tree,
  const std::string filename )
: tree( tree ), borderx( 20 ), bordery( 20 ), height( 800 ),
  width( height / tree.size[1] * tree.size[0] )
{
    /* Create Timestamp and rank strings for filenames */
    time_t t = time(0);   // get time now
    struct tm * now = localtime( &t );
    std::stringstream timestamp;
    timestamp << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
              << now->tm_mday << "_" << now->tm_hour << "-"
              << now->tm_min  << "_";
    std::string timesr = timestamp.str();
    this->out.open( timesr + filename + std::string(".svg"),
                     std::ofstream::out );
    this->imagesize  [0] = width  ; this->imagesize  [1] = height ;
    this->imageborder[0] = borderx; this->imageborder[1] = bordery;
    out << "<?xml version=\"1.0\" standalone=\"no\"?>"              "\n"
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\""       "\n"
        << "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">""\n"
        << "<svg"                                                   "\n"
        << " xmlns:xlink = \"http://www.w3.org/1999/xlink\""        "\n"
        << " xmlns   = \"http://www.w3.org/2000/svg\">"             "\n"
        << " width   = \"" << 2*borderx+width  << "px\""            "\n"
        << " height  = \"" << 2*bordery+height << "px\""            "\n"
        << " version = \"1.1\""                                     "\n"
        << ">"                                                      "\n"
        << "<defs>"                                                 "\n"
        << "  <style type=\"text/css\">"                            "\n"
        << "   rect {"                                              "\n"
        << "    fill : none;"                                       "\n"
        << "    stroke : white;"                                    "\n"
        << "    stroke-width : " << STROKE_WIDTH << "px"            "\n"
        << "   }"                                                   "\n"
        << "   circle {"                                            "\n"
        << "    fill : red;"                                        "\n"
        << "    stroke : none"                                      "\n"
        << "   }"                                                   "\n"
        << "  </style>"                                             "\n"
        << "</defs>"                                                "\n";
    /* Make Canvas Black */
    out << "<rect"                               "\n"
        << " x     =\"" << 0                << "\"\n"
        << " y     =\"" << 0                << "\"\n"
        << " width =\"" << 2*borderx+width  << "\"\n"
        << " height=\"" << 2*bordery+height << "\"\n"
        << " style =\"fill:black;stroke:none\"    \n"
        << "/>"                                  "\n";
}


template<typename T_DTYPE, int T_DIM>
void OctreeToSvg<T_DTYPE,T_DIM>::PrintGrid(void) {
    typedef struct{ int ichild; const Node* node; } tododata;
    std::stack<tododata> todo;
    tododata tmp = { /* ichild */ 0, /* node */ this->tree.root };
    todo.push( tmp );
    while( !todo.empty() ) {
        const Node * currentnode = todo.top().node;
        if ( todo.top().ichild == 0 ) {
            Keyvalues toBeStored = { /*ID*/ NBoxesDrawn++, /*visible*/ true };
            this->BoxesDrawn[ currentnode->center ] = toBeStored;
            /* Draw white stroked rectangle */
            /* We can just scale everything up, because Octree saves all data *
             * scaled down to a box at 0 with size 1 ! :)                     */
            VecD upperleft = VecD(0.5) + currentnode->center;
            upperleft[1] = 1 - upperleft[1]; // flip along y-axis
            upperleft -= 0.5*currentnode->size;
            upperleft *= imagesize;   // scale image up
            upperleft += imageborder;
            VecD rect_size = currentnode->size * imagesize;
            out << "<rect"                                   "\n"
                << " id    =\"" << toBeStored.id        << "\"\n"
                << " x     =\"" << upperleft[0] << "px" << "\"\n"
                << " y     =\"" << upperleft[1] << "px" << "\"\n"
                << " width =\"" << rect_size[0] << "px" << "\"\n"
                << " height=\"" << rect_size[1] << "px" << "\"\n"
                << "/>"                                      "\n";
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
}


template<typename T_DTYPE, int T_DIM>
void OctreeToSvg<T_DTYPE,T_DIM>::PrintPositions(void) {
    typedef struct{ int ichild; const Node* node; } tododata;
    std::stack<tododata> todo;
    tododata tmp = { 0, this->tree.root };
    todo.push( tmp );
    while( !todo.empty() ) {
        const Node * current = todo.top().node;
        if ( current->IsLeaf() ) {
            int i = 0;
            while ( current->getDataPtr(i) != NULL ) {
                /* Draw object into the world! */
                /* out << it->pos * imagesize << *(it->object) */
                VecD center( VecD(0.5) + current->getDataPtr(i)->pos );
                center[1] = 1 - center[1];
                center *= imagesize;   // scale image up
                center += imageborder;
                size_t id = reinterpret_cast<size_t>( current->getDataPtr(i)->object );
                out << "<circle"                          "\n"
                    << " id=\"" << id                << "\"\n"
                    << " cx=\"" << center[0] << "px" << "\"\n"
                    << " cy=\"" << center[1] << "px" << "\"\n"
                    << " r =\"3px\""                      "\n"
                    << "/>"                               "\n";
                ++i;
            }
            todo.pop();
        } else {
            if ( current->getChildPtr( todo.top().ichild ) != NULL ) {
                tmp.node   = current->getChildPtr( todo.top().ichild++ );
                tmp.ichild = 0;
                todo.push(tmp);
            } else {
                todo.pop();
            }
        }
    }
}

template<typename T_DTYPE, int T_DIM>
void OctreeToSvg<T_DTYPE,T_DIM>::AnimateUpdated( const Octreetype & newtree )
{
    /*std::cerr << "-BoxesDrawn-\n";
    typename std::map<VecD,Keyvalues>::iterator it = BoxesDrawn.begin();
    while( it != BoxesDrawn.end() ) {
        std::cerr << BoxesDrawnIt->first << " => ("
                  << BoxesDrawnIt->second.id << ","
                  << BoxesDrawnIt->second.visible << ")\n";
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
            VecD boxcenter = currentnode->center;
            BoxesDrawnIt   = BoxesDrawn.find( boxcenter );
            #if DEBUG_OCTREE_SVG >= 10
            std::cerr << boxcenter << " already drawn? "
                      << !(BoxesDrawnIt == BoxesDrawn.end() )
                      << " BoxesDrawnSize: " << BoxesDrawn.size() << std::endl;
            if ( BoxesDrawnIt != BoxesDrawn.end() )
                std::cerr << " Found: " << BoxesDrawnIt->first << " => ("
                          << BoxesDrawnIt->second.id << ","
                          << BoxesDrawnIt->second.visible << ")\n";
            #endif
            /* A new octant was created. We need to draw that first */
            if ( BoxesDrawnIt == BoxesDrawn.end() ) {
/* By setting visibility to false we ensure that in the If-case thereafter    *
 * the stroke will be set to white per <set/>. We do not set it right now to  *
 * white, because then it would be visible from the start !                   */
                Keyvalues toBeStored = { NBoxesDrawn++, /*visible*/ false };
                this->BoxesDrawn[ currentnode->center ] = toBeStored;
                BoxesDrawnIt = BoxesDrawn.find( boxcenter );

                #if DEBUG_OCTREE_SVG >= 10
                std::cerr << "Box around " << boxcenter << " was created !!! "
                          << "ID : " << toBeStored.id << "("
                          << BoxesDrawnIt->second.id << ")\n";
                #endif

                VecD upperleft = VecD(0.5) + currentnode->center;
                upperleft[1] = 1 - upperleft[1]; // flip along y-axis
                upperleft -= 0.5*currentnode->size;
                upperleft *= imagesize;   // scale image up
                upperleft += imageborder;
                VecD rect_size = currentnode->size * imagesize;

                out << "<rect"                                   "\n"
                    << " id    =\"" << toBeStored.id        << "\"\n"
                    << " x     =\"" << upperleft[0] << "px" << "\"\n"
                    << " y     =\"" << upperleft[1] << "px" << "\"\n"
                    << " width =\"" << rect_size[0] << "px" << "\"\n"
                    << " height=\"" << rect_size[1] << "px" << "\"\n"
                    << " style =\"stroke:none\""                 "\n"
                    << "/>"                                      "\n";
            }
            /* Some Box which was drawn, but is not visible, reappeard */
            if (   BoxesDrawnIt != BoxesDrawn.end() and
                 ! BoxesDrawnIt->second.visible )
            {
                #if DEBUG_OCTREE_SVG >= 10
                std::cerr << "Box around " << boxcenter << " reappeared !!! "
                          << "ID : " << BoxesDrawnIt->second.id << "\n";
                #endif

                out << "<set"                                            "\n"
                    << " xlink:href=\"#" << BoxesDrawnIt->second.id << "\"\n"
                    << " attributeName=\"stroke\""                       "\n"
                    << " begin=\"" << DUR*(currentTime+0.5) <<        "s\"\n"
                    << " to   =\"white\""                                "\n"
                    << "/>"                                              "\n";
                BoxesDrawnIt->second.visible = true;
            }
        }

        if ( currentnode->IsLeaf() ) {
            int i = 0;
            while ( currentnode->getDataPtr(i) != NULL ) {
                T_DTYPE * datum = currentnode->getDataPtr(i)->object;
                VecD newpos = currentnode->getDataPtr(i)->pos;
                VecD oldpos = this->tree.FindData( datum ) / this->tree.size;
                if ( oldpos != newpos ) {
                    #if DEBUG_OCTREE_SVG >= 9
                    std::cerr << "Data " << datum << " = ";
                    if (datum!=NULL)
                        std::cerr << *datum;
                    else
                        std::cerr << "NULL";
                    std::cerr << " at " << oldpos << " deleted or moved!\n";
                    #endif
                    if ( newpos[0] == newpos[0] ) {
                        newpos += VecD(0.5);
                        newpos[1] = 1 - newpos[1];
                        newpos *= imagesize;   // scale image up
                        newpos += imageborder;
                        oldpos += VecD(0.5);
                        oldpos[1] = 1 - oldpos[1];
                        oldpos *= imagesize;   // scale image up
                        oldpos += imageborder;
                        size_t id = reinterpret_cast<size_t>( datum );
                        out << "<animate"                            "\n"
                            << " xlink:href=\"#" << id  <<         "\"\n"
                            << " attributeName=\"cx\""               "\n"
                            << " fill =\"freeze\""                   "\n"
                            << " begin=\"" << DUR*currentTime <<  "s\"\n"
                            << " dur  =\"" << DUR             <<  "s\"\n"
                            << " from =\"" << oldpos[0]       << "px\"\n"
                            << " to   =\"" << newpos[0]       << "px\"\n"
                            << "/>"                                  "\n";
                        out << "<animate"                            "\n"
                            << " xlink:href=\"#" << id        <<   "\"\n"
                            << " attributeName=\"cy\""               "\n"
                            << " fill =\"freeze\""                   "\n"
                            << " begin=\"" << DUR*currentTime <<  "s\"\n"
                            << " dur  =\"" << DUR             <<  "s\"\n"
                            << " from =\"" << oldpos[1]       << "px\"\n"
                            << " to   =\"" << newpos[1]       << "px\"\n"
                            << "/>"                                  "\n";
                    } else { /* NAN means the datum wasn't found anymore */

                    }
                }
                ++i;
            }
            todo.pop();
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

    BoxesDrawnIt = BoxesDrawn.begin();
    while ( BoxesDrawnIt != BoxesDrawn.end() ) {
        VecD boxcenter = BoxesDrawnIt->first;
        int  id        = BoxesDrawnIt->second.id;
        if ( BoxesDrawnIt->second.visible and
             newtree.GetNodePtr( boxcenter ) == NULL )
        {
            #if DEBUG_OCTREE_SVG >= 10
            std::cerr << "Box around " << boxcenter << " vanished!!! "
                      << "ID: " << id << "\n";
            #endif
            out << "<set"                                     "\n"
                << " xlink:href=\"#" << id  <<              "\"\n"
                << " attributeName=\"stroke\""                "\n"
                << " begin=\"" << DUR*(currentTime+0.5) << "s\"\n"
                << " to   =\"none\""                          "\n"
                << "/>"                                       "\n";
            BoxesDrawnIt->second.visible = false;
        }
        ++BoxesDrawnIt;
    }

    currentTime += 1;
    this->tree = newtree;
}


} // namespace Octree


