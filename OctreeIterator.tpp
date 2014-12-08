#pragma once

#define DEBUG_OCTREE_ITERATOR 0


namespace Octree {


/* used by begin() */
template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::iterator::iterator( void ) {}

template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::iterator::iterator( Node * root ) {
    tododata toBeStored = { /* child index */ 0, root };
    todo.push( toBeStored );
}

template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::iterator::~iterator( void ) {}

template<typename T_DTYPE, int T_DIM>
Node<T_DTYPE,T_DIM>::iterator::iterator( const iterator & src ) {
    this->todo = src.todo;
}

template<typename T_DTYPE, int T_DIM>
void Node<T_DTYPE,T_DIM>::iterator::operator=( const iterator & src ) {
    this->todo = src.todo;
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::iterator& Node<T_DTYPE,T_DIM>::iterator::operator++( void ) {
    while( !todo.empty() ) {
        const Node * currentNode = todo.top().node;
        if ( currentNode->IsLeaf() ) {
            todo.pop();
        } else {
            if ( currentNode->getChildPtr( todo.top().ichild ) != NULL ) {
                tododata toBeStored = { /* child index */ 0,
                               const_cast<Node*>( currentNode->getChildPtr( todo.top().ichild ) ) };
                todo.top().ichild++;
                todo.push(toBeStored);
                return *this;
            } else
                todo.pop();
        }
    }
    return *this;
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::iterator Node<T_DTYPE,T_DIM>::iterator::operator++( int unused ) {
    iterator tmp( *this );
	++(*this);
	return tmp;
}

template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::iterator::operator==( const iterator & src ) {
    assert(false);
}

template<typename T_DTYPE, int T_DIM>
bool Node<T_DTYPE,T_DIM>::iterator::operator!=( const iterator & src ) {
    return (this->todo.size()) != (src.todo.size());
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::Node & Node<T_DTYPE,T_DIM>::iterator::operator*( void ) const {
    return *(this->todo.top().node);
}

template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::Node * Node<T_DTYPE,T_DIM>::iterator::operator->( void ) const {
    return this->todo.top().node;
}

/* returns iterator with only root-node in todo stack */
template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::iterator Node<T_DTYPE,T_DIM>::begin( void ) {
    iterator it( this );
    return it;
}

/* returns iterator with empty stack */
template<typename T_DTYPE, int T_DIM>
typename Node<T_DTYPE,T_DIM>::iterator Node<T_DTYPE,T_DIM>::end( void ) {
    iterator it;
    return it;
}



} // namespace Octree
