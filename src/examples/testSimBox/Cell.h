#pragma once

#include <iostream>

typedef struct CellDataStruct {
    int value;
    // standard copy assignment operator should suffice
} CellData;

/* Enables cout << Cell; This also works with fstream and therefore with tout */
std::ostream& operator<<( std::ostream& out, const CellData cell ) {
    out << cell.value;
    return out;
}
