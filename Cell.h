#pragma once

namespace SimulationBox {

typedef struct CellDataStruct {
    int value;
    // standard copy assignment operator should suffice
} CellData;


} // namespace SimulationBox

/* Enables cout << Cell; This also works with fstream and therefore with tout */
ostream& operator<<( ostream& out, const SimulationBox::CellData cell ) {
    out << cell.value;
    return out;
}
