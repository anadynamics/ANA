#ifndef ANA_PRIMITIVES_H
#define ANA_PRIMITIVES_H

#include <ANA/Includes.hpp>

// ANA definitions
using MD_Element = std::array<Point, 4>;
// Cell
using MD_Vector = std::vector<MD_Element>;
// Pocket
using MD_Matrix = std::vector<MD_Vector>;
// All voids
using NDD_Element = std::array<std::pair<Point, double>, 4>;
// Cell
using NDD_Vector = std::vector<NDD_Element>;
// Pocket
using NDD_Matrix = std::vector<NDD_Vector>;
// All voids
using NDD_IElement = std::array<unsigned int, 4>;
// Cell indices
using NDD_IVector = std::vector<NDD_IElement>;
// Pocket indices
using NDD_IMatrix = std::vector<NDD_IVector>;
// All voids indices
using NA_Vector = std::vector<Finite_cells_iterator>;
// Pocket
using NA_Matrix = std::vector<NA_Vector>;
// All voids
using Poly_Vector = std::vector<Polyhedron>;
// Pocket border cells
using Poly_Matrix = std::vector<Poly_Vector>;
// All null areas border cells
using ANA_molecule = std::vector<std::pair<Point, VertexInfo>>;

namespace ANA {

struct Molecule {
public:
    std::vector<std::pair<Point, VertexInfo>> data;
};

}

#endif // _H