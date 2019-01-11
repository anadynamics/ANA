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

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream) {
    std::string temp;
    double coord = 0;

    in_stream >> temp;
    try {
        coord = std::stof(temp);
    } catch (const std::invalid_argument &ia) {
        // some character present
        throw std::runtime_error("Can't parse input. There may be "
                                 "non-numerical characters. Aborting.");
    } catch (...) {
        // some other exception.
        throw std::runtime_error("Can't parse input. Aborting.");
    }

    return coord;
}

class Molecule {
public:
    Molecule() = default;

    Molecule(unsigned int const natoms, unsigned int const nres) :
        _natoms(natoms), _nres(nres) {
        _data.reserve(natoms);
        _alphaCarbons.reserve(nres);
    }

    unsigned int _natoms, _nres;
    std::vector<std::pair<Point, VertexInfo>> _data;
    // Alpha carbon indices
    std::vector<unsigned int> _alphaCarbons;
    // Hetero-atoms indices
    std::vector<unsigned int> _hetatms;
};

}

#endif // _H