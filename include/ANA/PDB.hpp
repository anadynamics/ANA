#ifndef ANA_PDB_H
#define ANA_PDB_H
#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

void draw(ConvexHull const &CH, std::string const &filename);

void draw(Triangle const &t, FILE *out_file, int &idx, int &resid);

void connect_triangle(FILE *out_file, int const first_t, int const last_t);

void draw(Cavity const &hueco, std::string const &filename);

void draw(
    Finite_cells_iterator const &cell, FILE *out_file, int &idx, int &resid);

void draw(Polyhedron const &polyhedron, FILE *out_file, int &idx, int &resid);

void draw(Point const &punto, FILE *out_file, int constidx, int const resid);

void connect_cell(FILE *out_file, int const first_cell, int const last_cell);

} // namespace PDB
// namespace ANA
#endif