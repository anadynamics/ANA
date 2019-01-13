#ifndef ANA_PRIMITIVES_H
#define ANA_PRIMITIVES_H

#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>

namespace ANA {

class Molecule {
public:
    Molecule() = default;

    Molecule(std::string const &filename, bool const atom_only);

    unsigned int _natoms, _nres;
    std::vector<std::pair<Point, VertexInfo>> _data;
    // Alpha carbon indices
    std::vector<unsigned int> _alphaCarbons;
    // Hetero-atoms indices
    std::vector<unsigned int> _hetatms;
};

class Cavity {
public:
    Cavity() = default;

    Cavity(Delaunay T) noexcept : _triangulation(T) {}

    Cavity(Delaunay &&T) noexcept : _triangulation(T) {}

    Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts);

    Delaunay _triangulation;
    // Cells.
    std::vector<Finite_cells_iterator> _all, _included, _inner, _intersecting,
        _joint;
    double _volume = 0;
};

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream);

// Calculate area of the cell's facets
inline auto get_facets_areas(Finite_cells_iterator const &cell_iterator)
    -> std::tuple<double, double, double, double> {

    double const f0 = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(1)->point(), cell_iterator->vertex(2)->point(),
        cell_iterator->vertex(3)->point()));
    double const f1 = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(2)->point(), cell_iterator->vertex(3)->point(),
        cell_iterator->vertex(0)->point()));
    double const f2 = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(0)->point(), cell_iterator->vertex(1)->point(),
        cell_iterator->vertex(3)->point()));
    double const f3 = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(0)->point(), cell_iterator->vertex(1)->point(),
        cell_iterator->vertex(2)->point()));

    return {f0, f1, f2, f3};
}
// Determine if any of the facets of the given cells has area larger than the
// criteria.
inline bool refine_cell_areas(
    Finite_cells_iterator const cell_iterator, double const top) {

    auto const [f0, f1, f2, f3] = get_facets_areas(cell_iterator);
    bool const larger = (f0 > top || f1 > top || f2 > top || f3 > top);

    return larger;
}

// Just calculate volume of the cell.
inline double cell_volume(Finite_cells_iterator const &cell_iterator) {
    return CGAL::to_double(CGAL::volume(cell_iterator->vertex(0)->point(),
        cell_iterator->vertex(1)->point(), cell_iterator->vertex(2)->point(),
        cell_iterator->vertex(3)->point()));
}

}

#endif // _H