#ifndef ANA_PRIMITIVES_H
#define ANA_PRIMITIVES_H

#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>

namespace ANA {

class Molecule {
public:
    Molecule() = default;

    Molecule(std::string const &filename, bool const atom_only);

    int _natoms, _nres;
    std::vector<std::pair<Point, VertexInfo>> _data;
    // Alpha carbon indices
    std::vector<int> _alphaCarbons;
    // Hetero-atoms indices
    std::vector<int> _hetatms;
};

class Cavity {
public:
    Cavity() = default;

    Cavity(Delaunay T) noexcept : _triangulation(T) {}

    Cavity(Delaunay &&T) noexcept : _triangulation(T) {}

    Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts);

    void add_border_tetra(Point const &p0, Point const &ip1, Point const &ip2,
        Point const &ip3, double const vdw0);

    void add_border_penta(Point const &p0, Point const &p1, Point const &ip2,
        Point const &ip3, Point const &ip4, Point const &ip5, double const vdw0,
        double const vdw1);

    void add_border_penta(Point const &p0, Point const &p1, Point const &p2,
        Point const &ip3, Point const &ip4, Point const &ip5, double const vdw0,
        double const vdw1, double const vdw2);

    friend void draw(Cavity const &hueco, std::string const &filename);

    Delaunay _triangulation;
    std::vector<Finite_cells_iterator> _all_cells, _included_cells;
    double _volume = 0, _outer_volume = 0;

private:
    // Border polyhedrons from the cells that intersected the convex hull.
    std::vector<Polyhedron> _border;
    // Number of polyhedron vertices. Needed for output.
    int _poly_vtx_cnt = 0;
};

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream);

// Calculate area of the cell's facets
inline auto get_facets_areas(Finite_cells_iterator const cell_iterator)
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

}

#endif // _H