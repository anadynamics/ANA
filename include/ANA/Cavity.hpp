#ifndef ANA_CAVITY_H
#define ANA_CAVITY_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

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

    friend void draw_lines(Cavity const &hueco, std::string const &filename);

    Delaunay _triangulation;
    std::vector<Finite_cells_iterator> _all_cells, _included_cells;
    double _volume = 0, _outer_volume = 0;

private:
    // Border polyhedrons from the cells that intersected the convex hull.
    std::vector<Tetrahedron> _tetra_border;
    std::vector<TriangularPrism> _penta_border;
    // Number of polyhedron vertices. Needed for output.
    std::vector<int> _poly_vtx_cnt;
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

} // namespace ANA

#endif // _H