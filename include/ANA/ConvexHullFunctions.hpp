#ifndef ANA_CONVEXHULLFUNCTIONS_H
#define ANA_CONVEXHULLFUNCTIONS_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CH_into_cavity(Cavity &hueco, ConvexHull const &CH);

// Discard cells without a vertex inside the specified convex hull. Lo
// precision. DEPRECATED.
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells);

// Discard cells without a vertex inside the specified convex hull. Hi
// precision.
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells, NA_Vector &out_intersecting_cells,
    std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &intersecting_total);

// Discard parts of cells outside the specified triangulation using
// intersections.
double discard_CH_1(NA_Vector const &in_intersecting_cells,
    Triang_Vector const &CH_triangs,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    const std::vector<int> &intersecting_total, Poly_Vector &border_poly,
    std::vector<std::array<double, 3>> &in_vtces_radii, int &atom_cnt_poly);

// Returns the intersection point between the segment and the convex hull.
inline Point get_intersection_point(Segment const &s, ConvexHull const &CH) {
    for (auto const &t : CH._triangles) {
        Object const inter_obj = CGAL::intersection(s, t);
        if (Point const *inter_point = CGAL::object_cast<Point>(&inter_obj)) {

            return *inter_point;
        }
    }
    throw("Fatal error, get_intersection_point() could not find an "
          "intersection.");
}

// Get intersection points between the cell and the included area.
auto get_vertices_3_out(Finite_cells_iterator const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, double>;

// Get intersection points between the cell and the included area.
auto get_vertices_2_out(Finite_cells_iterator const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double>;

// Get intersection points between the cell and the included area.
auto get_vertices_1_out(Finite_cells_iterator const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double,
        double>;

} // namespace ANA

#endif // _H