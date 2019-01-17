#ifndef ANA_CONVEX_HULL_FUNCTIONS_H
#define ANA_CONVEX_HULL_FUNCTIONS_H

#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CH_into_cavity(
    Molecule const &protein, Cavity &hueco, ConvexHull const &CH);

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
Point get_intersection_point(Segment const &s, ConvexHull const &CH);

// Returns the intersection points between the segments of the tetrahedron
// delimited by the points and the convex hull. 3 points out, 1 in.
auto get_intersection_points_3out(ConvexHull const &CH, Point const &in_p0,
    Point const &out_p1, Point const &out_p2, Point const &out_p3)
    -> std::tuple<Point, Point, Point>;

// Returns the intersection points between the segments of the tetrahedron
// delimited by the points and the convex hull. 2 points out, 2 in.
auto get_intersection_points_2out(ConvexHull const &CH, Point const &in_p0,
    Point const &in_p1, Point const &out_p2, Point const &out_p3)
    -> std::tuple<Point, Point, Point, Point>;

// Returns the intersection points between the segments of the tetrahedron
// delimited by the points and the convex hull. 1 point out, 3 in.
auto get_intersection_points_1out(ConvexHull const &CH, Point const &in_p0,
    Point const &in_p1, Point const &in_p2, Point const &out_p3)
    -> std::tuple<Point, Point, Point>;

} // namespace ANA

#endif // _H