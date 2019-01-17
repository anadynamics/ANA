#ifndef ANA_CGAL_UTILS_H
#define ANA_CGAL_UTILS_H

#include <ANA/Includes.hpp>

namespace ANA {

// Get the normal vector of the plane specified by p0, p1, p3.
inline Vector normal(Point const p0, Point const p1, Point const p2) {

    Vector const plane_vec_1 = p1 - p0;
    Vector const plane_vec_2 = p2 - p1;
    Vector plane_normal = CGAL::cross_product(plane_vec_1, plane_vec_2);
    plane_normal = plane_normal /
        std::sqrt(CGAL::to_double(plane_normal.squared_length()));

    return plane_normal;
}

// Just calculate volume of the cell.
inline double volume(Finite_cells_iterator const &cell_iterator) {
    return CGAL::to_double(CGAL::volume(cell_iterator->vertex(0)->point(),
        cell_iterator->vertex(1)->point(), cell_iterator->vertex(2)->point(),
        cell_iterator->vertex(3)->point()));
}

// Just calculate volume of the cell delimited by the 4 input points.
inline double volume(
    Point const &p0, Point const &p1, Point const &p2, Point const &p3) {
    return std::abs(CGAL::to_double(CGAL::volume(p0, p1, p2, p3)));
}

// Get the volume ocuppied by the sector of the sphere inscribed in the
// incident cell.
double sphere_sector_vol(Point const &p_0, Point const &p_1, Point const &p_2,
    Point const &p_3, double const radius);

} // namespace ANA

#endif // _H