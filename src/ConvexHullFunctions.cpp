#include <ANA/ConvexHullFunctions.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CH_into_cavity(Cavity &hueco, ConvexHull const &CH) {

    for (auto const &cell : hueco._all_cells) {
        std::vector<int> vertices_out, vertices_in;
        // Only inside points will give a >0 dot product against all normals.
        for (std::size_t i = 0; i <= 3; ++i) {
            Point const test_point(cell->vertex(i)->point());

            for (std::size_t j = 0; j < CH._triangles.size(); j++) {
                Vector test_vtor = test_point - CH._triangles[1][j];
                test_vtor = test_vtor /
                    std::sqrt(CGAL::to_double(test_vtor.squared_length()));
                double const test_dot_pdt =
                    CGAL::to_double(test_vtor * CH._normals[j]);
                if (test_dot_pdt < 0) {
                    vertices_out.push_back(i);
                    break;
                } else {
                    vertices_in.push_back(i);
                }
            }
        }
        auto const vertices_out_cnt = vertices_out.size();
        if (vertices_out_cnt == 0) {
            // cell is entirely contained in the included area.
            hueco._included_cells.push_back(cell);
        } else if (vertices_out_cnt == 4) {
            // cell is outside the included area.
            continue;
        } else {
            // Cell has 1-3 vertices inside the included area.
            // Get the points of the intersections between the convex hull
            // and the tetrahedron's segments that go from the inner
            // point(s) to the outer point(s). The polyhedron delimited by these
            // points will be added to the cavity.
            switch (vertices_out_cnt) {
            case 3: {
                auto const vtx_in_0 = cell->vertex(vertices_in[0]);
                Point const p_in_0(vtx_in_0->point());
                Point const p_out_1(cell->vertex(vertices_out[0])->point());
                Point const p_out_2(cell->vertex(vertices_out[1])->point());
                Point const p_out_3(cell->vertex(vertices_out[2])->point());

                auto const [ip1, ip2, ip3] = get_intersection_points_3out(
                    CH, p_in_0, p_out_1, p_out_2, p_out_3);

                hueco.add_border_tetra(
                    p_in_0, ip1, ip2, ip3, vtx_in_0->info()._radius);

                break;
            }
            case 2: {
                auto const vtx_in_0 = cell->vertex(vertices_in[0]);
                Point const p_in_0(vtx_in_0->point());
                auto const vtx_in_1 = cell->vertex(vertices_in[1]);
                Point const p_in_1(vtx_in_1->point());
                Point const p_out_2(cell->vertex(vertices_out[0])->point());
                Point const p_out_3(cell->vertex(vertices_out[1])->point());

                auto const [ip1, ip2, ip3, ip4] = get_intersection_points_2out(
                    CH, p_in_0, p_in_1, p_out_2, p_out_3);

                hueco.add_border_penta(p_in_0, p_in_1, ip1, ip2, ip3, ip4,
                    vtx_in_0->info()._radius, vtx_in_1->info()._radius);
                break;
            }
            case 1: {
                auto const vtx_in_0 = cell->vertex(vertices_in[0]);
                Point const p_in_0(vtx_in_0->point());
                auto const vtx_in_1 = cell->vertex(vertices_in[1]);
                Point const p_in_1(vtx_in_1->point());
                auto const vtx_in_2 = cell->vertex(vertices_in[2]);
                Point const p_in_2(vtx_in_2->point());
                Point const p_out_3(cell->vertex(vertices_out[0])->point());

                auto const [ip1, ip2, ip3] = get_intersection_points_1out(
                    CH, p_in_0, p_in_1, p_in_2, p_out_3);

                hueco.add_border_penta(p_in_0, p_in_1, p_in_2, ip1, ip2, ip3,
                    vtx_in_0->info()._radius, vtx_in_1->info()._radius,
                    vtx_in_2->info()._radius);
                break;
            }
            }
        }
    }

    return;
}

// Returns the intersection point between the segment and the convex hull.
Point get_intersection_point(Segment const &s, ConvexHull const &CH) {
    for (auto const &t : CH._triangles) {
        Object const inter_obj = CGAL::intersection(s, t);
        if (Point const *inter_point = CGAL::object_cast<Point>(&inter_obj)) {

            return *inter_point;
        }
    }
    throw("Fatal error, get_intersection() could not find an intersection.");
}

// Returns the intersection points between the segments of the tetrahedron
// delimited by the points and the convex hull.
auto get_intersection_points_3out(ConvexHull const &CH, Point const &in_p0,
    Point const &out_p1, Point const &out_p2, Point const &out_p3)
    -> std::tuple<Point, Point, Point> {

    auto ip1 = get_intersection_point(Segment(in_p0, out_p1), CH);
    auto ip2 = get_intersection_point(Segment(in_p0, out_p2), CH);
    auto ip3 = get_intersection_point(Segment(in_p0, out_p3), CH);

    // Check if included area and intersecting cell are sharing vertices,
    // creating degeneracies. I don't think this is necessary and it probably
    // doesn't work. Remove it. And then make ip1/2/3 `const`. TODO
    if (ip1 == ip2) {
        ip1 = ip1 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip1 == ip3) {
        ip1 = ip1 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip2 == ip3) {
        ip3 = ip3 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }

    return {ip1, ip2, ip3};
}

// Returns the intersection points between the segments of the tetrahedron
// delimited by the points and the convex hull. 2 points out, 2 in.
auto get_intersection_points_2out(ConvexHull const &CH, Point const &in_p0,
    Point const &in_p1, Point const &out_p2, Point const &out_p3)
    -> std::tuple<Point, Point, Point, Point> {

    auto ip1 = get_intersection_point(Segment(in_p0, out_p2), CH);
    auto ip2 = get_intersection_point(Segment(in_p0, out_p3), CH);
    auto ip3 = get_intersection_point(Segment(in_p1, out_p2), CH);
    auto ip4 = get_intersection_point(Segment(in_p1, out_p3), CH);

    // Check if included area and intersecting cell are sharing vertices,
    // creating degeneracies. I don't think this is necessary and it probably
    // doesn't work. Remove it. And then make ip1/2/3 `const`. TODO
    if (ip1 == ip2) {
        ip1 = ip1 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip1 == ip3) {
        ip1 = ip1 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip1 == ip4) {
        ip1 = ip1 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip2 == ip3) {
        ip2 = ip2 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip2 == ip4) {
        ip2 = ip2 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip3 == ip4) {
        ip3 = ip3 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }

    return {ip1, ip2, ip3, ip4};
}

auto get_intersection_points_1out(ConvexHull const &CH, Point const &in_p0,
    Point const &in_p1, Point const &in_p2, Point const &out_p3)
    -> std::tuple<Point, Point, Point> {

    auto ip1 = get_intersection_point(Segment(in_p0, out_p3), CH);
    auto ip2 = get_intersection_point(Segment(in_p1, out_p3), CH);
    auto ip3 = get_intersection_point(Segment(in_p2, out_p3), CH);

    // Check if included area and intersecting cell are sharing vertices,
    // creating degeneracies. I don't think this is necessary and it probably
    // doesn't work. Remove it. And then make ip1/2/3 `const`. TODO
    if (ip1 == ip2) {
        ip1 = ip1 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip1 == ip3) {
        ip1 = ip1 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }
    if (ip2 == ip3) {
        ip2 = ip2 + Vector(0.01, 0.01, 0.01);
        printf("Degeneracy! Q crees.\n");
    }

    return {ip1, ip2, ip3};
}

} // namespace ANA