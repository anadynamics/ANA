#include <ANA/NDDUtils.hpp>
#include <ANA/Utils.hpp>

namespace ANA {

// Cluster neighbouring cells. CGAL neighbouring method
void cluster_cells_cgal(NA_Vector const &input_cells, NA_Matrix &output_cells,
    unsigned int const min_cells_cluster) {

    unsigned int nbr_of_clusters = 0, cell_cnt = input_cells.size();
    std::vector<unsigned int> already_checked;
    Finite_cells_iterator fc_ite;
    NA_Vector neighbors;

    for (unsigned int i = 0; i <= cell_cnt; ++i) {
        if (std::find(already_checked.begin(), already_checked.end(), i) !=
            already_checked.end()) {
            continue;
        }
        already_checked.push_back(i);
        get_neighbors(input_cells, input_cells[i], already_checked, neighbors);
        if (neighbors.empty()) {
            continue;
        }
        // Now iterate over each of its neighbors and its neighbors's neighbors
        for (Finite_cells_iterator const &fc_ite2 : neighbors) {
            get_neighbors(input_cells, fc_ite2, already_checked, neighbors);
        }
        if (neighbors.size() >= min_cells_cluster) {
            ++nbr_of_clusters;
            output_cells.push_back(neighbors);
        }
        neighbors.clear();
    }
    return;
}

// Given a cell, get all neighbouring cells that haven't been discovered yet
void get_neighbors(NA_Vector const &input_cells,
    Finite_cells_iterator const &query_cell, std::vector<unsigned int> &except,
    NA_Vector &output_cells) {

    unsigned int i, cnt = 0, cell_cnt = input_cells.size();
    for (i = 0; i <= cell_cnt; ++i) {
        if (find(except.begin(), except.end(), i) != except.end()) {
            continue;
        }
        if (query_cell->has_neighbor(input_cells[i])) {
            except.push_back(i);
            output_cells.push_back(input_cells[i]);
            ++cnt;
        }
    }
    return;
}

// Functor to get pairs of intersecting boxes and group them in clusters.
class callback_for_boxes {
public:
    std::vector<unsigned int> *_id_vector_a, *_id_vector_b;
    NA_Vector *_neighbors_a, *_neighbors_b;

    callback_for_boxes(std::vector<unsigned int> &id_vector_a,
        std::vector<unsigned int> &id_vector_b, NA_Vector &neighbors_a,
        NA_Vector &neighbors_b) :
        _id_vector_a(&id_vector_a),
        _id_vector_b(&id_vector_b), _neighbors_a(&neighbors_a),
        _neighbors_b(&neighbors_b) {}

    void operator()(const Box &a, const Box &b) {
        // New boxes. Store their indices and their cells.
        _id_vector_a->push_back(a.id());
        _id_vector_b->push_back(b.id());
        _neighbors_a->push_back(a.handle());
        _neighbors_b->push_back(b.handle());

        return;
    }
};

// Cluster neighbouring cells. Iso-oriented boxes method
void cluster_cells_boxes(
    NA_Vector const &input_cells, NA_Matrix &output_cells) {

    unsigned int min_cells_cluster = 2;
    std::vector<Box> boxes;
    std::vector<unsigned int> id_vector_a, id_vector_b;
    std::vector<std::vector<unsigned int>> id_mtx;
    NA_Vector neighbors_a, neighbors_b;
    callback_for_boxes called(
        id_vector_a, id_vector_b, neighbors_a, neighbors_b);

    // Get cells bounding boxes.
    for (auto const &cell_ite : input_cells) {
        const Tetrahedron tetra(cell_ite->vertex(0)->point(),
            cell_ite->vertex(1)->point(), cell_ite->vertex(2)->point(),
            cell_ite->vertex(3)->point());
        boxes.push_back(Box(std::move(tetra.bbox()), cell_ite));
    }

    // Get the intersections
    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), called);

    // Get the indices that sort the id vectors
    std::vector<unsigned int> indices_a = sort_indices(id_vector_a);
    std::vector<unsigned int> indices_b = sort_indices(id_vector_b);

    // Get a sorted vector with all the ids and no repetitions. "map_aux" is
    // probably bigger than needed, since some cells won't have any
    // intersections
    // at all and won't appear in neither "indices_a" nor "indices_b".
    std::vector<unsigned int> map_aux(input_cells.size());
    map_aux[0] = id_vector_a[indices_a[0]];
    unsigned int cnt = 0, tmp = map_aux[cnt];
    // Get all the indices IDs from "id_vector_a".
    for (auto const &each : indices_a) {
        if (id_vector_a[each] != tmp) {
            ++cnt;
            map_aux[cnt] = id_vector_a[each];
            tmp = map_aux[cnt];
        }
    }
    // Now, if "id_vector_a" didn't have all the IDs, there will be 0s in the
    // "map_aux" vector.
    auto ite = map_aux.end() - 1;
    while (*ite == 0) {
        map_aux.pop_back();
        ite = map_aux.end() - 1;
    }

    cnt = 0;
    ite = std::lower_bound(
        map_aux.begin(), map_aux.end(), id_vector_b[indices_b[0]]);

    for (auto const &each : indices_b) {
        if (id_vector_b[each] != *ite) {
            // This may be a new ID or just a higher ID that's not were "ite"
            // points,
            // but forward.
            ite = std::lower_bound(
                map_aux.begin(), map_aux.end(), id_vector_b[each]);
            if (ite == map_aux.end()) {
                // Don't dereference the end iterator
                map_aux.push_back(id_vector_b[each]);
                ite = map_aux.end() - 1;
            } else if (id_vector_b[each] != *ite) {
                ite = map_aux.insert(ite, id_vector_b[each]);
            }
        }
    }
    // map_aux, together with map_bool serve as a "check and lookup" table.
    std::vector<bool> map_bool(map_aux.size());

    // Start mapping the clusters.
    unsigned int result = 0, new_id;
    for (size_t i = 0; i < map_aux.size(); i++) {

        if (map_bool[i]) {
            // Already added.
            continue;
        }

        std::vector<unsigned int> current, current_tmp;
        NA_Vector neighbors;
        // This is the starting cell of the cluster.
        current_tmp.push_back(map_aux[i]);

        if (lb_with_indices(id_vector_a, indices_a, map_aux[i], result) &&
            (id_vector_a[indices_a[result]] == map_aux[i])) {
            neighbors.push_back(neighbors_a[indices_a[result]]);

        } else if (lb_with_indices(
                       id_vector_b, indices_b, map_aux[i], result) &&
            (id_vector_b[indices_b[result]] == map_aux[i])) {
            neighbors.push_back(neighbors_b[indices_b[result]]);
        } else {
            std::cerr << "Iso-oriented boxes clustering method. Invalid "
                         "starting cell ID."
                      << '\n';
            std::exit(1);
        }

        while (current_tmp.size() != 0) {
            auto const box = current_tmp[0];

            // Search for intersections of this box in the "B" array.
            if (lb_with_indices(id_vector_a, indices_a, box, result)) {
                // Keep moving further until every intersections on this array
                // has
                // been checked.
                while (id_vector_a[indices_a[result]] == box) {

                    new_id = id_vector_b[indices_a[result]];
                    // Check if its a new box.
                    if ((std::find(current_tmp.begin(), current_tmp.end(),
                             new_id) == current_tmp.end()) &&
                        (std::find(current.begin(), current.end(), new_id) ==
                            current.end())) {
                        // Store it and its subsequent intersections.
                        current_tmp.push_back(new_id);
                        neighbors.push_back(
                            std::move(neighbors_b[indices_a[result]]));
                    }

                    if (++result >= indices_a.size()) {
                        // Don't overflow indices_a vector.
                        break;
                    }
                }
            }

            // Search for other intersections of this box in the "A" array.
            if (lb_with_indices(id_vector_b, indices_b, box, result)) {
                // Keep moving further until every intersections on this array
                // has
                // been checked.
                while (id_vector_b[indices_b[result]] == box) {

                    new_id = id_vector_a[indices_b[result]];
                    // Check if its a new box.
                    if ((std::find(current_tmp.begin(), current_tmp.end(),
                             new_id) == current_tmp.end()) &&
                        (std::find(current.begin(), current.end(), new_id) ==
                            current.end())) {
                        // Store it and its subsequent intersections.
                        current_tmp.push_back(new_id);
                        neighbors.push_back(
                            std::move(neighbors_a[indices_b[result]]));
                    }

                    if (++result >= indices_b.size()) {
                        // Don't overflow indices_a vector.
                        break;
                    }
                }
            }

            // Done with this box intersections. Move it to the list of already
            // checked boxes and check it.
            auto const done_box_ite = current_tmp.begin();
            std::move(
                done_box_ite, done_box_ite + 1, std::back_inserter(current));
            current_tmp.erase(done_box_ite, done_box_ite + 1);
            lb(map_aux, box, result);
            map_bool[result] = true;
        }

        // Done with this cluster.
        if (current.size() >= min_cells_cluster) {
            id_mtx.push_back(std::move(current));
            output_cells.push_back(std::move(neighbors));
        }
    }

    return;
}

// Keep cells that correspond to the included amino acids
void keep_included_aa_cells(NA_Vector const &input_cells,
    const std::vector<unsigned int> &aa_list,
    unsigned int const nbr_of_vertices_to_include, NA_Vector &output_cells) {
    unsigned int i, nbr_of_vertices;
    Point p1;

    for (Finite_cells_iterator const &fc_ite : input_cells) {
        nbr_of_vertices = 0;
        for (i = 0; i <= 3; ++i) {
            if (std::binary_search(aa_list.begin(), aa_list.end(),
                    fc_ite->vertex(i)->info()._resn)) {
                ++nbr_of_vertices;
                if (nbr_of_vertices >= nbr_of_vertices_to_include) {
                    output_cells.push_back(fc_ite);
                    break;
                }
            }
        }
    }
}

// Discard exposed cells.
void discard_ASA_dot_pdt_cm(Point const &cm,
    std::vector<Point> const &Calpha_xyz, double const min_dot,
    double const max_length, std::string const only_side_ASA,
    NA_Vector const &input_cells, NA_Vector &output_cells) {

    if (only_side_ASA == "inside") {
        for (Finite_cells_iterator const &cell_ite : input_cells) {

            Point test_point = CGAL::centroid(cell_ite->vertex(0)->point(),
                cell_ite->vertex(1)->point(), cell_ite->vertex(2)->point(),
                cell_ite->vertex(3)->point());
            Vector diff_cm = test_point - cm;
            diff_cm = diff_cm /
                (std::sqrt(CGAL::to_double(diff_cm.squared_length())));

            for (Point const &current_point : Calpha_xyz) {
                // if ( current_point == test_point ) continue; // only when
                // checking
                // vtces
                Vector diff_1 = current_point - test_point;
                double diff_norm =
                    std::sqrt(CGAL::to_double(diff_1.squared_length()));
                if (diff_norm > max_length)
                    continue;
                diff_1 = diff_1 / diff_norm;
                double dot_pdt = CGAL::to_double(diff_1 * diff_cm);
                if (dot_pdt > min_dot) {
                    output_cells.push_back(cell_ite);
                    break;
                }
            }
        }
    } else if (only_side_ASA == "outside") {
        for (Finite_cells_iterator const &cell_ite : input_cells) {

            double top_dot_pdt = -1.0;
            Point test_point = CGAL::centroid(cell_ite->vertex(0)->point(),
                cell_ite->vertex(1)->point(), cell_ite->vertex(2)->point(),
                cell_ite->vertex(3)->point());
            Vector diff_cm = test_point - cm;
            diff_cm = diff_cm /
                (std::sqrt(CGAL::to_double(diff_cm.squared_length())));

            for (Point const &current_point : Calpha_xyz) {
                // if ( current_point == test_point ) continue; // only when
                // checking
                // vtces
                Vector diff_1 = current_point - test_point;
                double diff_norm =
                    std::sqrt(CGAL::to_double(diff_1.squared_length()));
                if (diff_norm > max_length)
                    continue;
                diff_1 = diff_1 / diff_norm;
                double dot_pdt = CGAL::to_double(diff_1 * diff_cm);
                if (dot_pdt > top_dot_pdt) {
                    top_dot_pdt = dot_pdt;
                }
            }
            if (top_dot_pdt < min_dot) {
                output_cells.push_back(cell_ite);
            }
        }
    }

    return;
}

// Discard exposed cells. Calpha convex hull method
void discard_ASA_CACH(std::vector<Point> const &Calpha_xyz,
    std::string const only_side_ASA, NA_Vector const &input_cells,
    NA_Vector &output_cells) {

    Polyhedron CH;
    std::vector<Vector> CH_normals;
    std::vector<Point> CH_vtces;

    // Get the convex hull determined by the Calphas
    CGAL::convex_hull_3(Calpha_xyz.begin(), Calpha_xyz.end(), CH);
    // Triangle normals point inwards. Only inside points will give a positive
    // dot
    // product against all normals
    P_Facet_const_iterator f_end = CH.facets_end();
    for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
         ++f_ite) {
        P_Halfedge_around_facet_const_circulator he_ite = f_ite->facet_begin();
        // For some reason I have to store the poins first. There's some C++
        // programming aspect I'm missing here
        Point p0((++he_ite)->vertex()->point());
        Point p1((++he_ite)->vertex()->point());
        Point p2((++he_ite)->vertex()->point());
        Vector v1 = p1 - p0;
        Vector v2 = p2 - p1;
        Vector normal = CGAL::cross_product(v2, v1);
        normal = normal / std::sqrt(CGAL::to_double(normal.squared_length()));
        CH_normals.push_back(normal);
        CH_vtces.push_back(he_ite->vertex()->point());
    }

    // Now discard outside cells
    for (auto const &cell_ite : input_cells) {
        bool vtx_inside_bool = false;
        Point test_point = CGAL::centroid(cell_ite->vertex(0)->point(),
            cell_ite->vertex(1)->point(), cell_ite->vertex(2)->point(),
            cell_ite->vertex(3)->point());

        for (unsigned int j = 0; j < CH_vtces.size(); j++) {
            Vector test_vtor = test_point - CH_vtces[j];
            test_vtor = test_vtor /
                std::sqrt(CGAL::to_double(test_vtor.squared_length()));
            double test_dot_pdt = CGAL::to_double(test_vtor * CH_normals[j]);
            if (test_dot_pdt < 0) {
                // This centroid lies outside the Calpha convex hull
                vtx_inside_bool = true;
                break;
            }
        }
        if (!vtx_inside_bool && only_side_ASA == "inside") {
            // Cell has its centroid inside de convex hull determined by the
            // Calphas
            output_cells.push_back(cell_ite);
        } else if (vtx_inside_bool && only_side_ASA == "outside") {
            // Cell has its centroid outside de convex hull determined by the
            // Calphas
            output_cells.push_back(cell_ite);
        }
    }

    return;
}

// Discard exposed cells. axes method
void discard_ASA_dot_pdt_axes(std::vector<Point> const &Calpha_xyz,
    double const min_dot, double const max_length,
    std::string const only_side_ASA, NA_Vector const &input_cells,
    NA_Vector &output_cells) {

    double DP_xy = 0, DP_xz = 0, DP_yz = 0;

    if (only_side_ASA == "inside") {
        for (auto const &cell_ite : input_cells) {

            // These will indicate wether a half-axis intersects a Calpha atom
            bool bool_half_xy = 0, bool_half_yx = 0, bool_half_xz = 0,
                 bool_half_zx = 0, bool_half_yz = 0, bool_half_zy = 0;
            // Get centroid
            Point const ctd = CGAL::centroid(cell_ite->vertex(0)->point(),
                cell_ite->vertex(1)->point(), cell_ite->vertex(2)->point(),
                cell_ite->vertex(3)->point());
            const Vector Vctd = ctd - CGAL::ORIGIN;
            const Vector Vx = Vctd + Vector(1, 0, 0);
            const Vector Vy = Vctd + Vector(0, 1, 0);
            const Vector Vz = Vctd + Vector(0, 0, 1);
            // Get xy, xz and yz planes normal vectors
            Vector Vxy = CGAL::cross_product(Vx, Vy);
            Vector Vxz = CGAL::cross_product(Vx, Vz);
            Vector Vyz = CGAL::cross_product(Vy, Vz);
            // Normalize them
            Vxy = Vxy / (std::sqrt(CGAL::to_double(Vxy.squared_length())));
            Vxz = Vxz / (std::sqrt(CGAL::to_double(Vxz.squared_length())));
            Vyz = Vyz / (std::sqrt(CGAL::to_double(Vyz.squared_length())));

            // Now, find Calpha atoms that intersect (w/ a margin of error)
            // these
            // vectors
            for (auto const &CA : Calpha_xyz) {

                // Get difference vector from centroid to current Calpha
                Vector Vdiff = CA - ctd;
                double Vdiff_norm =
                    std::sqrt(CGAL::to_double(Vdiff.squared_length()));
                if (Vdiff_norm > max_length)
                    continue;
                Vdiff = Vdiff / Vdiff_norm;
                // Get the dot products
                if (!(bool_half_xy && bool_half_yx)) {
                    DP_xy = CGAL::to_double(Vdiff * Vxy);
                }
                if (!(bool_half_xz && bool_half_zx)) {
                    DP_xz = CGAL::to_double(Vdiff * Vxz);
                }
                if (!(bool_half_yz && bool_half_zy)) {
                    DP_yz = CGAL::to_double(Vdiff * Vyz);
                }

                // Determine if there is an intersection
                if (DP_xy >= min_dot)
                    bool_half_xy = true;
                if (DP_xy <= (-min_dot))
                    bool_half_yx = true;
                if (DP_xz >= min_dot)
                    bool_half_xz = true;
                if (DP_xz <= (-min_dot))
                    bool_half_zx = true;
                if (DP_yz >= min_dot)
                    bool_half_yz = true;
                if (DP_yz <= (-min_dot))
                    bool_half_zy = true;

                // Determine if there is enough of them to classify this cell as
                // being
                // inside
                unsigned int const sum_halves = bool_half_xy + bool_half_yx +
                    bool_half_xz + bool_half_zx + bool_half_yz + bool_half_zy;
                if (sum_halves < 4) {
                    // can't say yet
                    continue;
                } else if (sum_halves > 4) {
                    // inside
                    output_cells.push_back(std::move(cell_ite));
                    break;
                } else if (sum_halves == 4) {
                    // its in a "tunnel" or can't say yet
                    if ((!bool_half_xy && !bool_half_yx && bool_half_xz &&
                            bool_half_zx && bool_half_yz && bool_half_zy) ||
                        (bool_half_xy && bool_half_yx && !bool_half_xz &&
                            !bool_half_zx && bool_half_yz && bool_half_zy) ||
                        (bool_half_xy && bool_half_yx && bool_half_xz &&
                            bool_half_zx && !bool_half_yz && !bool_half_zy)) {
                        // inside
                        output_cells.push_back(std::move(cell_ite));
                        break;
                    }
                }
            }
        }
    } else if (only_side_ASA == "outside") {
        for (auto const &cell_ite : input_cells) {

            // These will indicate wether a half-axis intersects a Calpha atom
            bool bool_half_xy = 0, bool_half_yx = 0, bool_half_xz = 0,
                 bool_half_zx = 0, bool_half_yz = 0, bool_half_zy = 0,
                 bool_inside = 0;
            // Get centroid
            Point const ctd = CGAL::centroid(cell_ite->vertex(0)->point(),
                cell_ite->vertex(1)->point(), cell_ite->vertex(2)->point(),
                cell_ite->vertex(3)->point());
            const Vector Vctd = ctd - CGAL::ORIGIN;
            const Vector Vx = Vctd + Vector(1, 0, 0);
            const Vector Vy = Vctd + Vector(0, 1, 0);
            const Vector Vz = Vctd + Vector(0, 0, 1);
            // Get xy, xz and yz planes normal vectors
            Vector Vxy = CGAL::cross_product(Vx, Vy);
            Vector Vxz = CGAL::cross_product(Vx, Vz);
            Vector Vyz = CGAL::cross_product(Vy, Vz);
            // Normalize them
            Vxy = Vxy / (std::sqrt(CGAL::to_double(Vxy.squared_length())));
            Vxz = Vxz / (std::sqrt(CGAL::to_double(Vxz.squared_length())));
            Vyz = Vyz / (std::sqrt(CGAL::to_double(Vyz.squared_length())));

            // Now, find Calpha atoms that intersect (w/ a margin of error)
            // these
            // vectors
            for (auto const &CA : Calpha_xyz) {

                // Get difference vector from centroid to current Calpha
                Vector Vdiff = CA - ctd;
                double Vdiff_norm =
                    std::sqrt(CGAL::to_double(Vdiff.squared_length()));
                if (Vdiff_norm > max_length)
                    continue;
                Vdiff = Vdiff / Vdiff_norm;
                // Get the dot products
                if (!(bool_half_xy && bool_half_yx)) {
                    DP_xy = CGAL::to_double(Vdiff * Vxy);
                }
                if (!(bool_half_xz && bool_half_zx)) {
                    DP_xz = CGAL::to_double(Vdiff * Vxz);
                }
                if (!(bool_half_yz && bool_half_zy)) {
                    DP_yz = CGAL::to_double(Vdiff * Vyz);
                }

                // Determine if there is an intersection
                if (DP_xy >= min_dot)
                    bool_half_xy = true;
                if (DP_xy <= -min_dot)
                    bool_half_yx = true;
                if (DP_xz >= min_dot)
                    bool_half_xz = true;
                if (DP_xz <= -min_dot)
                    bool_half_zx = true;
                if (DP_yz >= min_dot)
                    bool_half_yz = true;
                if (DP_yz <= -min_dot)
                    bool_half_zy = true;

                // Determine if there is enough of them to classify this cell as
                // being
                // inside
                unsigned int const sum_halves = bool_half_xy + bool_half_yx +
                    bool_half_xz + bool_half_zx + bool_half_yz + bool_half_zy;

                if (sum_halves > 4) {
                    // inside
                    bool_inside = true;
                    break;
                } else if (sum_halves == 4) {
                    // its in a "tunnel" or can't say yet
                    if ((!bool_half_xy && !bool_half_yx && bool_half_xz &&
                            bool_half_zx && bool_half_yz && bool_half_zy) ||
                        (bool_half_xy && bool_half_yx && !bool_half_xz &&
                            !bool_half_zx && bool_half_yz && bool_half_zy) ||
                        (bool_half_xy && bool_half_yx && bool_half_xz &&
                            bool_half_zx && !bool_half_yz && !bool_half_zy)) {
                        // inside
                        bool_inside = true;
                        break;
                    }
                }
            }
            if (!bool_inside) {
                output_cells.push_back(std::move(cell_ite));
                bool_inside = false;
            }
        }
    }

    return;
}
// Fill the 2 input vectors with iterators for the outer and inner cells
// respectively
void partition_triangulation(
    Delaunay const &T, NA_Vector &outer_cells, NA_Vector &inner_cells) {

    Finite_cells_iterator fc_ite, fc_ite_end = T.finite_cells_end();

    for (fc_ite = T.finite_cells_begin(); fc_ite != fc_ite_end; ++fc_ite) {
        if (T.is_infinite(fc_ite->neighbor(0)) ||
            T.is_infinite(fc_ite->neighbor(1)) ||
            T.is_infinite(fc_ite->neighbor(2)) ||
            T.is_infinite(fc_ite->neighbor(3))) {
            outer_cells.push_back(fc_ite);
        } else {
            inner_cells.push_back(fc_ite);
        }
    }
    return;
}
// Calc volume and get the proper cells
double get_all_voids(Delaunay const &T, NA_Vector &big_cells,
    CellFilteringOptions const cell_opts) {

    double const min_cell_vol = (4 / 3) * M_PI * pow(cell_opts._minVR, 3);
    double const max_facet_area = M_PI * pow(cell_opts._maxSR, 2);
    double volume = 0;
    double current_cell_vol;
    Finite_cells_iterator fc_ite, fc_ite_end = T.finite_cells_end();

    for (fc_ite = T.finite_cells_begin(); fc_ite != fc_ite_end; fc_ite++) {
        current_cell_vol = cell_volume(fc_ite);
        if (current_cell_vol > min_cell_vol &&
            refine_cell_areas(fc_ite, max_facet_area) == 0) {
            big_cells.push_back(fc_ite);
            volume = volume + current_cell_vol;
        }
    }
    return volume;
}

// Substract the volume filled with the 4 atoms from the total volume of the
// corresponding cell
double refine_cell_volume(
    double const entire_cell_vol, Finite_cells_iterator const &cell_iterator) {
    double const rdW_0 = double(cell_iterator->vertex(0)->info()._radius);
    double const rdW_1 = double(cell_iterator->vertex(1)->info()._radius);
    double const rdW_2 = double(cell_iterator->vertex(2)->info()._radius);
    double const rdW_3 = double(cell_iterator->vertex(3)->info()._radius);
    Point const p_0 = cell_iterator->vertex(0)->point();
    Point const p_1 = cell_iterator->vertex(1)->point();
    Point const p_2 = cell_iterator->vertex(2)->point();
    Point const p_3 = cell_iterator->vertex(3)->point();
    double const vertex_0_sphere_sector_vol =
        sphere_sector_vol(p_0, p_1, p_2, p_3, rdW_0);
    double const vertex_1_sphere_sector_vol =
        sphere_sector_vol(p_3, p_0, p_1, p_2, rdW_3);
    double const vertex_2_sphere_sector_vol =
        sphere_sector_vol(p_2, p_3, p_0, p_1, rdW_2);
    double const vertex_3_sphere_sector_vol =
        sphere_sector_vol(p_1, p_2, p_3, p_0, rdW_1);
    return (entire_cell_vol - vertex_0_sphere_sector_vol -
        vertex_1_sphere_sector_vol - vertex_2_sphere_sector_vol -
        vertex_3_sphere_sector_vol);
}
// Get the volume ocuppied by the sector of the sphere inscribed in the
// incident
// cell
double sphere_sector_vol(Point const &p_0, Point const &p_1, Point const &p_2,
    Point const &p_3, double const radius) {
    // get 1st point of the mini tetrahedron
    Vector vec_1 = p_1 - p_0;
    vec_1 =
        vec_1 / (std::sqrt(CGAL::to_double(vec_1.squared_length()))) * radius;
    Point point_1 = p_0 + vec_1;
    // get 2nd point of the mini tetrahedron
    Vector vec_2 = p_2 - p_0;
    vec_2 =
        vec_2 / (std::sqrt(CGAL::to_double(vec_2.squared_length()))) * radius;
    Point point_2 = p_0 + vec_2;
    // get 3rd point of the mini tetrahedron
    Vector vec_3 = p_3 - p_0;
    vec_3 =
        vec_3 / (std::sqrt(CGAL::to_double(vec_3.squared_length()))) * radius;
    Point point_3 = p_0 + vec_3;

    // Now, get the volume of the sphere slice
    // get the normal vector
    Vector plane_vec_1 = p_2 - p_1;
    Vector plane_vec_2 = p_3 - p_2;
    Vector plane_normal = CGAL::cross_product(plane_vec_1, plane_vec_2);
    // normalize the normal vector
    plane_normal = plane_normal /
        (std::sqrt(CGAL::to_double(plane_normal.squared_length())));
    // get the distance between the delaunay vertex and the plane formed by p_1,
    // p_2 and p_3
    double dist_to_plane = CGAL::to_double(plane_normal * vec_1);
    double h = radius - dist_to_plane;
    double spherical_cap_volume =
        1 / 3 * M_PI * std::pow(h, 2) * (3 * radius - h);
    double mini_tetrahedron_vol =
        CGAL::to_double(CGAL::volume(p_0, point_1, point_2, point_3));

    // Add the 2 volumes that represent the space occupied by the atom with
    // coordinates p0
    double volume =
        std::abs(mini_tetrahedron_vol) + std::abs(spherical_cap_volume);
    if (isnan(volume)) {
        volume = 0;
        //    std::cerr << "Warning: sphere_sector_vol gave NaN value. Void
        //    calculation "
        //                 "keeps going."
        //              << '\n';
    }

    return volume;
}

// Discard cells without a vertex inside the specified convex hull. Lo
// precision.
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells) {

    std::vector<Vector> CH_normals;
    std::vector<Point> CH_vtces;
    // Triangle normals point inwards. Only inside points will give a positive
    // dot product against all normals
    for (auto const &triangle : CH_triangs) {
        Vector v1 = triangle.vertex(1) - triangle.vertex(0);
        Vector v2 = triangle.vertex(2) - triangle.vertex(1);
        Vector normal = CGAL::cross_product(v2, v1);
        normal = normal / std::sqrt(CGAL::to_double(normal.squared_length()));
        CH_normals.push_back(normal);
        CH_vtces.push_back(triangle.vertex(1));
    }

    // Now discard outside cells
    for (auto const &cell_ite : in_cells) {
        for (unsigned int i = 0; i <= 3; ++i) {
            bool vtx_inside_bool = false;
            Point test_point(cell_ite->vertex(i)->point());

            for (unsigned int j = 0; j < CH_vtces.size(); ++j) {
                Vector test_vtor = test_point - CH_vtces[j];
                test_vtor = test_vtor /
                    std::sqrt(CGAL::to_double(test_vtor.squared_length()));
                double test_dot_pdt =
                    CGAL::to_double(test_vtor * CH_normals[j]);

                if (test_dot_pdt < 0) {
                    vtx_inside_bool = true;
                    break;
                }
            }
            if (!vtx_inside_bool) {
                // Cell has at least 1 vertex inside included area. Keep it and
                // go on
                out_cells.push_back(cell_ite);
                break;
            }
        }
    }

    return;
}

// Discard cells without a vertex inside the specified convex hull. Hi
// precision
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells, NA_Vector &out_intersecting_cells,
    std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<unsigned int> &intersecting_total) {

    std::vector<Vector> CH_normals;
    std::vector<Point> CH_vtces;
    // Triangle normals point inwards. Only inside points will give a positive
    // dot product against all normals.
    for (auto const &triangle : CH_triangs) {
        Vector v1 = triangle.vertex(1) - triangle.vertex(0);
        Vector v2 = triangle.vertex(2) - triangle.vertex(1);
        Vector normal = CGAL::cross_product(v2, v1);
        normal = normal / std::sqrt(CGAL::to_double(normal.squared_length()));
        CH_normals.push_back(normal);
        CH_vtces.push_back(triangle.vertex(1));
    }

    // Now discard outside cells
    for (auto const &cell_ite : in_cells) {
        std::array<bool, 4> vtx_inside_bool = {false, false, false, false};
        unsigned int total = 0;
        // Get nbr of vertices that lie outside the include area
        for (unsigned int i = 0; i <= 3; ++i) {
            Point test_point(cell_ite->vertex(i)->point());

            for (unsigned int j = 0; j < CH_vtces.size(); j++) {
                Vector test_vtor = test_point - CH_vtces[j];
                test_vtor = test_vtor /
                    std::sqrt(CGAL::to_double(test_vtor.squared_length()));
                double test_dot_pdt =
                    CGAL::to_double(test_vtor * CH_normals[j]);
                if (test_dot_pdt < 0) {
                    vtx_inside_bool[i] = true;
                    ++total;
                    break;
                }
            }
        }

        if (total == 0) {
            // cell is entirely contained in the included area
            out_cells.push_back(cell_ite);
        } else if (total == 4) {
            // cell is outside the included area
            continue;
        } else {
            // cell instersects the included area
            out_intersecting_cells.push_back(cell_ite);
            intersecting_bool.push_back(vtx_inside_bool);
            intersecting_total.push_back(total);
        }
    }

    return;
}

// Discard parts of cells outside the specified triangulation using
// intersections
double discard_CH_1(NA_Vector const &in_intersecting_cells,
    Triang_Vector const &CH_triangs,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    const std::vector<unsigned int> &intersecting_total,
    Poly_Vector &border_poly,
    std::vector<std::array<double, 3>> &in_vtces_radii,
    unsigned int &atom_cnt_poly) {
    // Use this vector to obtain the indices of the vtces used to form the
    // intersecting segments
    std::vector<unsigned int> v = {0, 1, 2, 3};
    double volume = 0;

    for (unsigned int each = 0; each < in_intersecting_cells.size(); ++each) {

        switch (intersecting_total[each]) {
        case 3: {
            for (unsigned int i = 0; i < 4; i++) {
                if (!intersecting_bool[each][i]) {
                    // Got the handles of the included_area cell(s) that contain
                    // vertices of the intersecting cell
                    // 1 vertex inside included area
                    std::vector<Point> i_points;

                    // Get the points of the vtces of the intersecting cell
                    // p0 is the point inside the included area
                    Point p0(in_intersecting_cells[each]->vertex(i)->point());
                    double const VdW_radius_0 = double(
                        in_intersecting_cells[each]->vertex(i)->info()._radius);
                    Point p1(in_intersecting_cells[each]
                                 ->vertex(get_i_not_equal(v, i, 1))
                                 ->point());
                    Point p2(in_intersecting_cells[each]
                                 ->vertex(get_i_not_equal(v, i, 2))
                                 ->point());
                    Point p3(in_intersecting_cells[each]
                                 ->vertex(get_i_not_equal(v, i, 3))
                                 ->point());
                    // Store the VdW radii
                    std::array<double, 3> temp = {VdW_radius_0, 0, 0};
                    in_vtces_radii.push_back(temp);

                    // Get the intersecting segments
                    std::vector<Segment> segments = {
                        Segment(p0, p1), Segment(p0, p2), Segment(p0, p3)};

                    // Get the intersections
                    for (auto const &s : segments) {
                        for (auto const &t : CH_triangs) {
                            Object inter_obj = CGAL::intersection(s, t);
                            if (Point const *inter_point =
                                    CGAL::object_cast<Point>(&inter_obj)) {
                                i_points.push_back(*inter_point);
                                break;
                            }
                        }
                    }

                    // Check if included area and intersecting cell are sharing
                    // vertices, creating degeneracies.
                    bool degeneracies_bool = false;
                    for (unsigned int i = 0; i < i_points.size() - 1; i++) {
                        for (unsigned int k = i + 1; k < i_points.size(); k++) {
                            if (i_points[i] == i_points[k]) {
                                i_points[k] = i_points[k] +
                                    Vector(0.01 * i, 0.01 * k, 0.01);
                                degeneracies_bool = true;
                            }
                        }
                    }
                    if (degeneracies_bool == true) {
                        p0 = p0 - Vector(0.01, 0.01, 0.01);
                    }

                    // Construct the polyhedron and get its volume
                    Polyhedron CH;
                    CH.make_tetrahedron(
                        p0, i_points[0], i_points[1], i_points[2]);
                    border_poly.push_back(CH);
                    // This polyhedron contains 4 vertices
                    atom_cnt_poly += 4;
                    // Add the volume of this polyhedron
                    volume = volume +
                        std::abs(CGAL::to_double(CGAL::volume(
                            p0, i_points[0], i_points[1], i_points[2])));

                    // Substract the volume of the sphere sector
                    volume = volume -
                        sphere_sector_vol(p0, i_points[0], i_points[1],
                            i_points[2], VdW_radius_0);
                }
            }
            break;
        }
        case 2: {
            for (unsigned int i = 0; i < 4; i++) {
                if (!intersecting_bool[each][i]) {
                    for (unsigned int j = i + 1; j < 4; j++) {
                        if (!intersecting_bool[each][j]) {
                            // 2 vertices inside included area
                            std::vector<Point> i_points;

                            // Get the points of the vtces of the intersecting
                            // cell
                            Point p0(in_intersecting_cells[each]
                                         ->vertex(i)
                                         ->point());
                            Point p1(in_intersecting_cells[each]
                                         ->vertex(j)
                                         ->point());
                            // p0 p1 are the points inside the included area
                            double const VdW_radius_0 =
                                double(in_intersecting_cells[each]
                                           ->vertex(i)
                                           ->info()
                                           ._radius);
                            double const VdW_radius_1 =
                                double(in_intersecting_cells[each]
                                           ->vertex(j)
                                           ->info()
                                           ._radius);
                            std::vector<unsigned int> query_vec = {i, j};
                            Point p2(
                                in_intersecting_cells[each]
                                    ->vertex(get_i_not_equal(v, query_vec, 1))
                                    ->point());
                            Point p3(
                                in_intersecting_cells[each]
                                    ->vertex(get_i_not_equal(v, query_vec, 2))
                                    ->point());
                            // Store the VdW radii
                            std::array<double, 3> temp = {
                                VdW_radius_0, VdW_radius_1, 0};
                            in_vtces_radii.push_back(temp);

                            // Get the intersecting segments
                            std::vector<Segment> segments = {Segment(p0, p2),
                                Segment(p0, p3), Segment(p1, p2),
                                Segment(p1, p3)};

                            // Get the intersections
                            for (auto const &s : segments) {
                                for (auto const &t : CH_triangs) {
                                    Object inter_obj = CGAL::intersection(s, t);
                                    if (Point const *inter_point =
                                            CGAL::object_cast<Point>(
                                                &inter_obj)) {
                                        i_points.push_back(*inter_point);
                                        break;
                                    }
                                }
                            }

                            // Check if included area and intersecting cell are
                            // sharing
                            // vertices, creating degeneracies.
                            bool degeneracies_bool = false;
                            for (unsigned int i = 0; i < i_points.size() - 1;
                                 i++) {
                                for (unsigned int k = i + 1;
                                     k < i_points.size(); k++) {
                                    if (i_points[k] == i_points[k]) {
                                        i_points[k] = i_points[k] +
                                            Vector(0.01, 0.01, 0.01);
                                        degeneracies_bool = true;
                                    }
                                }
                            }
                            if (degeneracies_bool == true) {
                                p0 = p0 - Vector(0.01, 0.01, 0.01);
                                p1 = p1 - Vector(0.01, -0.01, 0.01);
                            }

                            // Construct the polyhedron and get its volume
                            Polyhedron CH;
                            CH.make_tetrahedron(
                                p0, i_points[0], i_points[1], p1);
                            CH.make_tetrahedron(
                                p1, i_points[0], i_points[1], i_points[2]);
                            CH.make_tetrahedron(
                                p1, i_points[1], i_points[2], i_points[3]);
                            border_poly.push_back(CH);
                            // This polyhedron contains 12 vertices
                            atom_cnt_poly += 12;
                            // Add the volume of this polyhedron
                            volume = volume +
                                std::abs(CGAL::to_double(CGAL::volume(
                                    p0, i_points[0], i_points[1], p1)));
                            volume = volume +
                                std::abs(CGAL::to_double(CGAL::volume(p1,
                                    i_points[0], i_points[1], i_points[2])));
                            volume = volume +
                                std::abs(CGAL::to_double(CGAL::volume(p1,
                                    i_points[1], i_points[2], i_points[3])));

                            // Substract the volume of the sphere sector
                            volume = volume -
                                sphere_sector_vol(p0, i_points[0], i_points[1],
                                    p1, VdW_radius_0);
                            volume = volume -
                                sphere_sector_vol(p1, i_points[0], i_points[1],
                                    i_points[2], VdW_radius_1);
                            volume = volume -
                                sphere_sector_vol(p1, i_points[1], i_points[2],
                                    i_points[3], VdW_radius_1);
                        }
                    }
                }
            }
            break;
        }
        case 1: {
            for (unsigned int i = 0; i < 4; i++) {
                if (!intersecting_bool[each][i]) {
                    for (unsigned int j = i + 1; j < 4; j++) {
                        if (!intersecting_bool[each][j]) {
                            for (unsigned int k = j + 1; k < 4; k++) {
                                if (!intersecting_bool[each][k]) {
                                    // 3 vertices inside the included area
                                    std::vector<Point> i_points;

                                    // Get the points of the vtces of the
                                    // intersecting cell
                                    Point p0(in_intersecting_cells[each]
                                                 ->vertex(i)
                                                 ->point());
                                    Point p1(in_intersecting_cells[each]
                                                 ->vertex(j)
                                                 ->point());
                                    Point p2(in_intersecting_cells[each]
                                                 ->vertex(k)
                                                 ->point());
                                    // p0 p1 p2 are the points inside the
                                    // included area
                                    double const VdW_radius_0 =
                                        double(in_intersecting_cells[each]
                                                   ->vertex(i)
                                                   ->info()
                                                   ._radius);
                                    double const VdW_radius_1 =
                                        double(in_intersecting_cells[each]
                                                   ->vertex(j)
                                                   ->info()
                                                   ._radius);
                                    double const VdW_radius_2 =
                                        double(in_intersecting_cells[each]
                                                   ->vertex(k)
                                                   ->info()
                                                   ._radius);
                                    std::vector<unsigned int> query_vec = {
                                        i, j, k};
                                    Point p3(in_intersecting_cells[each]
                                                 ->vertex(get_i_not_equal(
                                                     v, query_vec, 1))
                                                 ->point());
                                    // Store the VdW radii
                                    std::array<double, 3> temp = {VdW_radius_0,
                                        VdW_radius_1, VdW_radius_2};
                                    in_vtces_radii.push_back(temp);

                                    // Get the intersecting segments
                                    std::vector<Segment> segments = {
                                        Segment(p0, p3), Segment(p1, p3),
                                        Segment(p2, p3)};

                                    // Get the intersections
                                    for (auto const &s : segments) {
                                        for (auto const &t : CH_triangs) {
                                            Object inter_obj =
                                                CGAL::intersection(s, t);
                                            if (Point const *inter_point =
                                                    CGAL::object_cast<Point>(
                                                        &inter_obj)) {
                                                i_points.push_back(
                                                    *inter_point);
                                                break;
                                            }
                                        }
                                    }

                                    // Check if included area and intersecting
                                    // cell are sharing
                                    // vertices, creating degeneracies.
                                    bool degeneracies_bool = false;
                                    for (unsigned int i = 0;
                                         i < i_points.size() - 1; i++) {
                                        for (unsigned int k = i + 1;
                                             k < i_points.size(); k++) {
                                            if (i_points[k] == i_points[k]) {
                                                i_points[k] = i_points[k] +
                                                    Vector(0.01, 0.01, 0.01);
                                                degeneracies_bool = true;
                                            }
                                        }
                                    }
                                    if (degeneracies_bool == true) {
                                        p0 = p0 - Vector(0.01, 0.01, 0.01);
                                        p1 = p1 - Vector(0.01, -0.01, 0.01);
                                        p2 = p2 - Vector(0.01, -0.01, -0.01);
                                    }

                                    // Construct the polyhedron and get its
                                    // volume
                                    Polyhedron CH;
                                    CH.make_tetrahedron(
                                        p0, p1, p2, i_points[0]);
                                    CH.make_tetrahedron(
                                        p0, p2, i_points[0], i_points[1]);
                                    CH.make_tetrahedron(p2, i_points[0],
                                        i_points[1], i_points[2]);
                                    border_poly.push_back(CH);
                                    // This polyhedron contains 12 vertices
                                    atom_cnt_poly += 12;
                                    // Add the volume of this polyhedron
                                    volume = volume +
                                        std::abs(CGAL::to_double(CGAL::volume(
                                            p0, p1, p2, i_points[0])));
                                    volume = volume +
                                        std::abs(CGAL::to_double(CGAL::volume(
                                            p0, p2, i_points[0], i_points[1])));
                                    volume = volume +
                                        std::abs(CGAL::to_double(
                                            CGAL::volume(p2, i_points[0],
                                                i_points[1], i_points[2])));

                                    // Substract the volume of the sphere sector
                                    // 1st tetra
                                    volume = volume -
                                        sphere_sector_vol(p0, p1, p2,
                                            i_points[0], VdW_radius_0);
                                    volume = volume -
                                        sphere_sector_vol(p1, p2, p0,
                                            i_points[0], VdW_radius_1);
                                    volume = volume -
                                        sphere_sector_vol(p2, p0, p1,
                                            i_points[0], VdW_radius_2);
                                    // 2nd tetra
                                    volume = volume -
                                        sphere_sector_vol(p0, p2, i_points[0],
                                            i_points[1], VdW_radius_0);
                                    volume = volume -
                                        sphere_sector_vol(p2, p0, i_points[0],
                                            i_points[1], VdW_radius_2);
                                    // 3rd tetra
                                    volume = volume -
                                        sphere_sector_vol(p2, i_points[0],
                                            i_points[1], i_points[2],
                                            VdW_radius_2);
                                }
                            }
                        }
                    }
                }
            }
            break;
        }
        }
    }

    return volume;
}

// Extract vertices coordinates from the cells and store them in an "md_vector".
void na_vector_into_ndd_vector(
    NA_Vector const &in_cells, NDD_Vector &out_cells) {

    for (auto const &each : in_cells) {
        std::array<std::pair<Point, double>, 4> cell = {
            std::make_pair(
                each->vertex(0)->point(), each->vertex(0)->info()._radius),
            std::make_pair(
                each->vertex(1)->point(), each->vertex(1)->info()._radius),
            std::make_pair(
                each->vertex(2)->point(), each->vertex(2)->info()._radius),
            std::make_pair(
                each->vertex(3)->point(), each->vertex(3)->info()._radius)};

        out_cells.push_back(std::move(cell));
    }

    return;
}

// Tool for reading PDB to draw included area.
void tool_PDB_to_CH(
    std::string const &in_filename, std::string const &out_filename) {

    // Read molecule
    chemfiles::Trajectory in_traj(in_filename);
    chemfiles::Frame in_frame = in_traj.read();
    auto in_xyz = in_frame.positions();
    chemfiles::Topology in_top = in_frame.topology();

    std::ofstream outfile(out_filename);
    std::vector<Point> incl_area_points;
    Polyhedron CH;

    for (const chemfiles::Vector3D &atom_xyz : in_xyz) {
        incl_area_points.push_back(
            Point(atom_xyz[0], atom_xyz[1], atom_xyz[2]));
    }

    CGAL::convex_hull_3(incl_area_points.begin(), incl_area_points.end(), CH);
    P_Vertex_iterator v_end = CH.vertices_end();

    if (outfile.is_open()) {
        for (P_Vertex_iterator v_ite = CH.vertices_begin(); v_ite != v_end;
             ++v_ite) {
            outfile << v_ite->point() << '\n';
        }
        outfile.close();
    } else
        throw std::runtime_error("Cannot open output file. Aborting.");

    return;
}

// Tool for normalizing PDB, by renumbering its atoms and residues.
void tool_PDB_norm(
    std::string const &in_filename, std::string const &tool_pdb_norm) {

    std::vector<chemfiles::Residue> res_vec;
    std::vector<unsigned int> resid_vec;

    // Read PDB
    chemfiles::Trajectory input_pdb_traj(in_filename);
    chemfiles::Trajectory output_pdb_traj(tool_pdb_norm, 'w');
    auto nsteps = input_pdb_traj.nsteps();

    for (size_t i = 0; i < nsteps; i++) {
        auto input_pdb_frame = input_pdb_traj.read();
        auto input_pdb_top = input_pdb_frame.topology();
        auto in_xyz = input_pdb_frame.positions();
        chemfiles::Topology output_pdb_top;
        chemfiles::Frame output_pdb_frame;

        // Construct the new molecule.
        size_t j = 1;
        for (auto const &residuo : input_pdb_top.residues()) {
            chemfiles::Residue res(residuo.name(), j);
            for (auto const &k : residuo) {
                // Add atom to residue.
                res.add_atom(k);
                // Add atom to frame and topology.
                chemfiles::Atom atomo(
                    input_pdb_top[k].name(), input_pdb_top[k].type());
                atomo.set("is_hetatm",
                    input_pdb_top[k].get("is_hetatm").value_or(true));
                output_pdb_frame.add_atom(atomo, in_xyz[k]);
                output_pdb_top.add_atom(std::move(atomo));
            }
            output_pdb_top.add_residue(res);
            ++j;
        }
        output_pdb_frame.set_topology(output_pdb_top);

        // Write normalized PDB.
        output_pdb_traj.write(output_pdb_frame);
    }
    return;
}

} // namespace ANA
