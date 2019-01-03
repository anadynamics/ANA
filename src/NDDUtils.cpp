#include <ANA/NDDUtils.hpp>

namespace ANA {
namespace NDD {

    // Refine the provided list of amino acids. If its not present, then return
    // an array of two '0' elements.
    std::vector<unsigned int> adapt_AA_list(
        std::string &aa_list_proto, unsigned int const top) {
        std::vector<unsigned int> aa_list;

        if (aa_list_proto == "none") {
            aa_list.push_back(0);
            aa_list.push_back(0);
            return aa_list;
        } else {
            // borro el caracter ' pq me caga todo
            std::stringstream stream_aa(aa_list_proto);
            std::string temp_aa;
            while (!stream_aa.eof()) {
                unsigned int aa;
                stream_aa >> temp_aa;
                try {
                    aa = std::stoi(temp_aa);
                } catch (std::invalid_argument const &ia) {
                    // some character present. Doesn't matter, skip it and keep
                    // reading.
                    continue;
                } catch (std::out_of_range const oor) {
                    // int is too large to be represented by int
                    std::cerr << "Invalid residue number: " << temp_aa << '\n';
                    std::exit(1);
                } catch (...) {
                    // some other exception. Doesn't matter, skip it and keep
                    // reading.
                    std::cerr << "Invalid residue number: " << temp_aa
                              << ". Will keep reading." << '\n';
                    continue;
                }
                aa_list.push_back(aa);
            }
            // sort list of included amino acids
            std::sort(aa_list.begin(), aa_list.end());

            if (top != 0) {
                if (aa_list[aa_list.size() - 1] > top) {
                    std::cerr
                        << "Atom / residue list goes out of bounds. Check this "
                           "input list and your input PDB atom count. Quiting "
                           "now."
                        << '\n';
                    exit(0);
                }
            }

            return aa_list;
        }
    }

    // Read coordinates in pdb format using chemfiles.
    void read(const std::string &filename, bool const atom_only,
        std::string &include_CH_AA_proto, std::string &include_CH_atom_proto,
        std::string &sphere_proto, std::string &cylinder_proto,
        std::string &prism_proto, const std::string &include_CH_filename,
        ANA_molecule &molecule_points,
        std::vector<unsigned int> &include_CH_atoms, Triang_Vector &CH_triangs,
        std::vector<unsigned int> &hetatm_atoms) {

        std::vector<Point> incl_area_points;

        // Read PDB
        chemfiles::Trajectory input_pdb_traj(filename);
        chemfiles::Frame const input_pdb_frame = input_pdb_traj.read();
        auto in_xyz = input_pdb_frame.positions();
        chemfiles::Topology const input_pdb_top = input_pdb_frame.topology();
        unsigned int const natoms = input_pdb_top.natoms();
        // Residue selection takes precedence over atom selection, for included
        // area.
        auto const include_CH_AA = adapt_AA_list(include_CH_AA_proto);
        auto const CH_AA_end = include_CH_AA.cend();

        if (include_CH_AA[0] == 0) {
            include_CH_atoms = adapt_AA_list(include_CH_atom_proto, natoms);
        }
        bool listed_incl_CH = false;

        // If requested, use atoms indices to get included area.
        if (include_CH_atoms[0] == 0) {
            for (auto &each : include_CH_atoms) {
                // 0-index normalization.
                each = each - 1;
                incl_area_points.push_back(
                    Point(in_xyz[each][0], in_xyz[each][1], in_xyz[each][2]));
                listed_incl_CH = true;
            }
        }
        // Store the atoms points along with information in the variable
        // molecule_points.
        molecule_points.reserve(natoms);
        for (auto const &residuo : input_pdb_top.residues()) {
            auto res_name = residuo.name();
            for (auto const &i : residuo) {
                if (atom_only &&
                    input_pdb_top[i]
                        .get("is_hetatm")
                        .value_or(false)
                        .as_bool()) {
                    // Save the HETATM indices to discard them during MD and
                    // NDD runs.
                    hetatm_atoms.push_back(i);
                    continue;
                }
                // Iterate over each atom, construct the vertex info and get
                // coords.
                // GetVdwRad or GetCovalentRad?
                Vtx_info vi1;
                vi1._resi = res_name;
                vi1._index = i;
                auto const vdw = input_pdb_top[i].vdw_radius();
                if (vdw) {
                    vi1._radius = vdw.value();
                } else {
                    vi1._radius = 1.5;
                    std::cerr
                        << "Element for atom " << i + 1
                        << " not available. Using Van Der Walls radii of 1.5."
                        << '\n';
                }

                auto resid = residuo.id().value();
                vi1._resn = resid;
                Point p1(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
                molecule_points.push_back(std::make_pair(p1, vi1));

                if (input_pdb_top[i].name() == "CA") {
                    // Use it to draw the included convex hull, if requested.
                    if (include_CH_atoms[0] == 0) {
                        if (std::binary_search(
                                include_CH_AA.begin(), CH_AA_end, resid)) {
                            // The current residue was specified.
                            incl_area_points.push_back(p1);
                            include_CH_atoms.push_back(i);
                        }
                    }
                }
            }
        }

        // Get convex hull if requested.
        if (!(listed_incl_CH || include_CH_atoms[0] == 0) &&
            include_CH_filename != "none") {
            ANA::read_included_area(include_CH_filename, incl_area_points);
            if (incl_area_points.size() < 4) {
                throw std::runtime_error(
                    "Not possible to triangulate less than 4 points.");
            }

            Polyhedron CH;
            CGAL::convex_hull_3(
                incl_area_points.begin(), incl_area_points.end(), CH);
            P_Facet_const_iterator f_end = CH.facets_end();
            for (P_Facet_const_iterator f_ite = CH.facets_begin();
                 f_ite != f_end; ++f_ite) {
                // Fix around the weirdest CGAL bug.
                P_Halfedge_around_facet_const_circulator he_ite =
                    f_ite->facet_begin();
                auto const he_ite_0 = he_ite++;
                auto const he_ite_1 = he_ite++;
                auto const he_ite_2 = he_ite;

                CH_triangs.push_back(Triangle(he_ite_0->vertex()->point(),
                    he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
            }

        } else if (!(listed_incl_CH || (include_CH_atoms[0] == 0)) &&
            sphere_proto != "none") {
            std::stringstream stream_sphere(sphere_proto);
            double x = ANA::parse_double(stream_sphere);
            double y = ANA::parse_double(stream_sphere);
            double z = ANA::parse_double(stream_sphere);
            double r = ANA::parse_double(stream_sphere);
            double const cos_30 = sqrt(3) / 2;
            double const sin_30 = 0.5;

            Point center(x, y, z);
            incl_area_points.push_back(center + Vector(r, 0, 0));
            incl_area_points.push_back(center + Vector(0, r, 0));
            incl_area_points.push_back(center + Vector(0, 0, r));
            incl_area_points.push_back(center + Vector(-r, 0, 0));
            incl_area_points.push_back(center + Vector(0, -r, 0));
            incl_area_points.push_back(center + Vector(0, 0, -r));
            // X-Y plane
            incl_area_points.push_back(
                center + Vector(r * cos_30, r * sin_30, 0));
            incl_area_points.push_back(
                center + Vector(r * sin_30, r * cos_30, 0));
            incl_area_points.push_back(
                center + Vector(r * cos_30, -r * sin_30, 0));
            incl_area_points.push_back(
                center + Vector(r * sin_30, -r * cos_30, 0));
            incl_area_points.push_back(
                center + Vector(-r * cos_30, r * sin_30, 0));
            incl_area_points.push_back(
                center + Vector(-r * sin_30, r * cos_30, 0));
            incl_area_points.push_back(
                center + Vector(-r * cos_30, -r * sin_30, 0));
            incl_area_points.push_back(
                center + Vector(-r * sin_30, -r * cos_30, 0));
            // X-Z plane
            incl_area_points.push_back(
                center + Vector(r * cos_30, 0, r * sin_30));
            incl_area_points.push_back(
                center + Vector(r * sin_30, 0, r * cos_30));
            incl_area_points.push_back(
                center + Vector(r * cos_30, 0, -r * sin_30));
            incl_area_points.push_back(
                center + Vector(r * sin_30, 0, -r * cos_30));
            incl_area_points.push_back(
                center + Vector(-r * cos_30, 0, r * sin_30));
            incl_area_points.push_back(
                center + Vector(-r * sin_30, 0, r * cos_30));
            incl_area_points.push_back(
                center + Vector(-r * cos_30, 0, -r * sin_30));
            incl_area_points.push_back(
                center + Vector(-r * sin_30, 0, -r * cos_30));
            // Y-Z plane
            incl_area_points.push_back(
                center + Vector(0, r * cos_30, r * sin_30));
            incl_area_points.push_back(
                center + Vector(0, r * sin_30, r * cos_30));
            incl_area_points.push_back(
                center + Vector(0, r * cos_30, -r * sin_30));
            incl_area_points.push_back(
                center + Vector(0, r * sin_30, -r * cos_30));
            incl_area_points.push_back(
                center + Vector(0, -r * cos_30, r * sin_30));
            incl_area_points.push_back(
                center + Vector(0, -r * sin_30, r * cos_30));
            incl_area_points.push_back(
                center + Vector(0, -r * cos_30, -r * sin_30));
            incl_area_points.push_back(
                center + Vector(0, -r * sin_30, -r * cos_30));

            Polyhedron CH;
            CGAL::convex_hull_3(
                incl_area_points.begin(), incl_area_points.end(), CH);
            P_Facet_const_iterator f_end = CH.facets_end();
            for (P_Facet_const_iterator f_ite = CH.facets_begin();
                 f_ite != f_end; ++f_ite) {
                // Fix around the weirdest CGAL bug.
                P_Halfedge_around_facet_const_circulator he_ite =
                    f_ite->facet_begin();
                auto const he_ite_0 = he_ite++;
                auto const he_ite_1 = he_ite++;
                auto const he_ite_2 = he_ite;

                CH_triangs.push_back(Triangle(he_ite_0->vertex()->point(),
                    he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
            }

        } else if (!(listed_incl_CH || (include_CH_atoms[0] == 0)) &&
            cylinder_proto != "none") {
            std::stringstream stream_cylinder(cylinder_proto);

            double x1 = ANA::parse_double(stream_cylinder);
            double y1 = ANA::parse_double(stream_cylinder);
            double z1 = ANA::parse_double(stream_cylinder);
            double x2 = ANA::parse_double(stream_cylinder);
            double y2 = ANA::parse_double(stream_cylinder);
            double z2 = ANA::parse_double(stream_cylinder);
            double r = ANA::parse_double(stream_cylinder);

            double const cos_30 = sqrt(3) / 2;
            double const sin_30 = 0.5;

            Point center_1(x1, y1, z1);
            Point center_2(x2, y2, z2);
            Vector vdiff(center_2 - center_1);
            Vector n1(-vdiff.y(), vdiff.x(), 0);
            Vector n2 = CGAL::cross_product(vdiff, n1);
            n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
            n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

            // tap 1
            incl_area_points.push_back(center_1 + r * n1);
            incl_area_points.push_back(center_1 + r * n2);
            incl_area_points.push_back(center_1 - r * n1);
            incl_area_points.push_back(center_1 - r * n2);
            incl_area_points.push_back(
                center_1 + r * cos_30 * n1 + r * sin_30 * n2);
            incl_area_points.push_back(
                center_1 + r * sin_30 * n1 + r * cos_30 * n2);
            incl_area_points.push_back(
                center_1 + r * cos_30 * n1 - r * sin_30 * n2);
            incl_area_points.push_back(
                center_1 + r * sin_30 * n1 - r * cos_30 * n2);
            incl_area_points.push_back(
                center_1 - r * cos_30 * n1 + r * sin_30 * n2);
            incl_area_points.push_back(
                center_1 - r * sin_30 * n1 + r * cos_30 * n2);
            incl_area_points.push_back(
                center_1 - r * cos_30 * n1 - r * sin_30 * n2);
            incl_area_points.push_back(
                center_1 - r * sin_30 * n1 - r * cos_30 * n2);

            // tap 2
            incl_area_points.push_back(center_2 + r * n1);
            incl_area_points.push_back(center_2 + r * n2);
            incl_area_points.push_back(center_2 - r * n1);
            incl_area_points.push_back(center_2 - r * n2);
            incl_area_points.push_back(
                center_2 + r * cos_30 * n1 + r * sin_30 * n2);
            incl_area_points.push_back(
                center_2 + r * sin_30 * n1 + r * cos_30 * n2);
            incl_area_points.push_back(
                center_2 + r * cos_30 * n1 - r * sin_30 * n2);
            incl_area_points.push_back(
                center_2 + r * sin_30 * n1 - r * cos_30 * n2);
            incl_area_points.push_back(
                center_2 - r * cos_30 * n1 + r * sin_30 * n2);
            incl_area_points.push_back(
                center_2 - r * sin_30 * n1 + r * cos_30 * n2);
            incl_area_points.push_back(
                center_2 - r * cos_30 * n1 - r * sin_30 * n2);
            incl_area_points.push_back(
                center_2 - r * sin_30 * n1 - r * cos_30 * n2);

            Polyhedron CH;
            CGAL::convex_hull_3(
                incl_area_points.begin(), incl_area_points.end(), CH);
            P_Facet_const_iterator f_end = CH.facets_end();
            for (P_Facet_const_iterator f_ite = CH.facets_begin();
                 f_ite != f_end; ++f_ite) {
                // Fix around the weirdest CGAL bug.
                P_Halfedge_around_facet_const_circulator he_ite =
                    f_ite->facet_begin();
                auto const he_ite_0 = he_ite++;
                auto const he_ite_1 = he_ite++;
                auto const he_ite_2 = he_ite;

                CH_triangs.push_back(Triangle(he_ite_0->vertex()->point(),
                    he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
            }

        } else if (!(listed_incl_CH || (include_CH_atoms[0] == 0)) &&
            prism_proto != "none") {
            std::stringstream stream_prism(prism_proto);
            double x1 = ANA::parse_double(stream_prism);
            double y1 = ANA::parse_double(stream_prism);
            double z1 = ANA::parse_double(stream_prism);
            double x2 = ANA::parse_double(stream_prism);
            double y2 = ANA::parse_double(stream_prism);
            double z2 = ANA::parse_double(stream_prism);
            double width = ANA::parse_double(stream_prism) / 2;
            double height = ANA::parse_double(stream_prism) / 2;

            Point center_1(x1, y1, z1);
            Point center_2(x2, y2, z2);
            Vector vdiff(center_2 - center_1);
            Vector n1(-vdiff.y(), vdiff.x(), 0);
            Vector n2 = CGAL::cross_product(vdiff, n1);
            n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
            n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

            // tap 1
            incl_area_points.push_back(center_1 + width * n1 + height * n2);
            incl_area_points.push_back(center_1 + width * n1 - height * n2);
            incl_area_points.push_back(center_1 - width * n1 + height * n2);
            incl_area_points.push_back(center_1 - width * n1 - height * n2);

            // tap 2
            incl_area_points.push_back(center_2 + width * n1 + height * n2);
            incl_area_points.push_back(center_2 + width * n1 - height * n2);
            incl_area_points.push_back(center_2 - width * n1 + height * n2);
            incl_area_points.push_back(center_2 - width * n1 - height * n2);

            Polyhedron CH;
            CGAL::convex_hull_3(
                incl_area_points.begin(), incl_area_points.end(), CH);
            P_Facet_const_iterator f_end = CH.facets_end();
            for (P_Facet_const_iterator f_ite = CH.facets_begin();
                 f_ite != f_end; ++f_ite) {
                // Fix around the weirdest CGAL bug.
                P_Halfedge_around_facet_const_circulator he_ite =
                    f_ite->facet_begin();
                auto const he_ite_0 = he_ite++;
                auto const he_ite_1 = he_ite++;
                auto const he_ite_2 = he_ite;

                CH_triangs.push_back(Triangle(he_ite_0->vertex()->point(),
                    he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
            }

        } else if (listed_incl_CH || (include_CH_atoms[0] == 0)) {
            if (incl_area_points.size() < 4) {
                throw std::runtime_error(
                    "Not possible to triangulate less than 4 points.");
            }

            Polyhedron CH;
            CGAL::convex_hull_3(
                incl_area_points.begin(), incl_area_points.end(), CH);
            P_Facet_const_iterator f_end = CH.facets_end();
            for (P_Facet_const_iterator f_ite = CH.facets_begin();
                 f_ite != f_end; ++f_ite) {

                // Fix around the weirdest CGAL bug.
                P_Halfedge_around_facet_const_circulator he_ite =
                    f_ite->facet_begin();
                auto const he_ite_0 = he_ite++;
                auto const he_ite_1 = he_ite++;
                auto const he_ite_2 = he_ite;

                CH_triangs.push_back(Triangle(he_ite_0->vertex()->point(),
                    he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
            }
        }
        return;
    }

    // NDD Specific function for PDB input
    void ndd_read_PDB_get_cells(const std::string &filename,
        NDD_IVector const &in_void_cells_indices, NDD_Vector &output_cells) {

        // Read molecule
        chemfiles::Trajectory in_traj(filename);
        chemfiles::Frame in_frame = in_traj.read();
        auto in_xyz = in_frame.positions();
        chemfiles::Topology in_top = in_frame.topology();
        Vtx_info vi1;
        NDD_Element temp_cell_coords;

        for (NDD_IElement const &idx : in_void_cells_indices) {
            // Iterate over each cell of the current pocket
            for (unsigned int i = 0; i < 4; ++i) {
                // Iterate over each vertex of the current cell
                Point const p0 = Point(
                    in_xyz[idx[i]][0], in_xyz[idx[i]][1], in_xyz[idx[i]][2]);
                temp_cell_coords[i] = std::make_pair(
                    p0, in_top[idx[i]].vdw_radius().value_or(1.5));
            }
            output_cells.push_back(temp_cell_coords);
        }

        return;
    }

    // NDD Specific function for PDB input. Hi precision method
    void ndd_read_PDB_get_cells(const std::string &filename,
        NDD_IVector const &in_void_cells_indices,
        const std::vector<unsigned int> &include_CH_atoms,
        NDD_Vector &output_cells, Triang_Vector &CH_triangs) {

        // Read molecule
        chemfiles::Trajectory in_traj(filename);
        chemfiles::Frame in_frame = in_traj.read();
        auto in_xyz = in_frame.positions();
        chemfiles::Topology in_top = in_frame.topology();
        Vtx_info vi1;
        NDD_Element temp_cell_coords;
        std::vector<Point> incl_area_points;

        for (NDD_IElement const &idx : in_void_cells_indices) {
            // Iterate over each cell of the current pocket
            for (unsigned int i = 0; i < 4; ++i) {
                // Iterate over each vertex of the current cell
                Point const p0 = Point(
                    in_xyz[idx[i]][0], in_xyz[idx[i]][1], in_xyz[idx[i]][2]);
                temp_cell_coords[i] = std::make_pair(
                    p0, in_top[idx[i]].vdw_radius().value_or(1.5));
            }
            output_cells.push_back(temp_cell_coords);
        }

        // Actualize convex hull of the include area.
        for (auto const &each : include_CH_atoms) {
            incl_area_points.push_back(
                Point(in_xyz[each][0], in_xyz[each][1], in_xyz[each][2]));
        }

        Polyhedron CH;
        CGAL::convex_hull_3(
            incl_area_points.begin(), incl_area_points.end(), CH);
        P_Facet_const_iterator f_end = CH.facets_end();
        for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
             ++f_ite) {
            // Fix around the weirdest CGAL bug.
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            CH_triangs.push_back(Triangle(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
        }

        return;
    }

    // Analytical NDD.
    void ndd(NA_Vector const &cavity_void_cells,
        std::string const &modes_ndd_filename, std::string const &out_file) {

        std::vector<double> output_volumes;
        std::vector<unsigned int> all_indices;

        NDD_IVector const cells_indices = get_vertices(cavity_void_cells);

        Modes const modes(modes_ndd_filename);

        // ANA::NDD::ndd_write_out_file(output_volumes, out_file);

        [[maybe_unused]] auto borrame = out_file.end();
        return;
    }

    // Perform Non Delaunay Dynamics.
    void ndd_nondelaunay_dynamics_old(NA_Vector const &cavity_void_cells,
        const std::string &pdb_list, bool const precision,
        const std::vector<unsigned int> include_CH_atoms,
        const std::string &out_file) {
        std::vector<double> output_volumes;
        std::vector<unsigned int> all_indices;

        NDD_IVector cells_indices;
        std::ifstream input_pdbs_filename(pdb_list);
        std::string pdb_filename;

        if (input_pdbs_filename.is_open()) {
            // Get the atoms involved in each pocket.
            // new_cells_coordinates and cells_indices have similar structure.
            // They are vectors of vectors of arrays of 4 elements. Each array
            // refers
            // to a cell. Each vector of arrays refer to a pocket and all of
            // these pockets (vectors of arrays) are compiled into a vector.

            auto cells_indices = get_vertices(cavity_void_cells);

            if (precision == 1) {
                while (std::getline(input_pdbs_filename, pdb_filename)) {
                    NDD_Vector init_cells_coordinates, new_cells_coordinates,
                        cavity_void_coords, cavity_intersecting_coords;
                    Triang_Vector CH_triangs;
                    std::vector<std::array<bool, 4>> intersecting_bool;
                    std::vector<unsigned int> intersecting_total;
                    Poly_Vector border_poly;

                    // Get cell coordinates and the new include area triangles
                    ANA::NDD::ndd_read_PDB_get_cells(pdb_filename,
                        cells_indices, include_CH_atoms, new_cells_coordinates,
                        CH_triangs);
                    // Separate fully inside cells and the intersecting ones.
                    // Get the volume of this last group.
                    ndd_discard_CH_0(new_cells_coordinates, CH_triangs,
                        cavity_void_coords, cavity_intersecting_coords,
                        intersecting_bool, intersecting_total);
                    double poly_vol =
                        ndd_discard_CH_1(cavity_intersecting_coords, CH_triangs,
                            intersecting_bool, intersecting_total, border_poly);
                    // Store result
                    output_volumes.push_back(
                        poly_vol + ndd_get_void_volume(cavity_void_coords));
                }
            } else {
                while (std::getline(input_pdbs_filename, pdb_filename)) {
                    NDD_Vector init_cells_coordinates, new_cells_coordinates,
                        cavity_void_coords;
                    // Get cell coordinates
                    ANA::NDD::ndd_read_PDB_get_cells(
                        pdb_filename, cells_indices, new_cells_coordinates);
                    // Store result
                    output_volumes.push_back(
                        ndd_get_void_volume(new_cells_coordinates));
                }
            }

            ANA::NDD::ndd_write_out_file(output_volumes, out_file);
        } else
            std::cerr << "Unable to open " << pdb_list << " for NDD" << '\n';

        return;
    }

    // Get the indices of the atoms involved in the given cells
    NDD_IVector get_vertices(NA_Vector const &cavity_void_cells) {

        NDD_IVector cells_indices;
        cells_indices.reserve(cavity_void_cells.size() * 4);

        for (Finite_cells_iterator const &ac_ite : cavity_void_cells) {
            NDD_IElement temp{ac_ite->vertex(0)->info()._index,
                ac_ite->vertex(1)->info()._index,
                ac_ite->vertex(2)->info()._index,
                ac_ite->vertex(3)->info()._index};
            cells_indices.push_back(std::move(temp));
        }

        return cells_indices;
    }

    // Calc volume of the input cells. Reedited for array container.
    double ndd_get_void_volume(NDD_Vector const &cavity_void_cells) {

        double current_cell_vol, volume = 0;

        for (const NDD_Element &cell : cavity_void_cells) {
            // Iterate over each cell of and get the total volume of the
            // tetrahedron
            current_cell_vol = CGAL::to_double(CGAL::volume(
                cell[0].first, cell[1].first, cell[2].first, cell[3].first));
            // Substract the volume filled by the 4 atoms in the vtces
            current_cell_vol = ndd_refine_cell_volume(current_cell_vol, cell);
            volume = volume + current_cell_vol;
        }

        return volume;
    }
    // Discard cells without a vertex inside the specified convex hull. Hi
    // precision. NDD version
    void ndd_discard_CH_0(NDD_Vector const &in_coords,
        Triang_Vector const &CH_triangs, NDD_Vector &out_coords,
        NDD_Vector &out_intersecting_coords,
        std::vector<std::array<bool, 4>> &intersecting_bool,
        std::vector<unsigned int> &intersecting_total) {

        std::vector<Vector> CH_normals;
        std::vector<Point> CH_vtces;
        // Triangle normals point inwards. Only inside points will give a
        // positive dot product against all normals
        for (auto const &triangle : CH_triangs) {
            Vector v1 = triangle.vertex(1) - triangle.vertex(0);
            Vector v2 = triangle.vertex(2) - triangle.vertex(1);
            Vector normal = CGAL::cross_product(v2, v1);
            normal =
                normal / std::sqrt(CGAL::to_double(normal.squared_length()));
            CH_normals.push_back(normal);
            CH_vtces.push_back(triangle.vertex(1));
        }

        // Now discard outside cells
        for (auto const &ndd_array : in_coords) {
            std::array<bool, 4> vtx_inside_bool = {false, false, false, false};
            unsigned int total = 0;
            // Get nbr of vertices that lie outside the include area.
            for (unsigned int i = 0; i <= 3; ++i) {
                Point test_point(ndd_array[i].first);

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
                out_coords.push_back(ndd_array);
            } else if (total == 4) {
                // cell is outside the included area
                continue;
            } else {
                // cell instersects the included area
                out_intersecting_coords.push_back(ndd_array);
                intersecting_bool.push_back(vtx_inside_bool);
                intersecting_total.push_back(total);
            }
        }

        return;
    }

    // Discard parts of cells outside the specified triangulation using
    // intersecitons. NDD version
    double ndd_discard_CH_1(NDD_Vector const &in_intersecting_coords,
        Triang_Vector const &CH_triangs,
        const std::vector<std::array<bool, 4>> &intersecting_bool,
        const std::vector<unsigned int> &intersecting_total,
        Poly_Vector &border_poly) {
        // Use this vector to obtain the indices of the vtces used to form the
        // intersecting segments
        std::vector<unsigned int> v = {0, 1, 2, 3};
        double volume = 0;

        for (unsigned int each = 0; each < in_intersecting_coords.size();
             ++each) {

            switch (intersecting_total[each]) {
            case 3: {
                for (unsigned int i = 0; i < 4; i++) {
                    if (!intersecting_bool[each][i]) {
                        // Got the handles of the included_area cell(s) that
                        // contain vertices of the intersecting cell 1 vertex
                        // inside included area
                        std::vector<Point> i_points;

                        // Get the points of the vtces of the intersecting cell
                        // p0 is the point inside the included area
                        Point p0(in_intersecting_coords[each][i].first);
                        double const VdW_radius_0 =
                            double(in_intersecting_coords[each][i].second);
                        Point p1(in_intersecting_coords[each][get_i_not_equal(
                                                                  v, i, 1)]
                                     .first);
                        Point p2(in_intersecting_coords[each][get_i_not_equal(
                                                                  v, i, 2)]
                                     .first);
                        Point p3(in_intersecting_coords[each][get_i_not_equal(
                                                                  v, i, 3)]
                                     .first);

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

                        // Check if included area and intersecting cell are
                        // sharing vertices, creating degeneracies.
                        bool degeneracies_bool = false;
                        for (unsigned int i = 0; i < i_points.size() - 1; i++) {
                            for (unsigned int k = i + 1; k < i_points.size();
                                 k++) {
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

                        // Construct the polyhedron and get its volume.
                        Polyhedron CH;
                        CH.make_tetrahedron(
                            p0, i_points[0], i_points[1], i_points[2]);
                        border_poly.push_back(CH);
                        // Add the volume of this polyhedron.
                        volume = volume +
                            std::abs(CGAL::to_double(CGAL::volume(
                                p0, i_points[0], i_points[1], i_points[2])));

                        // Substract the volume of the sphere sector.
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
                                // Get the points of the vtces of the
                                // intersecting cell
                                Point p0(in_intersecting_coords[each][i].first);
                                Point p1(in_intersecting_coords[each][j].first);
                                // p0 p1 are the points inside the included area
                                double const VdW_radius_0 = double(
                                    in_intersecting_coords[each][i].second);
                                double const VdW_radius_1 = double(
                                    in_intersecting_coords[each][j].second);
                                std::vector<unsigned int> query_vec = {i, j};
                                Point p2(
                                    in_intersecting_coords[each]
                                                          [get_i_not_equal(
                                                               v, query_vec, 1)]
                                                              .first);
                                Point p3(
                                    in_intersecting_coords[each]
                                                          [get_i_not_equal(
                                                               v, query_vec, 2)]
                                                              .first);

                                // Get the intersecting segments
                                std::vector<Segment> segments = {
                                    Segment(p0, p2), Segment(p0, p3),
                                    Segment(p1, p2), Segment(p1, p3)};

                                // Get the intersections
                                for (auto const &s : segments) {
                                    for (auto const &t : CH_triangs) {
                                        Object inter_obj =
                                            CGAL::intersection(s, t);
                                        if (Point const *inter_point =
                                                CGAL::object_cast<Point>(
                                                    &inter_obj)) {
                                            i_points.push_back(*inter_point);
                                            break;
                                        }
                                    }
                                }

                                // Check if included area and intersecting cell
                                // are sharing vertices, creating degeneracies.
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
                                // Add the volume of this polyhedron
                                volume = volume +
                                    std::abs(CGAL::to_double(CGAL::volume(
                                        p0, i_points[0], i_points[1], p1)));
                                volume = volume +
                                    std::abs(CGAL::to_double(
                                        CGAL::volume(p1, i_points[0],
                                            i_points[1], i_points[2])));
                                volume = volume +
                                    std::abs(CGAL::to_double(
                                        CGAL::volume(p1, i_points[1],
                                            i_points[2], i_points[3])));

                                // Substract the volume of the sphere sector
                                volume = volume -
                                    sphere_sector_vol(p0, i_points[0],
                                        i_points[1], p1, VdW_radius_0);
                                volume = volume -
                                    sphere_sector_vol(p1, i_points[0],
                                        i_points[1], i_points[2], VdW_radius_1);
                                volume = volume -
                                    sphere_sector_vol(p1, i_points[1],
                                        i_points[2], i_points[3], VdW_radius_1);
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
                                        Point p0(in_intersecting_coords[each][i]
                                                     .first);
                                        Point p1(in_intersecting_coords[each][j]
                                                     .first);
                                        Point p2(in_intersecting_coords[each][k]
                                                     .first);
                                        // p0 p1 p3 are the points inside the
                                        // included area
                                        double const VdW_radius_0 = double(
                                            in_intersecting_coords[each][i]
                                                .second);
                                        double const VdW_radius_1 = double(
                                            in_intersecting_coords[each][j]
                                                .second);
                                        double const VdW_radius_2 = double(
                                            in_intersecting_coords[each][k]
                                                .second);
                                        std::vector<unsigned int> query_vec = {
                                            i, j, k};
                                        Point p3(in_intersecting_coords
                                                     [each][get_i_not_equal(v,
                                                                query_vec, 1)]
                                                         .first);

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
                                                        CGAL::object_cast<
                                                            Point>(
                                                            &inter_obj)) {
                                                    i_points.push_back(
                                                        *inter_point);
                                                    break;
                                                }
                                            }
                                        }

                                        // Check if included area and
                                        // intersecting cell are sharing
                                        // vertices, creating degeneracies.
                                        bool degeneracies_bool = false;
                                        for (unsigned int i = 0;
                                             i < i_points.size() - 1; i++) {
                                            for (unsigned int k = i + 1;
                                                 k < i_points.size(); k++) {
                                                if (i_points[k] ==
                                                    i_points[k]) {
                                                    i_points[k] = i_points[k] +
                                                        Vector(
                                                            0.01, 0.01, 0.01);
                                                    degeneracies_bool = true;
                                                }
                                            }
                                        }
                                        if (degeneracies_bool == true) {
                                            p0 = p0 - Vector(0.01, 0.01, 0.01);
                                            p1 = p1 - Vector(0.01, -0.01, 0.01);
                                            p2 =
                                                p2 - Vector(0.01, -0.01, -0.01);
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
                                        // Add the volume of this polyhedron
                                        volume = volume +
                                            std::abs(
                                                CGAL::to_double(CGAL::volume(
                                                    p0, p1, p2, i_points[0])));
                                        volume = volume +
                                            std::abs(CGAL::to_double(
                                                CGAL::volume(p0, p2,
                                                    i_points[0], i_points[1])));
                                        volume = volume +
                                            std::abs(CGAL::to_double(
                                                CGAL::volume(p2, i_points[0],
                                                    i_points[1], i_points[2])));

                                        // Substract the volume of the sphere
                                        // sector 1st tetra
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
                                            sphere_sector_vol(p0, p2,
                                                i_points[0], i_points[1],
                                                VdW_radius_0);
                                        volume = volume -
                                            sphere_sector_vol(p2, p0,
                                                i_points[0], i_points[1],
                                                VdW_radius_2);
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

} // namespace NDD
} // namespace ANA