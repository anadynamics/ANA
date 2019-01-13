#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(std::string const &in_filename, IncludedAreaOptions &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    // atom_cnt_poly is for MD only.
    unsigned int atom_cnt_poly = 0;
    double poly_vol = 0.;
    std::vector<Point> CAs_Points;
    ANA_molecule molecule_points;
    std::vector<unsigned int> CA_indices, include_CH_atoms;
    Triang_Vector CH_triangs;
    NA_Vector cavity_cells, cavity_included_cells, cavity_void_cells,
        cavity_intersecting_cells, cavity_joint_cells;
    NA_Matrix null_areas_mtx;
    Poly_Vector border_poly;
    std::vector<std::array<bool, 4>> intersecting_bool;
    std::vector<unsigned int> intersecting_total;
    std::vector<std::array<double, 3>> in_vtces_radii;

    Molecule const protein = ANA::Molecule(in_filename, atom_only);
    ANA::ConvexHull const CH(protein, IA_opts);
    ANA::Cavity hueco(protein, cell_opts);

    // Delaunay const T = ANA::triangulate(molecule_points);
    // ANA::get_all_voids(T, cavity_cells, cell_opts);

    discard_CH_0(cavity_cells, CH_triangs, cavity_void_cells,
        cavity_intersecting_cells, intersecting_bool, intersecting_total);
    poly_vol =
        discard_CH_1(cavity_intersecting_cells, CH_triangs, intersecting_bool,
            intersecting_total, border_poly, in_vtces_radii, atom_cnt_poly);

    null_areas_mtx.push_back(cavity_void_cells);

    printf("post discard \n");

    cavity_joint_cells.reserve(
        cavity_void_cells.size() + cavity_intersecting_cells.size());
    cavity_joint_cells.insert(cavity_joint_cells.end(),
        cavity_void_cells.begin(), cavity_void_cells.end());
    cavity_joint_cells.insert(cavity_joint_cells.end(),
        cavity_intersecting_cells.begin(), cavity_intersecting_cells.end());

    printf("post reserve\n");

    // ANA::NDD::ndd_nondelaunay_dynamics_old(cavity_joint_cells,
    //     pdbs_list_ndd_filename, precision, include_CH_atoms,
    //     out_ndd_filename);

    ANA::NDD::ndd(cavity_joint_cells, NDD_opts);

    printf("post ndds\n");

    ANA::write_output_volume(cavity_void_cells, poly_vol);

    return 0;
}

} // namespace ANA
