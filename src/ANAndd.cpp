#include "ANAndd.hpp"
namespace ANA {
int NDD_ANA(const std::string &in_filename, std::string &AA_indices_proto,
            const std::string &only_side_ASA, const std::string &ASA_method,
            std::string &exclude_ca_for_ASA_indices_proto,
            const std::string &list_wall,
            const std::string &list_wall_separator,
            std::string &include_CH_aa_proto,
            std::string &include_CH_atom_proto, std::string &sphere_proto,
            std::string &cylinder_proto, std::string &prism_proto,
            const std::string &include_CH_filename,
            const std::string &pdbs_list_ndd_filename,
            const std::string &out_ndd_filename, std::string &out_filename,
            const std::string &out_type,
            const bool triangulate_only_included_aas, const bool atom_only,
            const double minVR, const double maxSR, const double max_probe,
            const double max_probe_length, const double sphere_size,
            const unsigned int sphere_count,
            const unsigned int nbr_of_vertices_to_include,
            const unsigned int precision) {

  // atom_cnt_poly is for MD only.
  unsigned int atom_cnt_poly = 0;
  double poly_vol = 0;
  std::vector<Point> CAs_Points;
  ANA_molecule molecule_points;
  std::vector<unsigned int> AA_indices, CA_indices, include_CH_atoms;
  Point cm;
  Triang_Vector CH_triangs;
  NA_Vector cavity_cells, cavity_included_cells, cavity_void_cells,
      cavity_intersecting_cells, cavity_joint_cells;
  NA_Matrix null_areas_mtx;
  Poly_Vector border_poly;
  std::vector<std::array<bool, 4>> intersecting_bool;
  std::vector<unsigned int> intersecting_total;
  std::vector<unsigned int> wall_aa_idx, wall_atom_idx;
  std::vector<unsigned int> hetatm_atoms;
  std::vector<std::string> wall_aa_id;
  std::vector<std::array<double, 3>> in_vtces_radii;

  // Read original file
  const bool requested_CH = ANA::read_static(
      in_filename, triangulate_only_included_aas, atom_only, AA_indices_proto,
      exclude_ca_for_ASA_indices_proto, include_CH_aa_proto,
      include_CH_atom_proto, sphere_proto, cylinder_proto, prism_proto,
      include_CH_filename, molecule_points, cm, AA_indices, CA_indices,
      CAs_Points, include_CH_atoms, CH_triangs, hetatm_atoms);

  Delaunay T = ANA::triangulate(molecule_points);

  ANA::get_all_voids(T, cavity_cells, minVR, maxSR);

  if (!requested_CH) {

    if (triangulate_only_included_aas == false && AA_indices[0] != 0) {
      ANA::keep_included_aa_cells(cavity_cells, AA_indices,
                                  nbr_of_vertices_to_include,
                                  cavity_included_cells);
    } else {
      cavity_included_cells = cavity_cells;
    }

    if (ASA_method == "none") {
      cavity_void_cells = cavity_included_cells;
    } else if (ASA_method == "cm") {
      ANA::discard_ASA_dot_pdt_cm(cm, CAs_Points, max_probe, max_probe_length,
                                  only_side_ASA, cavity_included_cells,
                                  cavity_void_cells);
    } else if (ASA_method == "backbone") {
      ANA::discard_ASA_CACH(CAs_Points, only_side_ASA, cavity_included_cells,
                            cavity_void_cells);
    } else if (ASA_method == "axes") {
      ANA::discard_ASA_dot_pdt_axes(CAs_Points, max_probe, max_probe_length,
                                    only_side_ASA, cavity_included_cells,
                                    cavity_void_cells);
    }

  } else {
    if (precision == 1) {
      discard_CH_0(cavity_cells, CH_triangs, cavity_void_cells,
                   cavity_intersecting_cells, intersecting_bool,
                   intersecting_total);
      poly_vol = discard_CH_1(cavity_intersecting_cells, CH_triangs,
                              intersecting_bool, intersecting_total,
                              border_poly, in_vtces_radii, atom_cnt_poly);
    } else { // assume precision = 0
      discard_CH_0(cavity_cells, CH_triangs, cavity_void_cells);
    }
    null_areas_mtx.push_back(cavity_void_cells);
  }

  if (out_filename != "none") {
    if (out_type == "raw_pdb") {
      out_filename.append(".pdb");
      ANA::draw_raw_PDB(null_areas_mtx[0], border_poly, out_filename);
    } else if (out_type == "raw_cgo") {
      ANA::draw_raw_cgo(null_areas_mtx, border_poly, out_filename, in_filename,
                        precision);
    } else if (out_type == "grid_cgo") {
      ANA::draw_grid_cgo(null_areas_mtx, in_vtces_radii, intersecting_total,
                         border_poly, out_filename, in_filename, sphere_size,
                         sphere_count, precision);
    } else if (out_type == "grid_pdb") {
      out_filename.append(".pdb");
      ANA::draw_grid_pdb(null_areas_mtx[0], in_vtces_radii, intersecting_total,
                         border_poly, out_filename, sphere_count, precision);
    }
  }

  cavity_joint_cells.reserve(cavity_void_cells.size() +
                             cavity_intersecting_cells.size());
  cavity_joint_cells.insert(cavity_joint_cells.end(), cavity_void_cells.begin(),
                            cavity_void_cells.end());
  cavity_joint_cells.insert(cavity_joint_cells.end(),
                            cavity_intersecting_cells.begin(),
                            cavity_intersecting_cells.end());
  ANA::NDD::ndd_nondelaunay_dynamics(cavity_joint_cells, pdbs_list_ndd_filename,
                                     precision, include_CH_atoms,
                                     out_ndd_filename);

  if (list_wall == "atom") {
    unsigned int pock_cnt = 1;
    // Remove ".pdb"
    std::string filename = in_filename.substr(0, in_filename.size() - 4);
    filename.insert(0, "wall_");
    std::ofstream wall_out(filename);
    for (const NA_Vector &null_areas_vtor : null_areas_mtx) {
      ANA::wall_atom_output(wall_out, null_areas_vtor,
                            cavity_intersecting_cells, intersecting_bool,
                            requested_CH, precision, pock_cnt, 1,
                            list_wall_separator);
      ++pock_cnt;
    }
  } else if (list_wall == "residue") {
    unsigned int pock_cnt = 1;
    // Remove ".pdb"
    std::string filename = in_filename.substr(0, in_filename.size() - 4);
    filename.insert(0, "wall_");
    std::ofstream wall_out(filename);
    for (const NA_Vector &null_areas_vtor : null_areas_mtx) {
      ANA::wall_aa_output(wall_out, null_areas_vtor, cavity_intersecting_cells,
                          intersecting_bool, requested_CH, precision, pock_cnt,
                          1, list_wall_separator);
      ++pock_cnt;
    }
  }

  double volume = ANA::get_void_volume(cavity_void_cells);
  std::cout << volume + poly_vol << '\n';

  return 0;
}
} // namespace ANA
