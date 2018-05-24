#ifndef ANAWRITE
#define ANAWRITE
#include "ANAincludes.cpp"
#include "ANAutils.hpp"

namespace ANA {
extern std::ofstream out_vol_stream;
// Draw pockets in pymol CGO objects
void draw_raw_cgo(const NA_Matrix& list_of_pockets, const Poly_Vector& polys,
    std::string& out_script_template, const std::string pdb_filename,
    const unsigned int precision);
// Draw 1 pocket in pymol CGO objects
void draw_raw_cgo(std::ofstream& pymol_script, const NA_Vector& cells_to_draw);
// Draw a vector of polyhedrons
void draw_raw_polyhedrons(std::ofstream& pymol_script,
    const Poly_Vector& CH_vec, const std::string& model_template);
// Draw pockets as points in pymol CGO objects
void draw_grid_cgo(const NA_Matrix& list_of_pockets,
    const std::vector<std::array<double, 3>>& in_vtces_radii,
    const std::vector<unsigned int>& intersecting_total,
    const Poly_Vector& polys, std::string& out_script_template,
    const std::string pdb_filename, const double sphere_size,
    const unsigned int sphere_count, const unsigned int precision);
// Draw cells as points in pymol CGO objects
void draw_grid_cgo(std::ofstream& pymol_script, const NA_Vector& cells_to_draw,
    const double sphere_size, const unsigned int sphere_count);
// Draw a vector of polyhedrons as CGO pymol spheres
void draw_grid_cgo_polyhedrons(std::ofstream& pymol_script,
    const std::vector<std::array<double, 3>>& in_vtces_radii,
    const std::vector<unsigned int>& intersecting_total,
    const Poly_Vector& CH_vec, const double sphere_size,
    const double sphere_count);
// Draw pockets in .PDB format. Static version
void draw_raw_PDB(const NA_Vector& list_of_pockets, const Poly_Vector& polys,
    const std::string& out_filename);
// Draw pockets in .PDB format. MD version
void draw_raw_PDB(const NA_Vector& pocket, const Poly_Vector& polys,
    std::string out_filename, const unsigned int frame_cnt);
// Draw pockets in .PDB format. MD version, whole trajectory.
void draw_raw_PDB(const NDD_Matrix& list_of_pockets,
    const Poly_Matrix& list_of_polys, std::string& out_filename,
    const unsigned int max_atom_cnt);
// Draw pockets as dots in a PDB. MD version, whole trajectory.
void draw_grid_pdb(const NDD_Matrix& list_of_pockets,
    const Poly_Matrix& list_of_polys,
    const std::vector<std::array<double, 3>>& in_vtces_radii,
    const std::vector<std::vector<unsigned int>>& intersecting_total,
    const unsigned int sphere_count, const unsigned int precision,
    std::string& out_filename);
// Draw cells as dots in a PDB. Using "NDD_Vector" data structure.
unsigned int make_grid_pdb(const NDD_Vector& cells_to_draw,
    chemfiles::Topology& ana_void_top, chemfiles::Frame& ana_void_frame,
    const unsigned int sphere_count, unsigned int& res_cnt);
// Draw a vector of polyhedrons in a PDB.
unsigned int make_grid_pdb_polyhedrons(chemfiles::Topology& ana_void_top,
    chemfiles::Frame& ana_void_frame,
    const std::vector<std::array<double, 3>>& in_vtces_radii,
    const std::vector<unsigned int>& intersecting_total,
    const Poly_Vector& CH_vec, const double sphere_count, unsigned int atom_cnt,
    unsigned int& res_cnt);
// Draw pockets as dots in a PDB. Using CGAL Cell data
// structure. Static version.
void draw_grid_pdb(const NA_Vector& pocket,
    const std::vector<std::array<double, 3>>& in_vtces_radii,
    const std::vector<unsigned int>& intersecting_total,
    const Poly_Vector& polys, std::string& out_filename,
    const unsigned int sphere_count, const unsigned int precision);
// Draw cells as dots in a PDB. Using CGAL Cell data structure. Static version.
unsigned int make_grid_pdb(const NA_Vector& cells_to_draw,
    chemfiles::Topology& ana_void_top, chemfiles::Frame& ana_void_frame,
    const unsigned int sphere_count, unsigned int& res_cnt);

//////////
//
/////////

// Write pymol script header
void header_script_pymol(
    std::ofstream& pymol_script, const std::string pdb_filename);
// Write .PDB header
void header_PDB(
    const std::string& out_pdb_filename, const std::string& in_pdb_filename);
// Construct a residue object with 4 atoms starting at index cell_cnt
chemfiles::Residue make_cell_residue(
    const unsigned int cell_cnt, const unsigned int atom_cnt);
// Construct a residue object with atoms from 'first_atom' to 'atom_cnt'
inline chemfiles::Residue make_polyhedron_residue(const unsigned int res_cnt,
    const unsigned int first_atom, const unsigned int atom_cnt);
// Construct a residue object for grid output
chemfiles::Residue make_grid_residue(const unsigned int atom_cnt_old,
    const unsigned int atom_cnt, const unsigned int resi);
// Connect a residue object with the atoms starting at index cell_cnt
void connect_cell_residue(
    chemfiles::Topology& topology, const unsigned int atom_index);
// Connect a residue object with the atoms starting at index atom_index
inline void connect_facet(
    chemfiles::Topology& topology, const unsigned int atom_index);

//////////
//
/////////

// Draw a whole convex hull contained in a polyhedron. Static version.
void draw_CH(const Polyhedron& CH, std::string& out_filename);
// Draw a whole convex hull contained in vector of triangles. Static version.
void draw_CH(const Triang_Vector& CH_triang, std::string& out_filename);
// Draw a whole convex hull contained in vector of triangles. MD version.
void draw_CH(const Triang_Vector& CH_triang, chemfiles::Trajectory& out_traj);
// Draw a whole triangulation
void draw_triangulation(const Delaunay& T, std::string& script_filename);
// Write wall amino acids and atoms.
void wall_atom_output(std::ofstream& wall_out, const NA_Vector& in_cells,
    const NA_Vector& in_intersecting_cells,
    const std::vector<std::array<bool, 4>> intersecting_bool,
    const bool requested_CH, const unsigned int precision,
    const unsigned int pock_cnt, const unsigned int frame_cnt,
    const std::string& list_wall_separator);
// Write wall amino acids.
void wall_aa_output(std::ofstream& wall_out, const NA_Vector& in_cells,
    const NA_Vector& in_intersecting_cells,
    const std::vector<std::array<bool, 4>> intersecting_bool,
    const bool requested_CH, const unsigned int precision,
    const unsigned int pock_cnt, const unsigned int frame_cnt,
    const std::string& list_wall_separator);
// Get names and indices of participating atoms and amino acids.
void get_info_cell(const NA_Vector& null_areas_vtor,
    std::vector<unsigned int>& wall_aa_idx,
    std::vector<std::string>& wall_aa_id,
    std::vector<unsigned int>& wall_atom_idx);
// Get names and indices of participating atoms and amino acids for intersecting
// cells.
void get_info_cell(const NA_Vector& cavity_intersecting_cells,
    const std::vector<std::array<bool, 4>>& intersecting_bool,
    std::vector<unsigned int>& wall_aa_idx,
    std::vector<std::string>& wall_aa_id,
    std::vector<unsigned int>& wall_atom_idx);
// Get names and indices of participating amino acids.
void get_info_cell(const NA_Vector& null_areas_vtor,
    std::vector<unsigned int>& wall_aa_idx,
    std::vector<std::string>& wall_aa_id);
// Get names and indices of participating amino acids from intersecting cells.
void get_info_cell(const NA_Vector& cavity_intersecting_cells,
    const std::vector<std::array<bool, 4>>& intersecting_bool,
    std::vector<unsigned int>& wall_aa_idx,
    std::vector<std::string>& wall_aa_id);
// Write wall amino acids and atoms.
void write_wall_file(std::ofstream& pock_out_file,
    const std::string& pock_out_filename,
    const std::vector<unsigned int>& wall_aa_idx,
    const std::vector<std::string>& wall_aa_id,
    const std::vector<unsigned int>& wall_atom_idx,
    const unsigned int frame_cnt, const std::string& list_wall_separator);
// Write wall amino acids.
void write_wall_file(std::ofstream& pock_out_file,
    const std::string& pock_out_filename,
    const std::vector<unsigned int>& wall_aa_idx,
    const std::vector<std::string>& wall_aa_id, const unsigned int frame_cnt,
    const std::string& list_wall_separator);
// Final function to output volume. NA_Matrix (static) version.
void write_output_volume(NA_Matrix const& null_areas_vt_mt,
    const double poly_vol, const std::string& out_vol);
// Final function to output volume. NA_Vector (NDD) version.
void write_output_volume(NA_Vector const& null_areas_vt_mt,
    const double poly_vol, const std::string& out_vol);
// Final function to output volume. NA_Vector (MD) version.
void write_output_volume(NA_Vector const& null_areas_vtor,
    const double poly_vol, const unsigned int frame_cnt,
    const std::string& out_vol);
// Open output volume file, if requested.
void open_vol_file(std::string const& out_vol);
// Final function to output volume. NA_Vector (MD) version.
void write_output_volume(NA_Vector const& null_areas_vtor,
    double const poly_vol, unsigned int const frame_cnt);
}
namespace ANA {
namespace NDD {
// Write file with volumes of each input PDB
inline void ndd_write_out_file(
    const std::vector<double>& output_volumes, const std::string& out_file) {
    std::ofstream output(out_file);
    int i = 1;

    if (output.is_open()) {
        output << "Frame\tVolume" << '\n';
        for (double volume : output_volumes) {
            output << i << "\t" << volume << '\n';
            ++i;
        }
    } else
        throw std::runtime_error("Unable to open output file for NDD");
    return;
}
} // namespace NDD
} // namespace ANA
#endif
