#ifndef ANA_READ_H
#define ANA_READ_H

#include <ANA/Includes.hpp>

namespace ANA {
namespace NDD {
    // NDD Specific function for PDB input
    void ndd_read_PDB_get_cells(std::string const &filename,
        NDD_IVector const &in_cells_indices, NDD_Vector &output_cells);

    // NDD Specific function for PDB input. Hi precision method
    void ndd_read_PDB_get_cells(std::string const &filename,
        NDD_IVector const &in_void_cells_indices,
        const std::vector<unsigned int> &include_CH_atoms,
        NDD_Vector &output_cells, Triang_Vector &CH_triangs);
}
// Refine the provided list of amino acids. If its not present, then return an
// array of one 0 element.
template <class tipo>
bool adapt_aa_list(std::string &aa_list_proto, std::vector<tipo> &aa_list);
// Read coordinates in pdb format using chemfiles.
bool read_static(std::string const &filename,
    bool const triangulate_only_included_aas, bool const atom_only,
    std::string &aa_list_proto, std::string &exclude_ca_for_ASA_proto,
    std::string &include_CH_aa_proto, std::string &include_CH_atom_proto,
    std::string &sphere_proto, std::string &cylinder_proto,
    std::string &prism_proto, std::string const &include_CH_filename,
    ANA_molecule &molecule_points, Point &cm,
    std::vector<unsigned int> &aa_list, std::vector<unsigned int> &CA_indices,
    std::vector<Point> &CAs_Points, std::vector<unsigned int> &include_CH_aa,
    Triang_Vector &CH_triangs, std::vector<unsigned int> &hetatm_atoms);
// Read coordinates in netcdf format.
void read_MD(const chemfiles::Frame &in_frame, bool const requested_CH,
    std::string const &sphere_proto, std::string const &cylinder_proto,
    std::string const &prism_proto,
    const std::vector<unsigned int> &hetatm_atoms,
    std::vector<unsigned int> &include_CH_atoms,
    std::string const &include_CH_filename, Triang_Vector &CH_triangs,
    std::string const &ASA_method, const std::vector<unsigned int> &CA_indices,
    std::vector<Point> &CAs_points, ANA_molecule &molecule_points);
// Get the center of mass from a chemfiles molecule.
Point getCM(const chemfiles::span<chemfiles::Vector3D> &in_xyz,
    unsigned int const natoms);
// Read PDB to draw included area for MD
void read_included_area(
    std::string const &filename, std::vector<Point> &area_points);
// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream);
}
#endif // _H
