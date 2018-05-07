#ifndef ANAMD
#define ANAMD
namespace ANA {
// Main function for dynamic version of ANA.
int MD_ANA(const std::string &in_filename, const std::string &in_md_filename,
           std::string &AA_indices_proto, const std::string &only_side_ASA,
           const std::string &ASA_method,
           std::string &exclude_ca_for_ASA_indices_proto,
           const std::string &list_wall, const std::string &list_wall_separator,
           std::string &include_CH_aa_proto, std::string &include_CH_atom_proto,
           std::string &sphere_proto, std::string &cylinder_proto,
           std::string &prism_proto, const std::string &include_CH_filename,
           std::string &out_filename, const std::string &out_type,
           const bool triangulate_only_included_aas,
           const bool atom_only, const double minVR,
           const double maxSR, const double max_probe,
           const double max_probe_length, const double sphere_size,
           const unsigned int &sphere_count,
           const unsigned int nbr_of_vertices_to_include,
           const unsigned int precision, const unsigned int &md_start,
           const unsigned int &md_step, unsigned int &md_end);
}
#endif
