#ifndef ANASTATIC
#define ANASTATIC
#include <ANA/ANAincludes.hpp>
#include <ANA/ANAread.hpp>
#include <ANA/ANAutils.hpp>
#include <ANA/ANAwrite.hpp>
// prototype functions
namespace ANA {
// Main function for static version of ANA.
int static_ANA(std::string const &in_filename, std::string &AA_indices_proto,
    std::string const &ASA_method, std::string const &only_side_ASA,
    std::string &exclude_ca_for_ASA_indices_proto, std::string const &list_wall,
    std::string const &list_wall_separator, std::string const &clusters_method,
    std::string &include_CH_aa_proto, std::string &include_CH_atom_proto,
    std::string &sphere_proto, std::string &cylinder_proto,
    std::string &prism_proto, std::string const &include_CH_filename,
    std::string &out_filename, std::string const &out_vol,
    std::string const &out_type, bool const triangulate_only_included_aas,
    bool const atom_only, double const minVR, double const maxSR,
    double const max_probe, double const max_probe_length,
    double const sphere_size, unsigned int const sphere_count,
    unsigned int const nbr_of_vertices_to_include,
    unsigned int const clusters_min_size, unsigned int const precision);
} // namespace ANA
#endif
