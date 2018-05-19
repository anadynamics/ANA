#ifndef ANAPO
#define ANAPO
namespace ANA {
namespace PO = boost::program_options;
using namespace std;
const std::string help_header = "\t\t\t----- ANA -----";

int get_parameters(int ac, char* av[], std::string& input_struct_filename,
    std::string& input_md_filename, std::string& include_CH_filename,
    std::string& include_CH_aa_proto, std::string& include_CH_atom_proto,
    std::string& AA_indices_proto, bool& triangulate_only_included_aas,
    bool& atom_only, unsigned int& precision, unsigned int& clusters_min_size,
    unsigned int& nbr_of_vertices_to_include, unsigned int& md_start,
    unsigned int& md_step, unsigned int& md_end, double& minVR, double& maxSR,
    double& max_probe, double& max_probe_length, double& sphere_size,
    unsigned int& sphere_count, std::string& list_wall,
    std::string& list_wall_separator, std::string& clusters_method,
    std::string& only_side_ASA, std::string& ASA_method,
    std::string& exclude_ca_for_ASA, std::string& pdbs_list_ndd_filename,
    std::string& out_CH_filename, std::string& out_ndd_filename,
    std::string& out_filename, std::string& out_vol, std::string& output_type,
    std::string& tool_check_CH, std::string& tool_pdb_to_ch,
    std::string& sphere_proto, std::string& cylinder_proto,
    std::string& prism_proto, std::string& tool_pdb_norm,
    std::string& tool_aa_to_ca);
}
#endif
