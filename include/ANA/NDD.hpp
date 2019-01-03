#ifndef ANA_NDD_H
#define ANA_NDD_H

#include <ANA/Includes.hpp>
#include <ANA/Modes.hpp>
#include <ANA/NDDUtils.hpp>
#include <ANA/Read.hpp>
#include <ANA/Utils.hpp>
#include <ANA/Write.hpp>

namespace ANA {

// Main function for NDD version of ANA.
int NDD_ANA(std::string const &in_filename, std::string &include_CH_aa_proto,
    std::string &include_CH_atom_proto, std::string &sphere_proto,
    std::string &cylinder_proto, std::string &prism_proto,
    std::string const &include_CH_filename,
    std::string const &modes_ndd_filename,
    std::string const &pdbs_list_ndd_filename,
    std::string const &out_ndd_filename, bool const atom_only,
    double const minVR, double const maxSR);



} // namespace ANA

#endif // _H
