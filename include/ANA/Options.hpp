#ifndef ANA_OPTIONS_H
#define ANA_OPTIONS_H

#include <ANA/Includes.hpp>

namespace ANA {

struct IncludedAreaOptions {
public:
    std::string _include_CH_resn_proto = "none";
    std::string _include_CH_atom_proto = "none";
    std::string _sphere_proto;
    std::string _cylinder_proto;
    std::string _prism_proto;
    std::string _include_CH_filename;
};

struct NDDOptions {
public:
    std::string _modes_ndd_filename;
    std::string _pdbs_list_ndd_filename;
    std::string _out_ndd_filename;
};

struct CellFilteringOptions {
public:
    double _minVR;
    double _maxSR;
};

}

#endif // _H