#ifndef ANA_OPTIONS_H
#define ANA_OPTIONS_H

#include <ANA/Includes.hpp>

namespace ANA {

struct IncludedAreaOptions {
public:
    // Would like to have only 1 string, but that won't play nicely with Boost
    // Program Options. Hopefully I'll fix it someday. TODO
    // TODO remove "include_CH_" from the member variables names. It's
    // redundant.
    std::string _resn_proto = "none";
    std::string _atom_proto = "none";
    std::string _sphere_proto = "none";
    std::string _cylinder_proto = "none";
    std::string _prism_proto = "none";
    std::string _filename = "none";
    enum IAOption { none, residue, atom, sphere, cylinder, prism, file };
    IAOption _opt = IAOption::none;
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