#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(std::string const &in_filename, IncludedAreaOptions &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    Molecule const protein = ANA::Molecule(in_filename, atom_only);

    ANA::ConvexHull const CH(protein, IA_opts);

    ANA::Cavity hueco(protein, cell_opts);

    ANA::carve_CH_into_cavity(hueco, CH);

    // ANA::NDD::ndd(cavity_joint_cells, NDD_opts);

    // ANA::write_output_volume(cavity_void_cells, poly_vol);

    std::string sal{"sal.pdb"};
    ANA::draw(hueco, sal);

    return 0;
}

} // namespace ANA
