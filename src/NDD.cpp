#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(std::string const &in_filename, IncludedAreaOptions &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    Molecule const protein = ANA::Molecule(in_filename, atom_only);

    // IA_opts._resn_proto = "1 12 21 51 68 97 111 119 71 130 133 128 16 24 84";
    // IA_opts._opt = IncludedAreaOptions::residue;

    ANA::ConvexHull const CH = create_convex_hull(protein, IA_opts);

    ANA::Cavity hueco(protein, cell_opts);

    ANA::carve_CH_into_cavity(hueco, CH);

    // ANA::NDD::ndd(cavity_joint_cells, NDD_opts);

    // ANA::write_output_volume(cavity_void_cells, poly_vol);

    ANA::draw(hueco, "sal.pdb");
    ANA::draw(CH, "hull.pdb");

    return 0;
}

} // namespace ANA
