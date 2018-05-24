////////////////////////////////////////////////////////////////////////////////
//                                ANA                                         //
//                                                                            //
//                                                                            //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <ANA/ANAstatic.hpp>
#include <ANA/ANAmd.hpp>
#include <ANA/ANAndd.hpp>
#include <ANA/ANAPO.hpp>

int main(int argc, char* argv[]) {
    std::string in_filename, AA_indices_proto, exclude_ca_for_ASA_indices_proto,
        in_md_filename, out_type, pdbs_list_ndd_filename, out_CH_filename,
        out_ndd_filename, include_CH_filename,
        include_CH_aa_proto = "none", include_CH_atom_proto = "none",
        out_filename, out_vol, clusters_method = "facets",
        ASA_method = "dot_pdt", only_side_ASA, list_wall, list_wall_separator,
        tool_check_CH, tool_pdb_to_ch, sphere_proto, cylinder_proto,
        prism_proto, tool_pdb_norm, tool_aa_to_ca;

    bool triangulate_only_included_aas = false, atom_only = true;

    unsigned int nbr_of_vertices_to_include = 2, min_cells_cluster = 2,
                 estado = 1, md_start, md_step, md_end, precision, sphere_count;
    std::vector<unsigned int> include_CH_aa, hetatm_atoms;
    double minVR, maxSR, max_probe, max_probe_length, sphere_size;

    estado = ANA::get_parameters(argc, argv, in_filename, in_md_filename,
        include_CH_filename, include_CH_aa_proto, include_CH_atom_proto,
        AA_indices_proto, triangulate_only_included_aas, atom_only, precision,
        min_cells_cluster, nbr_of_vertices_to_include, md_start, md_step,
        md_end, minVR, maxSR, max_probe, max_probe_length, sphere_size,
        sphere_count, list_wall, list_wall_separator, clusters_method,
        only_side_ASA, ASA_method, exclude_ca_for_ASA_indices_proto,
        pdbs_list_ndd_filename, out_CH_filename, out_ndd_filename, out_filename,
        out_vol, out_type, tool_check_CH, tool_pdb_to_ch, sphere_proto,
        cylinder_proto, prism_proto, tool_pdb_norm, tool_aa_to_ca);
    if (estado == 1) {
        // Some error was found. Terminating execution.
        return 0;
    }

    //////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// tools
    ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    // Tool for CH included area visual check.
    if (tool_check_CH != "none") {
        Triang_Vector CH_triangs;
        std::vector<Point> CAs_Points;
        ANA_molecule molecule_points;
        std::vector<unsigned int> AA_indices, CA_indices, include_CH_atoms;
        Point cm;

        // Read input file
        const bool requested_CH = ANA::read_static(in_filename,
            triangulate_only_included_aas, atom_only, AA_indices_proto,
            exclude_ca_for_ASA_indices_proto, include_CH_aa_proto,
            include_CH_atom_proto, sphere_proto, cylinder_proto, prism_proto,
            include_CH_filename, molecule_points, cm, AA_indices, CA_indices,
            CAs_Points, include_CH_atoms, CH_triangs, hetatm_atoms);

        if (!requested_CH) {
            std::cerr << "No valid input for triangulation of included area. "
                         "Exiting now."
                      << '\n';
            return 0;
        }

        if (in_md_filename == "none") {

            ANA::draw_CH(CH_triangs, tool_check_CH);

        } else {

            chemfiles::Trajectory in_traj(in_md_filename);
            unsigned int frame_cnt = 1;
            // Set end.
            if (md_end == 0) {
                md_end = in_traj.nsteps();
            }
            // Set output file.
            tool_check_CH.append(".pdb");
            auto out_traj = chemfiles::Trajectory(tool_check_CH, 'w');

            while (!in_traj.done()) {
                // Set next step.
                const unsigned int current_step =
                    (frame_cnt - 1) * md_step + (md_start - 1);
                if (current_step >= md_end) {
                    // Done.
                    break;
                }
                // Read next frame.
                chemfiles::Frame in_frame = in_traj.read_step(current_step);
                // Update CH
                ANA::read_MD(in_frame, requested_CH, sphere_proto,
                    cylinder_proto, prism_proto, hetatm_atoms, include_CH_atoms,
                    include_CH_filename, CH_triangs, ASA_method, CA_indices,
                    CAs_Points, molecule_points);
                // New CH model.
                ANA::draw_CH(CH_triangs, out_traj);
                ++frame_cnt;
            }
        }
        return 0;
    }

    // Tool to get convex hull of input PDB and writing its vertices to a file.
    if (tool_pdb_to_ch != "none") {
        ANA::tool_PDB_to_CH(in_filename, tool_pdb_to_ch);
        return 0;
    }

    // Tool for normalizing PDB, by renumbering its atoms and residues.
    if (tool_pdb_norm != "none") {
        tool_pdb_norm.append(".pdb");
        ANA::tool_PDB_norm(in_filename, tool_pdb_norm);
        return 0;
    }

    // Tool for getting the indices of the Calpha atoms from the
    // included_area_residues config.
    if (tool_aa_to_ca != "none") {

        Triang_Vector CH_triangs;
        std::vector<Point> CAs_Points;
        ANA_molecule molecule_points;
        std::vector<unsigned int> AA_indices, CA_indices, include_CH_atoms;
        Point cm;

        // Read input file
        ANA::read_static(in_filename, triangulate_only_included_aas, atom_only,
            AA_indices_proto, exclude_ca_for_ASA_indices_proto,
            include_CH_aa_proto, include_CH_atom_proto, sphere_proto,
            cylinder_proto, prism_proto, include_CH_filename, molecule_points,
            cm, AA_indices, CA_indices, CAs_Points, include_CH_atoms,
            CH_triangs, hetatm_atoms);
        std::sort(include_CH_atoms.begin(), include_CH_atoms.end());

        std::cout << "\t\t/// Calpha indices ///" << '\n';
        for (const auto& each : include_CH_atoms) {
            std::cout << each + 1 << " ";
        }
        std::cout << '\n';

        return 0;
    }

    //////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// end tools
    ///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    // Get output volume file ready, if requested.
    ANA::open_vol_file(out_vol);

    if (in_md_filename == "none") {
        if (pdbs_list_ndd_filename == "none") {
            // Static PDB
            estado = ANA::static_ANA(in_filename, AA_indices_proto, ASA_method,
                only_side_ASA, exclude_ca_for_ASA_indices_proto, list_wall,
                list_wall_separator, clusters_method, include_CH_aa_proto,
                include_CH_atom_proto, sphere_proto, cylinder_proto,
                prism_proto, include_CH_filename, out_filename, out_vol,
                out_type, triangulate_only_included_aas, atom_only, minVR,
                maxSR, max_probe, max_probe_length, sphere_size, sphere_count,
                nbr_of_vertices_to_include, min_cells_cluster, precision);
            if (estado == 1) {
                return 0;
            }
        } else {
            // NDD
            estado = ANA::NDD_ANA(in_filename, AA_indices_proto, only_side_ASA,
                ASA_method, exclude_ca_for_ASA_indices_proto, list_wall,
                list_wall_separator, include_CH_aa_proto, include_CH_atom_proto,
                sphere_proto, cylinder_proto, prism_proto, include_CH_filename,
                pdbs_list_ndd_filename, out_ndd_filename, out_filename, out_vol,
                out_type, triangulate_only_included_aas, atom_only, minVR,
                maxSR, max_probe, max_probe_length, sphere_size, sphere_count,
                nbr_of_vertices_to_include, precision);
            if (estado == 1) {
                return 0;
            }
        }

    } else {
        // MD run.
        estado = ANA::MD_ANA(in_filename, in_md_filename, AA_indices_proto,
            ASA_method, only_side_ASA, exclude_ca_for_ASA_indices_proto,
            list_wall, list_wall_separator, include_CH_aa_proto,
            include_CH_atom_proto, sphere_proto, cylinder_proto, prism_proto,
            include_CH_filename, out_filename, out_vol, out_type,
            triangulate_only_included_aas, atom_only, minVR, maxSR, max_probe,
            max_probe_length, sphere_size, sphere_count,
            nbr_of_vertices_to_include, precision, md_start, md_step, md_end);
        if (estado == 1) {
            return 0;
        }
    }
    return 0;
}
