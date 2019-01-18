#include <ANA/Primitives.hpp>

namespace ANA {

Molecule::Molecule(std::string const &filename, bool const atom_only) {

    // Read PDB
    chemfiles::Trajectory input_pdb_traj(filename);
    auto const input_pdb_frame = input_pdb_traj.read();
    auto const in_xyz = input_pdb_frame.positions();
    auto const input_pdb_top = input_pdb_frame.topology();

    int const natoms = input_pdb_top.natoms();
    int const nres = input_pdb_top.residues().size();

    _data.reserve(natoms);
    _alphaCarbons.reserve(nres);

    for (auto const &residuo : input_pdb_top.residues()) {
        auto const res_name = residuo.name();
        for (auto const &i : residuo) {
            if (atom_only &&
                input_pdb_top[i].get("is_hetatm").value_or(false).as_bool()) {
                // Save the HETATM indices to discard them during MD and
                // NDD runs.
                _hetatms.push_back(i);
                continue;
            }

            auto const vdw_opt = input_pdb_top[i].vdw_radius();
            double vdw = 1.5;
            if (vdw_opt) {
                vdw = vdw_opt.value();
            } else {
                printf("Element for atom %i not available. Using Van Der Walls "
                       "radius of 1.5.\n",
                    static_cast<int>(i) + 1);
            }

            VertexInfo const vi1(i, vdw, residuo.id().value(), res_name);
            Point const p1(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
            _data.emplace_back(p1, vi1);

            if (input_pdb_top[i].name() == "CA") {
                // Will use to draw Convex Hull or whatever. It's always
                // useful.
                _alphaCarbons.push_back(i);
            }
        }
    }

    return;
}

Cavity::Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts) {

    _triangulation.insert(molecule._data.begin(), molecule._data.end());

    Finite_cells_iterator fc_ite_end = _triangulation.finite_cells_end();
    for (auto fc_ite = _triangulation.finite_cells_begin();
         fc_ite != fc_ite_end; ++fc_ite) {

        double const vol = volume(fc_ite);
        if (vol > cell_opts._min_CV) {
            _all_cells.push_back(fc_ite);
            _volume = _volume + vol;
        }
    }
}

void Cavity::add_border_tetra(Point const &p0, Point const &ip1,
    Point const &ip2, Point const &ip3, double const vdw0) {

    double const vol =
        volume(p0, ip1, ip2, ip3) - sphere_sector_vol(p0, ip1, ip2, ip3, vdw0);
    _outer_volume += vol;

    Polyhedron P;
    P.make_tetrahedron(p0, ip1, ip2, ip3);
    _border.push_back(std::move(P));
    _poly_vtx_cnt += 4;
    return;
}

void Cavity::add_border_penta(Point const &p0, Point const &p1,
    Point const &ip2, Point const &ip3, Point const &ip4, Point const &ip5,
    double const vdw0, double const vdw1) {

    double const vol1 = volume(p0, p1, ip2, ip3) -
        sphere_sector_vol(p0, p1, ip2, ip3, vdw0) -
        sphere_sector_vol(p1, p0, ip2, ip3, vdw1);

    double const vol2 =
        volume(p1, ip2, ip3, ip4) - sphere_sector_vol(p1, ip2, ip3, ip4, vdw1);

    double const vol3 =
        volume(p1, ip3, ip4, ip5) - sphere_sector_vol(p1, ip3, ip4, ip5, vdw1);

    _outer_volume += vol1 + vol2 + vol3;

    Polyhedron P;
    P.make_tetrahedron(p0, p1, ip2, ip3);
    P.make_tetrahedron(p1, ip2, ip3, ip4);
    P.make_tetrahedron(p1, ip3, ip4, ip5);
    _border.push_back(std::move(P));
    _poly_vtx_cnt += 12;
    return;
}

void Cavity::add_border_penta(Point const &p0, Point const &p1, Point const &p2,
    Point const &ip3, Point const &ip4, Point const &ip5, double const vdw0,
    double const vdw1, double const vdw2) {

    double const vol1 = volume(p0, p1, p2, ip3) -
        sphere_sector_vol(p0, p1, p2, ip3, vdw0) -
        sphere_sector_vol(p1, p0, p2, ip3, vdw1) -
        sphere_sector_vol(p2, p1, p0, ip3, vdw2);

    double const vol2 = volume(p1, p2, ip3, ip4) -
        sphere_sector_vol(p1, p2, ip3, ip4, vdw1) -
        sphere_sector_vol(p2, p1, ip3, ip4, vdw2);

    double const vol3 =
        volume(p2, ip3, ip4, ip5) - sphere_sector_vol(p2, ip3, ip4, ip5, vdw2);

    _outer_volume += vol1 + vol2 + vol3;

    Polyhedron P;
    P.make_tetrahedron(p0, p1, p2, ip3);
    P.make_tetrahedron(p1, p2, ip3, ip4);
    P.make_tetrahedron(p2, ip3, ip4, ip5);
    _border.push_back(std::move(P));
    _poly_vtx_cnt += 12;
    return;
}

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream) {
    std::string temp;
    double coord = 0;

    in_stream >> temp;
    try {
        coord = std::stof(temp);
    } catch (const std::invalid_argument &ia) {
        // some character present
        throw std::runtime_error("Can't parse input. There may be "
                                 "non-numerical characters. Aborting.");
    } catch (...) {
        // some other exception.
        throw std::runtime_error("Can't parse input. Aborting.");
    }

    return coord;
}

} // namespace ANA