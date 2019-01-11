#include <ANA/ConvexHull.hpp>

namespace ANA {

// Refine the provided list of amino acids. Throws a lot.
std::vector<unsigned int> string_to_list(
    std::string const &list_proto, unsigned int const top) {

    std::vector<unsigned int> list;

    std::stringstream stream_aa(list_proto);
    std::string temp_aa;
    while (!stream_aa.eof()) {
        unsigned int aa;
        stream_aa >> temp_aa;
        try {
            aa = std::stoi(temp_aa);
        } catch (std::out_of_range const oor) {
            // int is too large to be represented by int
            throw std::out_of_range(
                "Invalid atom / residue number in config file. Aborting.");
        } catch (...) {
            // some other exception.
            throw std::runtime_error(
                "Invalid atom / residue number. Aborting.");
        }
        list.push_back(aa);
    }
    // sort list of included amino acids
    std::sort(list.begin(), list.end());

    if (top != 0) {
        if (list[list.size() - 1] > top) {
            throw std::runtime_error(
                "Atom / residue list goes out of bounds. Check this input list "
                "and your input PDB atom / residue count. Aborting.");
        }
    }

    return list;
}

ConvexHull::ConvexHull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts) {

    switch (IA_opts._opt) {
    case IncludedAreaOptions::IAOption::residue:
        try {
            ConvexHull(protein, IA_opts._resn_proto, ResidueTag);
        } catch (...) {
            throw;
        }
        break;
    case IncludedAreaOptions::IAOption::atom:
        try {
            ConvexHull(protein, IA_opts._atom_proto, AtomTag);
        } catch (...) {
            throw;
        }
        break;
    case IncludedAreaOptions::IAOption::sphere:
        try {
            ConvexHull(protein, IA_opts._sphere_proto, SphereTag);
        } catch (...) {
            throw;
        }
        break;
    case IncludedAreaOptions::IAOption::cylinder:
        try {
            ConvexHull(protein, IA_opts._cylinder_proto, CylinderTag);
        } catch (...) {
            throw;
        }
        break;
    case IncludedAreaOptions::IAOption::prism:
        try {
            ConvexHull(protein, IA_opts._prism_proto, PrismTag);
        } catch (...) {
            throw;
        }
        break;
    case IncludedAreaOptions::IAOption::file:
        try {
            ConvexHull(protein, IA_opts._filename, FileTag);
        } catch (...) {
            throw;
        }
        break;
    }
}

ConvexHull::ConvexHull(
    Molecule const &protein, std::string const &resn_proto, ResidueTag) {

    auto const residues = string_to_list(resn_proto, protein._nres);

    std::vector<unsigned int> incl_area_points;
    incl_area_points.reserve(residues.size());

    for (auto const i : residues) {
        incl_area_points.push_back(protein._data[i - 1].first);
    }

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(
    Molecule const &protein, std::string const &atom_proto, AtomTag) {

    auto const atoms = string_to_list(atom_proto, protein._natoms);

    std::vector<unsigned int> incl_area_points;
    incl_area_points.reserve(atoms.size());

    for (auto const i : atoms) {
        // 0-index normalization.
        incl_area_points.push_back(protein._data[i - 1].first);
    }

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(std::string const &sphere_proto, SphereTag) {

    std::stringstream stream_sphere(sphere_proto);
    double const x = parse_double(stream_sphere);
    double const y = parse_double(stream_sphere);
    double const z = parse_double(stream_sphere);
    double const r = parse_double(stream_sphere);
    double constexpr cos_30 = 0.86602540378;
    double constexpr sin_30 = 0.5;

    Point const center(x, y, z);
    // Drawin a pseudo-sphere. The first 6 correspondon the XYZ axes, the next 8
    // to the X-Y plane, then 8 more for the X-Z plane and 8 for the Y-Z plane.
    std::array<Point, 30> const incl_area_points{center + Vector(r, 0, 0),
        center + Vector(0, r, 0), center + Vector(0, 0, r),
        center + Vector(-r, 0, 0), center + Vector(0, -r, 0),
        center + Vector(0, 0, -r), center + Vector(r * cos_30, r * sin_30, 0),
        center + Vector(r * sin_30, r * cos_30, 0),
        center + Vector(r * cos_30, -r * sin_30, 0),
        center + Vector(r * sin_30, -r * cos_30, 0),
        center + Vector(-r * cos_30, r * sin_30, 0),
        center + Vector(-r * sin_30, r * cos_30, 0),
        center + Vector(-r * cos_30, -r * sin_30, 0),
        center + Vector(-r * sin_30, -r * cos_30, 0),
        center + Vector(r * cos_30, 0, r * sin_30),
        center + Vector(r * sin_30, 0, r * cos_30),
        center + Vector(r * cos_30, 0, -r * sin_30),
        center + Vector(r * sin_30, 0, -r * cos_30),
        center + Vector(-r * cos_30, 0, r * sin_30),
        center + Vector(-r * sin_30, 0, r * cos_30),
        center + Vector(-r * cos_30, 0, -r * sin_30),
        center + Vector(-r * sin_30, 0, -r * cos_30),
        center + Vector(0, r * cos_30, r * sin_30),
        center + Vector(0, r * sin_30, r * cos_30),
        center + Vector(0, r * cos_30, -r * sin_30),
        center + Vector(0, r * sin_30, -r * cos_30),
        center + Vector(0, -r * cos_30, r * sin_30),
        center + Vector(0, -r * sin_30, r * cos_30),
        center + Vector(0, -r * cos_30, -r * sin_30),
        center + Vector(0, -r * sin_30, -r * cos_30)};

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(std::string const &cylinder_proto, CylinderTag) {

    std::stringstream stream_cylinder(cylinder_proto);

    double const x1 = parse_double(stream_cylinder);
    double const y1 = parse_double(stream_cylinder);
    double const z1 = parse_double(stream_cylinder);
    double const x2 = parse_double(stream_cylinder);
    double const y2 = parse_double(stream_cylinder);
    double const z2 = parse_double(stream_cylinder);
    double const r = parse_double(stream_cylinder);
    double constexpr cos_30 = 0.86602540378;
    double const sin_30 = 0.5;

    Point const center_1(x1, y1, z1), center_2(x2, y2, z2);

    Vector const vdiff(center_2 - center_1);
    Vector n1(-vdiff.y(), vdiff.x(), 0);
    Vector n2 = CGAL::cross_product(vdiff, n1);
    n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
    n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

    // The first 12 correspond to the first tap, the other half correspond to
    // the 2nd tap.
    std::array<Point, 24> incl_area_points{center_1 + r * n1, center_1 + r * n2,
        center_1 - r * n1, center_1 - r * n2,
        center_1 + r * cos_30 * n1 + r * sin_30 * n2,
        center_1 + r * sin_30 * n1 + r * cos_30 * n2,
        center_1 + r * cos_30 * n1 - r * sin_30 * n2,
        center_1 + r * sin_30 * n1 - r * cos_30 * n2,
        center_1 - r * cos_30 * n1 + r * sin_30 * n2,
        center_1 - r * sin_30 * n1 + r * cos_30 * n2,
        center_1 - r * cos_30 * n1 - r * sin_30 * n2,
        center_1 - r * sin_30 * n1 - r * cos_30 * n2, center_2 + r * n1,
        center_2 + r * n2, center_2 - r * n1, center_2 - r * n2,
        center_2 + r * cos_30 * n1 + r * sin_30 * n2,
        center_2 + r * sin_30 * n1 + r * cos_30 * n2,
        center_2 + r * cos_30 * n1 - r * sin_30 * n2,
        center_2 + r * sin_30 * n1 - r * cos_30 * n2,
        center_2 - r * cos_30 * n1 + r * sin_30 * n2,
        center_2 - r * sin_30 * n1 + r * cos_30 * n2,
        center_2 - r * cos_30 * n1 - r * sin_30 * n2,
        center_2 - r * sin_30 * n1 - r * cos_30 * n2};

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(std::string const &prism_proto, PrismTag) {

    std::stringstream stream_prism(prism_proto);
    double const x1 = parse_double(stream_prism);
    double const y1 = parse_double(stream_prism);
    double const z1 = parse_double(stream_prism);
    double const x2 = parse_double(stream_prism);
    double const y2 = parse_double(stream_prism);
    double const z2 = parse_double(stream_prism);
    double const width = parse_double(stream_prism) / 2;
    double const height = parse_double(stream_prism) / 2;

    Point const center_1(x1, y1, z1), center_2(x2, y2, z2);
    Vector const vdiff(center_2 - center_1);
    Vector n1(-vdiff.y(), vdiff.x(), 0);
    Vector n2 = CGAL::cross_product(vdiff, n1);
    n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
    n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

    // 8 vertices of a prism.
    std::array<Point, 8> incl_area_points{center_1 + width * n1 + height * n2,
        center_1 + width * n1 - height * n2,
        center_1 - width * n1 + height * n2,
        center_1 - width * n1 - height * n2,
        center_2 + width * n1 + height * n2,
        center_2 + width * n1 - height * n2,
        center_2 - width * n1 + height * n2,
        center_2 - width * n1 - height * n2};

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(
    Molecule const &protein, std::string const &filename, FileTag) {
    throw("Convex hull construction from input file not supported yet. "
          "Aborting.");
}

}

// Run the actual Convex Hull algorithm with CGAL.
void run_convex_hull(std::vector<unsigned int> const &points) {
    if (points.size() < 4) {
        throw std::runtime_error(
            "Not possible to triangulate less than 4 points. Aborting.");
    }

    Polyhedron CH;
    CGAL::convex_hull_3(points.begin(), points.end(), CH);
    P_Facet_const_iterator f_end = CH.facets_end();
    for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
         ++f_ite) {
        // Fix around the weirdest CGAL bug.
        P_Halfedge_around_facet_const_circulator he_ite = f_ite->facet_begin();
        auto const he_ite_0 = he_ite++;
        auto const he_ite_1 = he_ite++;
        auto const he_ite_2 = he_ite;

        _data.emplace_back(he_ite_0->vertex()->point(),
            he_ite_1->vertex()->point(), he_ite_2->vertex()->point());
    }
    return;
}
}
