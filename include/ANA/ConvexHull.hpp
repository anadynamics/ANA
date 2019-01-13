#ifndef ANA_CONVEX_HULL_H
#define ANA_CONVEX_HULL_H

#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

// Refine the provided list of amino acids. Throws a lot.
std::vector<unsigned int> adapt_AA_list(
    std::string &aa_list_proto, unsigned int const top = 0);

struct ResidueTag {};
struct AtomTag {};
struct SphereTag {};
struct CylinderTag {};
struct PrismTag {};
struct FileTag {};

struct ConvexHull {
public:
    ConvexHull() = default;

    ConvexHull(Molecule const &protein, IncludedAreaOptions const &IA_opts);

    // Atom and residue constructors are identical for now. I'll probably do
    // something different in the future.
    ConvexHull(
        Molecule const &protein, std::string const &resn_proto, ResidueTag);

    ConvexHull(Molecule const &protein, std::string const &atom_proto, AtomTag);

    ConvexHull(std::string const &sphere_proto, SphereTag);

    ConvexHull(std::string const &cylinder_proto, CylinderTag);

    ConvexHull(std::string const &prism_proto, PrismTag);

    ConvexHull(std::string const &filename, FileTag);

    // Run the actual Convex Hull algorithm with CGAL.
    // Should only work with containers. Will fix w/ c++20. TODO
    template <class T>
    void run_convex_hull(T const &points) {
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
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            _data.emplace_back(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point());
        }
        return;
    }

    std::vector<Triangle> _data;
};

}

#endif // _H