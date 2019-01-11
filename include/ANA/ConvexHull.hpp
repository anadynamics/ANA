#ifndef ANA_CONVEX_HULL_H
#define ANA_CONVEX_HULL_H

#include <ANA/Includes.hpp>

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
    void run_convex_hull(std::vector<unsigned int> const &points);

    std::vector<Triangle> _data;
};

}

#endif // _H