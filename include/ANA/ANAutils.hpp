#ifndef ANAUTILS
#define ANAUTILS
#include <ANA/ANAincludes.hpp>
#include <ANA/ANAread.hpp>
#include <ANA/ANAwrite.hpp>
// prototype functions
namespace ANA {
// Just calculate volume of the cell.
inline double cell_volume(Finite_cells_iterator const &cell_iterator) {
    return CGAL::to_double(CGAL::volume(cell_iterator->vertex(0)->point(),
        cell_iterator->vertex(1)->point(), cell_iterator->vertex(2)->point(),
        cell_iterator->vertex(3)->point()));
}
// Substract the volume filled with the 4 atoms from the total volume of
// the corresponding cell
double refine_cell_volume(
    double const entire_cell_vol, Finite_cells_iterator const &cell_iterator);
// Get the volume occupied by the sector of the sphere inscribed in the
// incident cell
double sphere_sector_vol(Point const &p_0, Point const &p_1, Point const &p_2,
    Point const &p_3, double const radius);
// Cluster neighbouring cells
void cluster_cells_cgal(NA_Vector const &input_cells, NA_Matrix &output_cells,
    unsigned int const min_cells_cluster);
// Given a cell, get all neighbouring cells that haven't been discovered yet
void get_neighbors(NA_Vector const &input_cells,
    Finite_cells_iterator const &query_cell, std::vector<unsigned int> &except,
    NA_Vector &output_cells);
// Cluster neighbouring cells. Iso-oriented boxes method
void cluster_cells_boxes(NA_Vector const &input_cells, NA_Matrix &output_cells);
// Keep cells that correspond to the included amino acids
void keep_included_aa_cells(NA_Vector const &input_cells,
    const std::vector<unsigned int> &aa_list,
    unsigned int const nbr_of_vertices_to_include, NA_Vector &output_cells);
// Discard exposed cells
void discard_ASA_dot_pdt_cm(Point const &cm,
    std::vector<Point> const &Calpha_xyz, double const min_dot,
    double const max_length, std::string const only_side_ASA,
    NA_Vector const &input_cells, NA_Vector &output_cells);
// Discard exposed cells. Calpha convex hull method
void discard_ASA_CACH(std::vector<Point> const &Calpha_xyz,
    std::string const only_side_ASA, NA_Vector const &input_cells,
    NA_Vector &output_cells);
// Discard exposed cells. Alpha shape method
void discard_ASA_dot_pdt_axes(std::vector<Point> const &Calpha_xyz,
    double const min_dot, double const max_length,
    std::string const only_side_ASA, NA_Vector const &input_cells,
    NA_Vector &output_cells);

// Triangulate only the specified amino acids or the whole molecule.
inline Delaunay triangulate(ANA_molecule const &molecule_points) {
    Delaunay T;
    T.insert(molecule_points.begin(), molecule_points.end());
    return T;
}
// Calc volume of the input cells.
inline double get_void_volume(NA_Vector const &input_cells) {
    double volume = 0;
    double current_cell_vol;
    for (Finite_cells_iterator const &fc_ite : input_cells) {
        current_cell_vol = cell_volume(fc_ite);
        current_cell_vol = refine_cell_volume(current_cell_vol, fc_ite);
        volume = volume + current_cell_vol;
    }
    return volume;
}
// Disregard lone cells.
inline void disregard_lone_cells(NA_Vector const &input_cells,
    NA_Vector &output_cells, unsigned int const lone_threshold) {

    unsigned int i, j, cell_cnt = input_cells.size();
    ;
    Finite_cells_iterator cell_ite1, cell_ite2;
    for (i = 0; i != cell_cnt; ++i) {
        cell_ite1 = input_cells[i];
        for (j = 1; j != cell_cnt; ++j) {
            cell_ite2 = input_cells[j];
            if (cell_ite1->has_neighbor(cell_ite2)) {
                output_cells.push_back(cell_ite1);
                break;
            }
        }
    }
}
// Calculate area of the cell's facets
inline void cell_facets_areas(Finite_cells_iterator const &cell_iterator,
    double &facet_0_area, double &facet_1_area, double &facet_2_area,
    double &facet_3_area) {
    facet_0_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(1)->point(), cell_iterator->vertex(2)->point(),
        cell_iterator->vertex(3)->point()));
    facet_1_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(2)->point(), cell_iterator->vertex(3)->point(),
        cell_iterator->vertex(0)->point()));
    facet_2_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(0)->point(), cell_iterator->vertex(1)->point(),
        cell_iterator->vertex(3)->point()));
    facet_3_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(0)->point(), cell_iterator->vertex(1)->point(),
        cell_iterator->vertex(2)->point()));
    return;
}
// Determine if any of the facets of the given cells has area larger
// than criteria
inline int refine_cell_areas(
    Finite_cells_iterator const cell_iterator, double const criteria) {

    double f0_area;
    double f1_area;
    double f2_area;
    double f3_area;

    cell_facets_areas(cell_iterator, f0_area, f1_area, f2_area, f3_area);
    if (f0_area > criteria || f1_area > criteria || f2_area > criteria ||
        f3_area > criteria) {
        return 1;
    } else
        return 0;
}
// Fill the 2 input vectors with iterators for the outer and inner cells
// respectively
void partition_triangulation(
    Delaunay const &T, NA_Vector &outer_cells, NA_Vector &inner_cells);
// Calc volume and get the proper cells
double get_all_voids(Delaunay const &T, NA_Vector &big_cells,
    double const min_vol_radius, double const max_area_radius);
// Discard cells without a vertex inside the specified convex hull. Lo
// precision.
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells);
// Discard cells without a vertex inside the specified convex hull. Hi
// precision.
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells, NA_Vector &out_intersecting_cells,
    std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<unsigned int> &intersecting_total);
// Discard parts of cells outside the specified triangulation using
// intersecitons.
double discard_CH_1(NA_Vector const &in_intersecting_cells,
    Triang_Vector const &CH_triangs,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    const std::vector<unsigned int> &intersecting_total,
    Poly_Vector &border_poly,
    std::vector<std::array<double, 3>> &in_vtces_radii,
    unsigned int &atom_cnt_poly);
// Discard cells without a vertex inside the specified convex hull
void discard_CH(
    NA_Vector const &in_cells, Polyhedron &CH, NA_Vector &out_cells);
// Extract vertices coordinates from the cells and store them in an "md_vector".
void na_vector_into_ndd_vector(
    NA_Vector const &in_cells, NDD_Vector &out_cells);
// Tool for reading PDB to draw included area.
void tool_PDB_to_CH(
    std::string const &in_filename, std::string const &out_filename);
// Tool for normalizing PDB, by renumbering its atoms and residues.
void tool_PDB_norm(
    std::string const &in_filename, std::string const &tool_pdb_norm);
// Helper function for inserting new elements in ordered vectors.
template <class Vector, class T_to_insert>
void insert_into_ord_vtor(Vector &v, const T_to_insert &to_insert);
// Helper function for getting the indices that sort a vector.
template <typename T>
std::vector<unsigned int> sort_indices(const std::vector<T> &v);
// Helper function to use binary search to find the lowest bound of a query in
// an unsorted vector and vector of indices that sorts it in ascending order. It
// returns the index of the element that satisfies the lower bound condition.
template <typename v_type>
bool lb_with_indices(const std::vector<v_type> &v1,
    const std::vector<unsigned int> &indices, const v_type q1,
    unsigned int &first);
// Helper function to use binary search to find the lowest bound of a query in
// a sorted vector in ascending order. It returns true if a match is found, and
// stores the index of the element that satisfies the lower bound condition in
// the variable "first".
template <typename v_type>
bool lb(const std::vector<v_type> &v1, const v_type q1, unsigned int &first);
// Helper function for taking the "i" number in the vector that doesn't match
// the query
template <class query_type>
query_type get_i_not_equal(const std::vector<query_type> &in_vec,
    const query_type &query, unsigned int const i);
// Helper function for taking the "i" number in the 'in_vec'' that doesn't
// match the query vector
template <class query_type>
query_type get_i_not_equal(const std::vector<query_type> &in_vec,
    const std::vector<query_type> &query_vec, unsigned int const i);
} // namespace ANA

namespace ANA {
namespace NDD {
    // Analytical NDD
    void ndd(NA_Vector const &cavity_void_cells,
        std::string const &modes_ndd_filename, std::string const &out_file);
    // Perform Non Delaunay Dynamics.
    void ndd_nondelaunay_dynamics_old(NA_Vector const &cavity_void_cells,
        std::string const &pdb_list, bool const precision,
        const std::vector<unsigned int> include_CH_atoms,
        std::string const &out_file);
    // Get the indices of the atoms involved in the given cells
    NDD_IVector get_vertices(NA_Vector const &cavity_void_cells);
    // Calc volume of the input cells. Reedited for array container.
    double ndd_get_void_volume(NDD_Vector const &cavity_void_cells);
    // Substract the volume filled with the 4 atoms from the total volume of
    // the corresponding cell. Reedited for array container.
    inline double ndd_refine_cell_volume(
        double const entire_cell_vol, const NDD_Element &cell) {
        Point const p_0 = cell[0].first;
        Point const p_1 = cell[1].first;
        Point const p_2 = cell[2].first;
        Point const p_3 = cell[3].first;

        double const vtx_0_sphere_sector =
            sphere_sector_vol(p_0, p_1, p_2, p_3, cell[0].second);
        double const vtx_1_sphere_sector =
            sphere_sector_vol(p_3, p_0, p_1, p_2, cell[3].second);
        double const vtx_2_sphere_sector =
            sphere_sector_vol(p_2, p_3, p_0, p_1, cell[2].second);
        double const vtx_3_sphere_sector =
            sphere_sector_vol(p_1, p_2, p_3, p_0, cell[1].second);
        return (entire_cell_vol - vtx_0_sphere_sector - vtx_1_sphere_sector -
            vtx_2_sphere_sector - vtx_3_sphere_sector);
    }
    // Discard cells without a vertex inside the specified convex hull. Hi
    // precision, NDD version
    void ndd_discard_CH_0(NDD_Vector const &in_cells,
        Triang_Vector const &CH_triangs, NDD_Vector &out_cells,
        NDD_Vector &out_intersecting_cells,
        std::vector<std::array<bool, 4>> &intersecting_bool,
        std::vector<unsigned int> &intersecting_total);
    // Discard parts of cells outside the specified triangulation using
    // intersecitons
    double ndd_discard_CH_1(NDD_Vector const &in_intersecting_coords,
        Triang_Vector const &CH_triangs,
        const std::vector<std::array<bool, 4>> &intersecting_bool,
        const std::vector<unsigned int> &intersecting_total,
        Poly_Vector &border_poly);
} // namespace NDD
} // namespace ANA
#endif
