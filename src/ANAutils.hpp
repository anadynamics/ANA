#ifndef ANA_UTILS_H
#define ANA_UTILS_H
// prototype functions
namespace ANA {
namespace NDD {
// Perform Non Delaunay Dynamics. Lo precision
void ndd_nondelaunay_dynamics(const NA_Vector &cavity_void_cells,
                              const std::string &pdb_list, const bool precision,
                              const std::vector<unsigned int> include_CH_atoms,
                              const std::string &out_file);
// Get the indices of the atoms involved in the given cells
void ndd_get_involved_vertices(const NA_Vector &cavity_void_cells,
                               NDD_IVector &cells_indices);
// Calc volume of the input cells. Reedited for array container.
double ndd_get_void_volume(const NDD_Vector &cavity_void_cells);
// Substract the volume filled with the 4 atoms from the total volume of
// the corresponding cell. Reedited for array container.
double ndd_refine_cell_volume(double entire_cell_vol, const NDD_Element &cell);
// Discard cells without a vertex inside the specified convex hull. Hi
// precision, NDD version
void ndd_discard_CH_0(const NDD_Vector &in_cells,
                      const Triang_Vector &CH_triangs, NDD_Vector &out_cells,
                      NDD_Vector &out_intersecting_cells,
                      std::vector<std::array<bool, 4>> &intersecting_bool,
                      std::vector<unsigned int> &intersecting_total);
// Discard parts of cells outside the specified triangulation using
// intersecitons
double
ndd_discard_CH_1(const NDD_Vector &in_intersecting_coords,
                 const Triang_Vector &CH_triangs,
                 const std::vector<std::array<bool, 4>> &intersecting_bool,
                 const std::vector<unsigned int> &intersecting_total,
                 Poly_Vector &border_poly);
}

// Cluster neighbouring cells
void cluster_cells_cgal(const NA_Vector &input_cells, NA_Matrix &output_cells,
                        const unsigned int min_cells_cluster);
// Given a cell, get all neighbouring cells that haven't been discovered yet
void get_neighbors(const NA_Vector &input_cells,
                   const Finite_cells_iterator &query_cell,
                   std::vector<unsigned int> &except, NA_Vector &output_cells);
// Cluster neighbouring cells. Iso-oriented boxes method
void cluster_cells_boxes(const NA_Vector &input_cells, NA_Matrix &output_cells);
// Keep cells that correspond to the included amino acids
void keep_included_aa_cells(const NA_Vector &input_cells,
                            const std::vector<unsigned int> &aa_list,
                            const unsigned int nbr_of_vertices_to_include,
                            NA_Vector &output_cells);
// Discard exposed cells
void discard_ASA_dot_pdt_cm(const Point &cm,
                            const std::vector<Point> &Calpha_xyz,
                            const double min_dot, const double max_length,
                            const std::string only_side_ASA,
                            const NA_Vector &input_cells,
                            NA_Vector &output_cells);
// Discard exposed cells. Calpha convex hull method
void discard_ASA_CACH(const std::vector<Point> &Calpha_xyz,
                      const std::string only_side_ASA,
                      const NA_Vector &input_cells, NA_Vector &output_cells);
// Discard exposed cells. Alpha shape method
void discard_ASA_dot_pdt_axes(const std::vector<Point> &Calpha_xyz,
                              const double min_dot, const double max_length,
                              const std::string only_side_ASA,
                              const NA_Vector &input_cells,
                              NA_Vector &output_cells);
// Triangulate only the specified amino acids or the whole molecule
Delaunay triangulate(ANA_molecule &molecule_points);
// Calc volume of the input cells
double get_void_volume(const NA_Vector &input_cells);
// Disregard lone cells
void disregard_lone_cells(const NA_Vector &input_cells, NA_Vector &output_cells,
                          unsigned int lone_threshold);
// Determine if any of the facets of the given cells has area larger
// than criteria
int refine_cell_areas(const Finite_cells_iterator cell_iterator,
                      const double criteria);
// Calculate area of the cell's facets
void cell_facets_areas(const Finite_cells_iterator &cell_iterator,
                       double &facet_0_area, double &facet_1_area,
                       double &facet_2_area, double &facet_3_area);
// Fill the 2 input vectors with iterators for the outer and inner cells
// respectively
void partition_triangulation(const Delaunay &T, NA_Vector &outer_cells,
                             NA_Vector &inner_cells);
// Calc volume and get the proper cells
double get_all_voids(const Delaunay &T, NA_Vector &big_cells,
                     const double min_vol_radius, const double max_area_radius);
// Substract the volume filled with the 4 atoms from the total volume of
// the corresponding cell
double refine_cell_volume(const double entire_cell_vol,
                          const Finite_cells_iterator &cell_iterator);
// Just calculate volume of the cell
double cell_volume(const Finite_cells_iterator &cell_iterator);
// Get the volume ocuppied by the sector of the sphere inscribed in the
// incident cell
double sphere_sector_vol(const Point &p_0, const Point &p_1, const Point &p_2,
                         const Point &p_3, const double radius);
// Discard cells without a vertex inside the specified convex hull. Lo
// precision.
void discard_CH_0(const NA_Vector &in_cells, const Triang_Vector &CH_triangs,
                  NA_Vector &out_cells);
// Discard cells without a vertex inside the specified convex hull. Hi
// precision.
void discard_CH_0(const NA_Vector &in_cells, const Triang_Vector &CH_triangs,
                  NA_Vector &out_cells, NA_Vector &out_intersecting_cells,
                  std::vector<std::array<bool, 4>> &intersecting_bool,
                  std::vector<unsigned int> &intersecting_total);
// Discard parts of cells outside the specified triangulation using
// intersecitons.
double discard_CH_1(const NA_Vector &in_intersecting_cells,
                    const Triang_Vector &CH_triangs,
                    const std::vector<std::array<bool, 4>> &intersecting_bool,
                    const std::vector<unsigned int> &intersecting_total,
                    Poly_Vector &border_poly,
                    std::vector<std::array<double, 3>> &in_vtces_radii,
                    unsigned int &atom_cnt_poly);
// Discard cells without a vertex inside the specified convex hull
void discard_CH(const NA_Vector &in_cells, Polyhedron &CH,
                NA_Vector &out_cells);
// Extract vertices coordinates from the cells and store them in an "md_vector".
void na_vector_into_ndd_vector(const NA_Vector &in_cells,
                               NDD_Vector &out_cells);
// Tool for reading PDB to draw included area.
void tool_PDB_to_CH(const std::string &in_filename,
                    const std::string &out_filename);
// Tool for normalizing PDB, by renumbering its atoms and residues.
void tool_PDB_norm(const std::string &in_filename,
                   const std::string &tool_pdb_norm);
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
                           const query_type &query, const unsigned int i);
// Helper function for taking the "i" number in the 'in_vec'' that doesn't
// match the query vector
template <class query_type>
query_type get_i_not_equal(const std::vector<query_type> &in_vec,
                           const std::vector<query_type> &query_vec,
                           const unsigned int i);
}
#endif