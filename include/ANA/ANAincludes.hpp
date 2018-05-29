////////////////////////////////////////////////////////////////////////////////
//                                 ANA                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ANAINCLUDES
#define ANAINCLUDES
#include  <iostream>
#include  <fstream>
#include  <assert.h>
#include  <vector>
#include  <numeric>
#include  <cmath>
#include  <algorithm>
#include  <string>
#include  <sstream>
#include  <iterator>
#include  <typeinfo>
#include  <boost/program_options.hpp>
#include  "chemfiles.hpp"

#include  <CGAL/Origin.h>
#include  <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include  <CGAL/Triangulation_3.h>
#include  <CGAL/Delaunay_triangulation_3.h>
#include  <CGAL/Triangulation_vertex_base_with_info_3.h>
#include  <CGAL/Tetrahedron_3.h>
#include  <CGAL/Polyhedron_3.h>
#include  <CGAL/Triangle_3.h>
#include  <CGAL/convex_hull_3.h>
#include  <CGAL/box_intersection_d.h>

// clang-format off
class Vtx_info  {
public:
// Return the atom index
	const unsigned int& GetIndex() const    { return index_; }
// Return the atom radii
	const double& GetRadii() const  { return radii_; }
// Return the atom's amino acid name
	const std::string& GetAa() const    { return residue_; }
// Return the atom's residue number
	const unsigned int& GetResn() const { return resn_; }

// Assign the atom index
	void AssignIndex(unsigned int input_index) { index_ = input_index; }
// Assign the atom radii
	void AssignRadii(double input_radii) { radii_ = input_radii; }
// Assign the atom's amino acid name
	void AssignAa(std::string input_residue) { residue_ = input_residue; }
// Assign the atom's residue number
	void AssignResn(unsigned int input_resn) { resn_ = input_resn; }

private:
unsigned int index_;
double radii_;
// atom's residue name in 3 letter format
std::string residue_;
// atom's residue number
unsigned int resn_;
};

// Colours for pymol CGO objects
const unsigned int col_nbr = 16;
const std::array<float, col_nbr> red = {1.0, 1.0, 0.0,   0.0, 1.0, 1.0, 0.0, 0.4,  0.4,
 0.0,  0.75, 0.75,  0.0, 0.875, 0.875,   0.0};
const std::array<float, col_nbr> green = {1.0, 0.0, 1.0,   0.0, 1.0, 0.0, 1.0, 0.4,  0.0,
 0.4,  0.75,  0.0, 0.75, 0.875,   0.0, 0.875};
const std::array<float, col_nbr> blue = {1.0, 0.0, 0.0,   1.0, 0.0, 1.0, 1.0,   0,  0.4,
 0.4,   0.0, 0.75, 0.75,   0.0, 0.875, 0.875};

// Basic definitions
 using  EPIC                                    = CGAL::Exact_predicates_inexact_constructions_kernel;
 using  Vb                                      = CGAL::Triangulation_vertex_base_with_info_3<Vtx_info, EPIC>;
 using  Tds                                     = CGAL::Triangulation_data_structure_3<Vb>;
 using  Delaunay                                = CGAL::Delaunay_triangulation_3<EPIC, Tds>;

// Definitions for Delaunay triangulation
 using Vector                                   = CGAL::Vector_3<EPIC>;
 using Segment                                  = CGAL::Segment_3<EPIC>;
 using Point                                    = Delaunay::Point;
 using Vertex_handle                            = Delaunay::Vertex_handle;
 using Cell_handle                              = Delaunay::Cell_handle;
 using Vertex_iterator                          = Delaunay::Vertex_iterator;
 using Finite_vertices_iterator                 = Delaunay::Finite_vertices_iterator;
 using All_vertices_iterator                    = Delaunay::All_vertices_iterator;
 using Edge_iterator                            = Delaunay::Edge_iterator;
 using Facet_iterator                           = Delaunay::Facet_iterator;
 using Finite_facets_iterator                   = Delaunay::Finite_facets_iterator;
 using Cell_iterator                            = Delaunay::Cell_iterator;
 using Finite_cells_iterator                    = Delaunay::Finite_cells_iterator;
 using All_cells_iterator                       = Delaunay::All_cells_iterator;
 using Cell_circulator                          = Delaunay::Cell_circulator;

// Definitions for convex hull
using Polyhedron                               = CGAL::Polyhedron_3<EPIC>;
using P_Facet_iterator                         = Polyhedron::Facet_iterator;
using P_Facet_const_iterator                   = Polyhedron::Facet_const_iterator;
using P_Edge_iterator                          = Polyhedron::Edge_iterator;
using P_Edge_const_iterator                    = Polyhedron::Edge_const_iterator;
using P_Halfedge_around_facet_circulator       = Polyhedron::Halfedge_around_facet_circulator;
using P_Halfedge_around_facet_const_circulator = Polyhedron::Halfedge_around_facet_const_circulator;
using P_Vertex_iterator                        = Polyhedron::Vertex_iterator;
using P_Vertex_const_iterator                  = Polyhedron::Vertex_const_iterator;

// ANA definitions
using MD_Element                               = std::array<Point, 4>; // Cell
using MD_Vector                                = std::vector<MD_Element>; // Pocket
using MD_Matrix                                = std::vector<MD_Vector>; // All voids
using NDD_Element                              = std::array<std::pair<Point, double>, 4>; // Cell
using NDD_Vector                               = std::vector<NDD_Element>; // Pocket
using NDD_Matrix                               = std::vector<NDD_Vector>; // All voids
using NDD_IElement                             = std::array<unsigned int, 4>; // Cell indices
using NDD_IVector                              = std::vector<NDD_IElement>; // Pocket indices
using NDD_IMatrix                              = std::vector<NDD_IVector>; // All voids indices
using NA_Vector                                = std::vector<Finite_cells_iterator>; // Pocket
using NA_Matrix                                = std::vector<NA_Vector>; // All voids
using Poly_Vector                              = std::vector<Polyhedron>; // Pocket border cells
using Poly_Matrix                              = std::vector<Poly_Vector>; // All null areas border cells
using ANA_molecule                             = std::vector<std::pair<Point, Vtx_info>>;

// Miscellaneous definitions
using Object                                   = CGAL::Object;
using Triangle                                 = CGAL::Triangle_3<EPIC>;
using Triang_Vector                            = std::vector<Triangle>;
using Segment                                  = EPIC::Segment_3;
using Tetrahedron                              = CGAL::Tetrahedron_3<EPIC>;
using Tetra_Vector                             = std::vector<Tetrahedron>;
using Box 																		 = CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Finite_cells_iterator>;
#endif
