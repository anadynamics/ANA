#ifndef ANAINCLUDES
#define ANAINCLUDES

#include "chemfiles.hpp"
#include <algorithm>
#include <array>
#include <assert.h>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <typeinfo>
#include <utility>
#include <vector>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Origin.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/convex_hull_3.h>

struct VertexInfo {
public:
    VertexInfo() = default;

    VertexInfo(unsigned int const index, double const radius,
        unsigned int const resn, std::string_view const resi) :
        _index(index),
        _radius(_radius), _resn(resn), _resi(resi) {}

    unsigned int _index;
    double _radius;
    // atom's residue number
    unsigned int _resn;
    // atom's residue name in 3 letter format
    std::string _resi;
};

// clang-format off
// Basic definitions
using EPIC = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, EPIC>;
using Tds = CGAL::Triangulation_data_structure_3<Vb>;
using Delaunay = CGAL::Delaunay_triangulation_3<EPIC, Tds>;

// Definitions for Delaunay triangulation
using Vector = CGAL::Vector_3<EPIC>;
using Segment = CGAL::Segment_3<EPIC>;
using Point = Delaunay::Point;
using Vertex_handle = Delaunay::Vertex_handle;
using Cell_handle = Delaunay::Cell_handle;
using Vertex_iterator = Delaunay::Vertex_iterator;
using Finite_vertices_iterator = Delaunay::Finite_vertices_iterator;
using All_vertices_iterator = Delaunay::All_vertices_iterator;
using Edge_iterator = Delaunay::Edge_iterator;
using Facet_iterator = Delaunay::Facet_iterator;
using Finite_facets_iterator = Delaunay::Finite_facets_iterator;
using Cell_iterator = Delaunay::Cell_iterator;
using Finite_cells_iterator = Delaunay::Finite_cells_iterator;
using All_cells_iterator = Delaunay::All_cells_iterator;
using Cell_circulator = Delaunay::Cell_circulator;

// Definitions for convex hull
using Polyhedron = CGAL::Polyhedron_3<EPIC>;
using P_Facet_iterator = Polyhedron::Facet_iterator;
using P_Facet_const_iterator = Polyhedron::Facet_const_iterator;
using P_Edge_iterator = Polyhedron::Edge_iterator;
using P_Edge_const_iterator = Polyhedron::Edge_const_iterator;
using P_Halfedge_around_facet_circulator =
    Polyhedron::Halfedge_around_facet_circulator;
using P_Halfedge_around_facet_const_circulator =
    Polyhedron::Halfedge_around_facet_const_circulator;
using P_Vertex_iterator = Polyhedron::Vertex_iterator;
using P_Vertex_const_iterator = Polyhedron::Vertex_const_iterator;

// Miscellaneous definitions
using Object = CGAL::Object;
using Triangle = CGAL::Triangle_3<EPIC>;
using Triang_Vector = std::vector<Triangle>;
using Segment = EPIC::Segment_3;
using Tetrahedron = CGAL::Tetrahedron_3<EPIC>;
using Tetra_Vector = std::vector<Tetrahedron>;
using Box =
    CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Finite_cells_iterator>;

#endif
