#ifndef CONTOURS_H
#define CONTOURS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Epick_d.h>

// convex hull headers
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>

// Alpha shape headers
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <fstream>
#include <iostream>
#include <list>

#include <vector>
#include <numeric>

// CGAL for points
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_2                                           Point_2;
// typedef K::Segment_d                                         Segment;

// Convex hull typedef
typedef CGAL::Convex_hull_traits_adapter_2<K,
          CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;

// typedef CGAL::Convex_hull_traits_adapter_2<K,
//           CGAL::Pointer_property_map<Point_d>::type > Convex_hull_traits_2;


typedef K::FT                                                FT;
// Alpha shape typedef
typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>          Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                 Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;

/**
 * @brief Compute the convex hull in out
 * @param out store the indices of the points in the convex hull. Edges?
 */
void get_vertices_convex_hull_2D(std::vector<Point_2> points, std::vector<std::size_t> out);

/**
 * @brief Find the alpha shape of a set of points and store the segments
 * 
 * @tparam OutputIterator : A list or a vector for example
 * @param A 
 * @param out store the segments
 */
template <class OutputIterator>
void alpha_edges( const Alpha_shape_2& A, OutputIterator out);


#endif // CONTOURS_H