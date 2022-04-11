# include "contourEdges.h"

void get_vertices_convex_hull_2D(std::vector<Point_2> points, std::vector<std::size_t> out)
{
  std::vector<std::size_t> indices(points.size());
  std::iota(indices.begin(), indices.end(),0);
  std::cout << "Number of points = " << indices.size() << "\n"; 
  CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
                      Convex_hull_traits_2(CGAL::make_property_map(points)));
}

template <class OutputIterator>
void alpha_edges( const Alpha_shape_2& A, OutputIterator out)
{
  Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
                             end = A.alpha_shape_edges_end();
  for( ; it!=end; ++it)
    *out++ = A.segment(*it);
}

// template <class OutputIterator>
// bool file_input(OutputIterator out)
// {
//   std::ifstream is("data/fin", std::ios::in);
//   if(is.fail())
//   {
//     std::cerr << "unable to open file for input" << std::endl;
//     return false;
//   }
//   int n;
//   is >> n;
//   std::cout << "Reading " << n << " points from file" << std::endl;
//   std::copy_n(std::istream_iterator<Point_2>(is), n, out);
//   return true;
// }
// // Reads a list of points and returns a list of segments
// // corresponding to the Alpha shape.
// int main()
// {
//   std::list<Point> points;
//   if(! file_input(std::back_inserter(points)))
//     return -1;
//   Alpha_shape_2 A(points.begin(), points.end(),
//                   FT(10000),
//                   Alpha_shape_2::GENERAL);
//   std::vector<Segment> segments;
//   alpha_edges(A, std::back_inserter(segments));
//   std::cout << "Alpha Shape computed" << std::endl;
//   std::cout << segments.size() << " alpha shape edges" << std::endl;
//   std::cout << "Optimal alpha: " << *A.find_optimal_alpha(1)<<std::endl;
//   return 0;
// }