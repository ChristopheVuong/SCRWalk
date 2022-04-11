#include <chrono>


#include <iostream>
#include <fstream>

#include <gudhi/pick_n_random_points.h>
 
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <vector>
#include <iterator>


#include <gudhi/Rips_complex.h>
// #include <gudhi/Sparse_rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <string>

#include <limits>  // for std::numeric_limits

#include "rwchains.h"
#include "contourEdges.h"


// #define CGAL_EIGEN3_ENABLED


int main(int argc, char **argv) {   

    // ----------------------------------------------------------------------------
    // Details of complex and paths to results
    // ----------------------------------------------------------------------------
    const int N_V = 300; // Number of vertices before sampling
    const double HOLE_RADIUS = 0.25; // radius of the holes in the point clouds in [-1, 1]
    // const double THRESHOLD = 0.6;
    const double THRESHOLD = 0.3;
  
    const std::string POINTS_FILENAME = "Points-Rips";
    const std::string EDGES_FILENAME = "Edges-Rips";
    const std::string TRIANGLES_FILENAME = "Triangles-Rips";
    const std::string RESULT_FILENAME = "Rips";
    const std::string SARESULT_FILENAME = "Rips-SA";
    // const std::string RESULT_FILENAME = "Rips_SA";
    // using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;    
    // using Point = std::vector<double>;

    // ----------------------------------------------------------------------------
    // Parameters of walk with or without the simulated annealing
    // ----------------------------------------------------------------------------
    const int N_STEPS = 10000; // Number of steps of the random walk
    const float T0 = 10000000;
    const float alpha = 0.998;
    // N_SA_steps = int(-log(T0) / log(alpha)) + 1, here equal to int(8050.98) + 1
    
    using Simplex_tree = Gudhi::Simplex_tree<>;


    // // ----------------------------------------------------------------------------
    // // Init Zhihan's simplex tree
    // // ----------------------------------------------------------------------------
    // Simplex_tree stree;
    // stree.insert_simplex_and_subfaces({4, 5, 6}, 0.);
    // stree.insert_simplex_and_subfaces({1, 5, 6}, 0.);
    // stree.insert_simplex_and_subfaces({1, 4, 6}, 0.);
    // stree.insert_simplex_and_subfaces({1, 4, 5}, 0.);
    // stree.insert_simplex_and_subfaces({3, 4, 6}, 0.);
    // stree.insert_simplex_and_subfaces({2, 3, 6}, 0.);
    // stree.insert_simplex_and_subfaces({0, 1}, 0.);
    // stree.insert_simplex_and_subfaces({0, 2}, 0.);
    // stree.insert_simplex_and_subfaces({0, 2, 7}, 0.);
    // stree.insert_simplex_and_subfaces({0, 1, 7}, 0.);
    // stree.insert_simplex_and_subfaces({1, 6, 7}, 0.);
    // stree.insert_simplex_and_subfaces({2, 6, 7}, 0.);
    // Generate a Rips complex
    using Filtration_value = Simplex_tree::Filtration_value;
    using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
    

    // ----------------------------------------------------------------------------
    // Init of a Rips complex from random points in [-1, 1]
    // ----------------------------------------------------------------------------
    // std::vector<Point> points;
    // points.push_back({1.0, 1.0});
    // points.push_back({7.0, 0.0});
    // points.push_back({4.0, 6.0});
    // points.push_back({9.0, 6.0});
    // points.push_back({0.0, 14.0});
    // points.push_back({2.0, 19.0});
    // points.push_back({9.0, 17.0});

    typedef CGAL::Epick_d<CGAL::Dimension_tag<2> > K;
    typedef typename K::Point_d Point_d;
    
    CGAL::Random rd;
    
    std::vector<Point_d> points; // for the Rips complex

    std::vector<Point_2> points2D; // for the convex hull

    std::ofstream file_points(POINTS_FILENAME); // rewrites existing content
    for (int i = 0; i < N_V; ++i)
    {
        // points.push_back(Point_d(rd.get_double(-1., 1), rd.get_double(-1., 1),
        //                         rd.get_double(-1., 1), rd.get_double(-1., 1)));
        double x = rd.get_double(-1., 1);
        double y = rd.get_double(-1., 1);
        // 2 holes: center (-0.5, -0.45) and center (0.35, 0.4)
        if (((x + 0.5)*(x + 0.5) + (y + 0.45)*(y + 0.45) >= HOLE_RADIUS*HOLE_RADIUS) && ((x - 0.35) * (x - 0.35) + (y - 0.4)*(y - 0.4) >= HOLE_RADIUS*HOLE_RADIUS))
        {
            // do not sample inside some predefined circles
            // points.push_back(Point_d(x, y));
            points.push_back(Point_d(x, y));
            points2D.push_back(Point_2(x, y)); // for the convex hull
            file_points << x << " " << y << std::endl; 
        }   
    }
    
    file_points.close();
    // K k;
    // std::vector<Point_d> results;
    // Gudhi::subsampling::pick_n_random_points(points, 100, std::back_inserter(results));
    // std::clog << "Before sparsification: " << points.size() << " points.\n";
    // std::clog << "After  sparsification: " << results.size() << " points.\n";
    
    // ----------------------------------------------------------------------------
    // Init a Rips complex from points
    // ----------------------------------------------------------------------------
    // double threshold = 12.0;
    Rips_complex rips_complex_from_points(points, THRESHOLD, Gudhi::Euclidean_distance());
    
    Simplex_tree stree;
    rips_complex_from_points.create_complex(stree, 2);
    
    // ----------------------------------------------------------------------------
    // Display information about the one skeleton Rips complex
    // ----------------------------------------------------------------------------
    std::clog << "Rips complex is of dimension " << stree.dimension() <<
                " - " << stree.num_simplices() << " simplices - " <<
                stree.num_vertices() << " vertices." << std::endl;

    // std::clog << "Iterator on Rips complex simplices in the filtration order, with [filtration value]:" <<
    //             std::endl;
    // for (auto f_simplex : stree.filtration_simplex_range()) {
    //     std::clog << "   ( ";
    //     for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
    //     std::clog << vertex << " ";
    //     }
    //     std::clog << ") -> " << "[" << stree.filtration(f_simplex) << "] ";
    //     std::clog << "Cofaces : ";
    //     for (auto cf : stree.cofaces_simplex_range(f_simplex, 1))
    //     {
    //         std::clog << "   ( ";
    //         for (auto vc : stree.simplex_vertex_range(cf))
    //         {
    //             std::clog << vc << " ";
    //         }
    //         std::clog << "   ) ";
    //     }
    //     std::clog << std::endl;
    // }


    // Then enumerate the edges that contains those points

    // ----------------------------------------------------------------------------
    // Store the edges and the triangles in a text file
    // ----------------------------------------------------------------------------
    std::ofstream file_edges(EDGES_FILENAME); // rewrites existing content
    std::ofstream file_triangles(TRIANGLES_FILENAME); // rewrites existing content
    for(auto sh : stree.skeleton_simplex_range(2))
    {
        if (stree.dimension(sh) == 1)
        {
            file_edges << "( ";
            
            for (auto vertex : stree.simplex_vertex_range(sh))
            {
                file_edges << vertex << " "; // space between the vertices
            }
            file_edges << ")" << std::endl;
        }


        if (stree.dimension(sh) == 2)
        {
            file_triangles << "( ";
            for (auto vertex : stree.simplex_vertex_range(sh))
            {
                file_triangles << vertex << " "; // space between the vertices
            }
            file_triangles << ")" << std::endl;
        }

    }
    file_edges.close();
    file_triangles.close();

    // ----------------------------------------------------------------------------
    // Init a random walk object
    // ----------------------------------------------------------------------------
    // find a subset of edges
    std::set<Simplex_tree::Simplex_handle> c0 = {};
    
    // ----------------------------------------------------------------------------
    // Init a random walk object with the intersection of convex hulls and edges
    // ----------------------------------------------------------------------------

    // std::vector<std::size_t> vertices_hull;
    // get_vertices_convex_hull_2D(points2D, vertices_hull);
    // std::cout << "The convex hull of size " << vertices_hull.size() << " contains vertices:"; 
    // for (auto v : vertices_hull) 
    // {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;
    // int inHull;
    // ----------------------------------------------------------------------------
    // Build the initial chain
    // ----------------------------------------------------------------------------
    int num_edges = 0;
    for(auto e : stree.skeleton_simplex_range(1))
    {
        if (num_edges == 5) 
        {
            break;
        }
        
        if ((stree.dimension(e) == 1))
        {
        //     inHull = 0;
        //     for (auto vertex : stree.simplex_vertex_range(e))
        //     {
        //         if (std::count(vertices_hull.begin(), vertices_hull.end(), vertex)) {
        //             inHull ++;
        //         }
        //     }
        //     if (inHull == 2)
        //     {
        //         c0.insert(e);
        //     }
            c0.insert(e);
            num_edges++;
        }
    }

    std::cout << c0.size() << " edges in the initial chain." << std::endl;

    // ----------------------------------------------------------------------------
    // Walk
    // ----------------------------------------------------------------------------
    // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();


    RWZ2chains rwZ2 (stree, c0);
    rwZ2.run(N_STEPS, RESULT_FILENAME);
    
    auto stop = std::chrono::high_resolution_clock::now();
    
    // Get duration
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
 
    std::cout << "Time taken by the random walk: "
         << duration.count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    
    RWZ2chains rwZ2SA (stree, c0);
    rwZ2SA.runSA(T0, alpha, SARESULT_FILENAME);

    stop = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by the SA random walk: "
         << duration.count() << " seconds" << std::endl;

    return 0;
}