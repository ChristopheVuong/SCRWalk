#include <gudhi/pick_n_random_points.h>
 
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>
#include <iostream>
#include <vector>
#include <iterator>


#include <gudhi/Rips_complex.h>
// #include <gudhi/Sparse_rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <string>

#include <limits>  // for std::numeric_limits

// #include"cnpy.h"

#include "rwchains.h"

#include <chrono>

// #define CGAL_EIGEN3_ENABLED


int main(int argc, char **argv) {   

    // ----------------------------------------------------------------------------
    // Parameters of complex and paths to results
    // ----------------------------------------------------------------------------
    const int N_V = 500; // Number of vertices before sampling
    const std::string RESULT_FILENAME = "experiments/Rips_SA";
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

    typedef CGAL::Epick_d<CGAL::Dimension_tag<4> > K;
    typedef typename K::Point_d Point_d;
    
    CGAL::Random rd;
    
    std::vector<Point_d> points;
    for (int i = 0; i < N_V; ++i)
        points.push_back(Point_d(rd.get_double(-1., 1), rd.get_double(-1., 1),
                                rd.get_double(-1., 1), rd.get_double(-1., 1)));
    
    K k;
     //save it to file
    cnpy::npy_save("points.npy",&points[0],{N_V}, "w");
    // export points to a .xyz file (xy here), then process in Python with pandas 
    // std::vector<Point_d> results;
    // Gudhi::subsampling::pick_n_random_points(points, 100, std::back_inserter(results));
    // std::clog << "Before sparsification: " << points.size() << " points.\n";
    // std::clog << "After  sparsification: " << results.size() << " points.\n";
    
    // ----------------------------------------------------------------------------
    // Init a Rips complex from points
    // ----------------------------------------------------------------------------
    // double threshold = 12.0;
    double threshold = 0.5;
    Rips_complex rips_complex_from_points(points, threshold, Gudhi::Euclidean_distance());
    
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
    // ----------------------------------------------------------------------------
    // Init a random walk object
    // ----------------------------------------------------------------------------
    // find a subset of edges
    std::set<Simplex_tree::Simplex_handle> c0 = {};
    int num_edges = 0;
    for(auto e : stree.skeleton_simplex_range(1))
    {
        if ((stree.dimension(e) == 1) && (num_edges < 5))
        {
            c0.insert(e);
            num_edges++;
        }
    }
    // c0.insert(stree.find({3, 2})); 
    // c0.insert(stree.find({4, 3})); 
    // c0.insert(stree.find({6, 4})); 
    // c0.insert(stree.find({6, 2})); 
    //// Coface information
    // std::cout << "Initialization done with " << c0.size() << " edges!\n";  
    // for (auto s : c0) 
    // {
    //     std::clog << "   ( ";
    //     for (auto vertex : stree.simplex_vertex_range(s)) 
    //     {
    //         std::clog << vertex << " ";
    //     }

    //     std::clog << "Cofaces : ";
    //     auto cofaces = stree.cofaces_simplex_range(s, 1);
    //     for (auto cf : cofaces)
    //     {
    //         std::clog << "   ( ";
    //         for (auto vc : stree.simplex_vertex_range(cf))
    //         {
    //             std::clog << vc << " ";
    //         }
    //         std::clog << "   )";
    //     }
    //     std::clog << "Size = " << cofaces.size();
    //     std::clog << std::endl;
    // }
    // std::clog << std::endl; 

    // ----------------------------------------------------------------------------
    // Walk
    // ----------------------------------------------------------------------------
    // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();


    RWZ2chains rwZ2 (stree, c0);
    // rwZ2.run(N_STEPS, RESULT_FILENAME);
    rwZ2.runSA(T0, alpha, RESULT_FILENAME);
    // auto stop = std::chrono::high_resolution_clock::now();
    
    // Get duration. Substart timepoints to
    // get duration. To cast it to proper unit
    // use duration cast method
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
 
    // std::cout << "Time taken by the random walk: "
    //      << duration.count() << " microseconds" << std::endl;

    return 0;
}