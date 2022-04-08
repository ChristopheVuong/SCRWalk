#include <gudhi/pick_n_random_points.h>
 
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>


#include <gudhi/Rips_complex.h>
// #include <gudhi/Sparse_rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <string>

#include <limits>  // for std::numeric_limits

#include "rwchains.h"
#include "laplacian_flint.h"

#include <chrono>

// #define CGAL_EIGEN3_ENABLED


int main(int argc, char **argv) {   

    // ----------------------------------------------------------------------------
    // Parameters of complex and paths to results
    // ----------------------------------------------------------------------------
    const int N_V = 100; // Number of vertices before sampling
    const double HOLE_RADIUS = 0.3; // radius of the holes in the point clouds in [-1, 1]
    const std::string POINTS_FILENAME = "Points-Rips-SA-Lap";
    const std::string RESULT_FILENAME = "Rips-SA-Lap";
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
    
    std::vector<Point_d> points;

    std::ofstream file_points(POINTS_FILENAME);
    for (int i = 0; i < N_V; ++i)
    {
        // points.push_back(Point_d(rd.get_double(-1., 1), rd.get_double(-1., 1),
        //                         rd.get_double(-1., 1), rd.get_double(-1., 1)));
        double x = rd.get_double(-1., 1);
        double y = rd.get_double(-1., 1);
        // 2 holes: center (-0.5, -0.5) and center (0.5, 0.4)
        if (((x + 0.5)*(x + 0.5) + (y + 0.5)*(y + 0.5) >= HOLE_RADIUS*HOLE_RADIUS) && ((x - 0.5) * (x - 0.5) + (y - 0.4)*(y - 0.4) >= HOLE_RADIUS*HOLE_RADIUS))
        {
            // do not sample inside some predefined circles
            points.push_back(Point_d(x, y));
            file_points << x << " " << y << std::endl; 
        }
    }
    K k;
    file_points.close();

    // Gudhi::subsampling::pick_n_random_points(points, 100, std::back_inserter(results));
    // std::clog << "Before sparsification: " << points.size() << " points.\n";
    // std::clog << "After  sparsification: " << results.size() << " points.\n";
    
    // ----------------------------------------------------------------------------
    // Init a Rips complex from points
    // ----------------------------------------------------------------------------
    // double threshold = 12.0;
    double threshold = 0.6;
    Rips_complex rips_complex_from_points(points, threshold, Gudhi::Euclidean_distance());
    
    Simplex_tree stree;
    rips_complex_from_points.create_complex(stree, 2);
    
    // ----------------------------------------------------------------------------
    // Display information about the one skeleton Rips complex
    // ----------------------------------------------------------------------------
    std::clog << "Rips complex is of dimension " << stree.dimension() <<
                " - " << stree.num_simplices() << " simplices - " <<
                stree.num_vertices() << " vertices." << std::endl;

    // ----------------------------------------------------------------------------
    // Init a random walk object with null basis of Laplacian
    // ----------------------------------------------------------------------------
    // Compute the null space of the Laplacian
    auto start = std::chrono::high_resolution_clock::now();

    fmpz_mat_t null_basis;
    get_null_space_laplacian1(stree, null_basis); 
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
 
    std::cout << "Time taken for computing the Laplacian null space with Flint: "
         << duration.count() << " microseconds" << std::endl;
    
    // fmpz_mat_print_pretty(null_basis);
    // Choose one column of the matrix of null basis
    const int IDX_ELEMENT = 0;
    // Create a initial chain of edges
    std::map<ST::Simplex_handle, long> c0 = {};

    int idx = 0;
    for (auto sh : stree.skeleton_simplex_range(1))
    {
        if (stree.dimension(sh) == 1)
        {
            if (*fmpz_mat_entry(null_basis, idx, IDX_ELEMENT) != 0) 
            {
                c0.insert({sh, *fmpz_mat_entry(null_basis, idx, IDX_ELEMENT)});
            }
            idx ++; // beware exceed size of matrix
        }
    }
    // for (auto const &w: c.weights) {
    //     std::cout << w << std::endl;
    // }

    fmpz_mat_clear(null_basis);

    // Initialize the random walk with that element and the weights corresponding to the coefficients in the column 

    
    // ----------------------------------------------------------------------------
    // Walk
    // ----------------------------------------------------------------------------
    // Get starting timepoint
    start = std::chrono::high_resolution_clock::now();


    RWZchains rwZ (stree, c0);
    // rwZ.run(N_STEPS, RESULT_FILENAME);
    rwZ.runSA(T0, alpha, RESULT_FILENAME);
    stop = std::chrono::high_resolution_clock::now();
    
    // Get duration. Substart timepoints to
    // get duration. To cast it to proper unit
    // use duration cast method
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
 
    std::cout << "Time taken by the random walk: "
         << duration.count() << " microseconds" << std::endl;

    return 0;
}