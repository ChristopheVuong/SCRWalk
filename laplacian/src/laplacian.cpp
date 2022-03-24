/**
   A small naive module (completely inefficient for now) to compute the Laplacian of a Gudhi simplex tree complex.
   Using the basic multidimensional array in C++ (may be an issue for space complexity)
   Can also give a shot to std::array<int, type size> 
   A sparse version might be suited for applications

   Update : Use Linbox sparse matrix that we fill in with the outputs of Gudhi by giving the set of edges, triangles, etc.
*/
// #include <linbox/sparse-coo-matrix.h>


#include <gudhi/Simplex_tree.h>
using ST = Gudhi::Simplex_tree<>;



//// Without Linbox
// /**
//  * @brief Computes k-th boundary matrix of a simplicial encoded in a Gudhi simplex tree
//  * 
//  * @param stree the simplex tree
//  * @param k the order
//  * @return int* a pointer to the boundary matrix
//  */
// int* compute_boundary_matrix(ST stree, int k);


// int* compute_boundary_matrix(ST stree, int k)
// {
//    // Map the k-simplices and (k-1)-simplices to the set of integers from 0 to n_k and 0 to n_{k-1} respectively
//    std::map<ST::Simplex_handle, int> k_simplices;
//    std::map<ST::Simplex_handle, int> k_1_simplices;
//    int j = 0; // index of the k-simplex (column)
//    int i = 0; // index of the (k-1)-simplex (line)
//    auto skeleton_range_k = stree.skeleton_simplex_range(k);
//    for(auto sh : skeleton_range_k)
//    {
//       if (stree.dimension(sh) == k) 
//       {
//          k_simplices.insert(std::pair<ST::Simplex_handle, int>(sh, j));
//          j++; //postfix increment
//       }
//       else 
//       {
//          if (stree.dimension(sh) == k-1) 
//          {
//             k_1_simplices.insert(std::pair<ST::Simplex_handle, int>(sh, i));
//             i++;
//          }
         
//       }
//    }
   
//    // Create an array for boundary matrix
//    int *boundary_matrix = new int[i][j]; // intialized at 0, pointer that cannot go out of scope
//    // Loop over every k-simplex in the simplex tree
//    int column = 0; 
//    for(auto sh : skeleton_range_k)
//    {
      
//       if (stree.dimension(sh) == k) 
//       {
//          bool orient = false; // sign to alternate
//          for (auto b : stree.boundary_simplex_range(sh))
//          {
//             orient = !orient;
//             int idx = k_simplices.find(b)->second; // among the (k-1)-simplices
//             *boundary_matrix[idx][column] = 2 * orient - 1;
//          }
//       }
      

//    }
//    // do not forger to delete[] boundary matrix once it is over after the computation
//    return boundary_matrix;
// }


// Givaro::ZRing<Integer> ZZ;


// LinBox::SparseMatrix<Givaro::ZRing<Integer>> compute_boundary_matrix_matrix(ST stree, int k);

// LinBox::SparseMatrix<Givaro::ZRing<Integer>> compute_boundary_matrix_matrix(ST stree, int k)
// {
//    // Map the k-simplices and (k-1)-simplices to the set of integers from 0 to n_k and 0 to n_{k-1} respectively
//    std::map<ST::Simplex_handle, int> k_simplices;
//    std::map<ST::Simplex_handle, int> k_1_simplices;
//    int j = 0; // index of the k-simplex (column)
//    int i = 0; // index of the (k-1)-simplex (line)
//    auto skeleton_range_k = stree.skeleton_simplex_range(k);
//    for(auto sh : skeleton_range_k)
//    {
//       if (stree.dimension(sh) == k) 
//       {
//          k_simplices.insert(std::pair<ST::Simplex_handle, int>(sh, j));
//          j++; //postfix increment
//       }
//       else 
//       {
//          if (stree.dimension(sh) == k-1) 
//          {
//             k_1_simplices.insert(std::pair<ST::Simplex_handle, int>(sh, i));
//             i++;
//          }
         
//       }
//    }
//    // Construct a fixed sparse array
//    LinBox::SparseMatrix<Givaro::ZRing<Integer>> B(ZZ, i, j);
//    // Fill the matrix with Gudhi information
//    return B;
// }



// Compute the coboundary matrix

// Get the null space of sparse matrix