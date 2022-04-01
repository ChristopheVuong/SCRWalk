/**
   A small module  to compute the Laplacian of a Gudhi simplex tree complex.
   Using the basic multidimensional array in C++ (may be an issue for space complexity)
   Can also give a shot to std::array<int, type size> 
   

   Note: A sparse version might be suited for applications using Linbox.
*/
// #include <linbox/sparse-coo-matrix.h>

#include "laplacian.h"
#include <boost/bimap.hpp>


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



std::pair<int,int> compute_boundary_matrix(ST stree, const int k, fmpz_mat_t boundary)
{
   // Map the k-simplices and (k-1)-simplices to the set of integers from 0 to n_k and 0 to n_{k-1} respectively
   std::map<slong, ST::Simplex_handle> k_1_simplices;
   std::map<ST::Simplex_handle, slong> k_simplices;
   slong ncols = 0; // index of the k-simplex (column)
   slong nrows = 0; // index of the (k-1)-simplex (line)
   auto skeleton_range_k = stree.skeleton_simplex_range(k+1);
   for(auto sh : skeleton_range_k)
   {
      if (stree.dimension(sh) == k+1) 
      {
         k_1_simplices.insert(std::pair<int, ST::Simplex_handle>(ncols, sh));
         ncols++; //postfix increment
      }
      else 
      {
         if (stree.dimension(sh) == k) 
         {
            k_simplices.insert(std::pair<ST::Simplex_handle, int>(sh, nrows));
            nrows++;
         }
         
      }
   }
   // Resize the boundary matrix according to the obtained size
   fmpz_mat_init(boundary, nrows, ncols);
   std::pair <int,int> shape (nrows, ncols);
   // Fill the matrix column by column
   for (slong j = 0; j < ncols; j++) 
   {
      bool isNegative = false;
      for (auto s : stree.boundary_simplex_range(k_1_simplices[j]))
      {
         isNegative = !isNegative;
         slong i = k_simplices[s];
         fmpz* entry = fmpz_mat_entry(boundary, i, j);
         *entry = 2 * isNegative - 1;
      }

   }
   return shape;
   
}

int compute_boundary_matrices(ST stree, const int k, fmpz_mat_t boundarykplus, fmpz_mat_t boundaryk)
{
   typedef boost::bimap< slong, ST::Simplex_handle > bmapsh;
   // Enumerate the (k+2)-simplices, (k+1)-simplices and k-simplices
   std::map<slong, ST::Simplex_handle>   kplus2_simplices;
   bmapsh kplus1_simplices;
   std::map<ST::Simplex_handle, slong> k_simplices;
   slong nkplus2 = 0; // index of the (k+2)-simplex (column)
   slong nkplus1 = 0; // index of the (k+1)-simplex (line)
   slong nk = 0; // index of the k-simplex (line)
   auto skeleton_range_kplus2 = stree.skeleton_simplex_range(k+2);
   for(auto sh : skeleton_range_kplus2)
   {
      if (stree.dimension(sh) == k+2) 
      {
         kplus2_simplices.insert(std::pair<int,ST::Simplex_handle>(nkplus2, sh));
         nkplus2++; //postfix increment
      }
      else if (stree.dimension(sh) == k+1)
      {
         kplus1_simplices.insert({nkplus1, sh});
         nkplus1++; //postfix increment
      }
      else 
      {
         if (stree.dimension(sh) == k) 
         {
            k_simplices.insert(std::pair<ST::Simplex_handle, int>(sh, nk));
            nk++;
         }
         
      }
   }
   // Resize the boundary matrix according to the obtained size
   fmpz_mat_init(boundarykplus, nkplus1, nkplus2);
   fmpz_mat_init(boundaryk, nk, nkplus1);
   // Fill the first matrix line by line and the second one column by column
   for (slong j = 0; j < nkplus1; j++) 
   {
      bool isNegativek = false;
      for (auto s : stree.boundary_simplex_range(kplus1_simplices.left.at(j)))
      {
         isNegativek = !isNegativek;
         slong i = k_simplices[s];
         fmpz_set_ui(fmpz_mat_entry(boundaryk, i, j), 2 * isNegativek - 1);
      }

   }

   for (slong j = 0; j < nkplus2; j++) 
   {
      bool isNegativekplus = false;
      for (auto s : stree.boundary_simplex_range(kplus2_simplices[j]))
      {
         isNegativekplus = !isNegativekplus;
         slong i = kplus1_simplices.right.at(s);
         fmpz_set_ui(fmpz_mat_entry(boundaryk, i, j), 2 * isNegativekplus - 1);
      }

   }
   return nkplus1;
   
}

void set_laplaciank(fmpz_mat_t boundaryk, fmpz_mat_t boundarykplus, fmpz_mat_t laplacian)
{
   fmpz_mat_t C1;
   fmpz_mat_t C2;
   fmpz_mat_t boundaryk_T;
   fmpz_mat_t boundarykplus_T;
   fmpz_mat_transpose(boundaryk_T, boundaryk);
   fmpz_mat_mul(C1, boundaryk_T, boundaryk);
   fmpz_mat_transpose(boundarykplus_T, boundarykplus);
   fmpz_mat_mul(C2, boundarykplus, boundarykplus_T);
   fmpz_mat_add(laplacian, C1, C2);
}

// Get the null space of Laplacian matrix
void get_null_space_laplacian1(ST stree, fmpz_mat_t null_basis)
{
   fmpz_mat_t boundary1;
   fmpz_mat_t boundary2;
   fmpz_mat_t laplacian;
   int n = compute_boundary_matrices(stree, 1, boundary1, boundary2);
   // init the laplacian with the appropriate size
   // std::cout << "Ok boundaries !";
   fmpz_mat_print_pretty(boundary1);
   fmpz_mat_init(laplacian, n, n);
   set_laplaciank(boundary1, boundary2, laplacian);
   fmpz_mat_init(null_basis, n, n); // enumeration over the k-th simplex
   fmpz_mat_nullspace(null_basis, laplacian);
   fmpz_mat_clear(boundary1);
   fmpz_mat_clear(boundary2);
   fmpz_mat_clear(laplacian);

}