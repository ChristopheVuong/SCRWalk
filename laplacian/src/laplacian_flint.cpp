/**
   A small module  to compute the Laplacian of a Gudhi simplex tree complex.
   Using the basic multidimensional array in C++ (may be an issue for space complexity)
   Can also give a shot to std::array<int, type size>


   Note: A sparse version might be suited for applications using Linbox.
*/
// #include <linbox/sparse-coo-matrix.h>
#include "laplacian_flint.h"

//// Without Linbox


// Compute the coboundary matrix

// Question: rather use the map interface of Gudhi at https://gudhi.inria.fr/doc/2.3.0/struct_filtered_complex.html

std::pair<int, int> compute_boundary_matrix(ST stree, const int k, fmpz_mat_t boundary)
{
   // Map the k-simplices and (k-1)-simplices to the set of integers from 0 to n_k and 0 to n_{k-1} respectively
   std::map<slong, ST::Simplex_handle> k_1_simplices;
   std::map<ST::Simplex_handle, slong> k_simplices;
   slong ncols = 0; // index of the k-simplex (column)
   slong nrows = 0; // index of the (k-1)-simplex (line)
   auto skeleton_range_k = stree.skeleton_simplex_range(k + 1);
   for (auto sh : skeleton_range_k)
   {
      if (stree.dimension(sh) == k + 1)
      {
         k_1_simplices.insert(std::pair<int, ST::Simplex_handle>(ncols, sh));
         ncols++; // postfix increment
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
   std::pair<int, int> shape(nrows, ncols);
   // Fill the matrix column by column
   for (slong j = 0; j < ncols; j++)
   {
      bool isNegative = false;
      for (auto s : stree.boundary_simplex_range(k_1_simplices[j]))
      {
         isNegative = !isNegative;
         slong i = k_simplices[s];
         fmpz_set_si(fmpz_mat_entry(boundary, i, j), 2 * isNegative - 1);
      }
   }
   return shape;
}

std::tuple<int, int, int> compute_boundary_matrices(ST stree, const int k, fmpz_mat_t boundarykplus, fmpz_mat_t boundaryk)
{
   typedef boost::bimap<slong, ST::Simplex_handle> bmapsh;
   // Enumerate the (k+2)-simplices, (k+1)-simplices and k-simplices
   std::map<slong, ST::Simplex_handle> kplus2_simplices;
   bmapsh kplus1_simplices;
   std::map<ST::Simplex_handle, slong> k_simplices;
   slong nkplus2 = 0; // index of the (k+2)-simplex (column)
   slong nkplus1 = 0; // index of the (k+1)-simplex (line)
   slong nk = 0;      // index of the k-simplex (line)
   auto skeleton_range_kplus2 = stree.skeleton_simplex_range(k + 2);
   for (auto sh : skeleton_range_kplus2)
   {
      if (stree.dimension(sh) == k + 2)
      {
         kplus2_simplices.insert(std::pair<int, ST::Simplex_handle>(nkplus2, sh));
         nkplus2++; // postfix increment
      }
      else if (stree.dimension(sh) == k + 1)
      {
         kplus1_simplices.insert({nkplus1, sh});
         nkplus1++; // postfix increment
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
   // std::cout << nk << " " << nkplus1 << " " << nkplus2;
   // Resize the boundary matrix according to the obtained size
   fmpz_mat_init(boundarykplus, nkplus1, nkplus2);
   fmpz_mat_init(boundaryk, nk, nkplus1);
   // Fill the 2 matrices column by column
   for (int j = 0; j < nkplus1; j++)
   {
      bool isNegativek = false;
      for (auto s : stree.boundary_simplex_range(kplus1_simplices.left.at(j)))
      {
         isNegativek = !isNegativek;
         slong i = k_simplices[s];
         fmpz_set_si(fmpz_mat_entry(boundaryk, i, j), (int)(2 * isNegativek - 1));
      }
   }

   for (int j = 0; j < nkplus2; j++)
   {
      bool isNegativekplus = false;
      for (auto s : stree.boundary_simplex_range(kplus2_simplices[j]))
      {
         isNegativekplus = !isNegativekplus;
         slong i = kplus1_simplices.right.at(s);
         fmpz_set_si(fmpz_mat_entry(boundarykplus, i, j), 2 * isNegativekplus - 1);
      }
   }
   return std::make_tuple(nk, nkplus1, nkplus2);
}

void set_laplaciank(fmpz_mat_t boundaryk, fmpz_mat_t boundarykplus, std::tuple<int, int, int> shape, fmpz_mat_t laplacian)
{
   fmpz_mat_t C1;
   fmpz_mat_t C2;
   fmpz_mat_t boundaryk_T;
   fmpz_mat_t boundarykplus_T;

   int n = std::get<1>(shape);
   fmpz_mat_init(laplacian, n, n);
   fmpz_mat_init(C1, n, n);
   fmpz_mat_init(C2, n, n);
   // allocate the shape for the transpose matrix
   // fmpz_mat_init(boundarykT, fmpz_mat_t::boundaryk.cols() , boundaryk.rows());
   fmpz_mat_init(boundaryk_T, n, std::get<0>(shape));
   fmpz_mat_transpose(boundaryk_T, boundaryk);
   fmpz_mat_mul(C1, boundaryk_T, boundaryk);
   fmpz_mat_init(boundarykplus_T, std::get<2>(shape), n);
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
   auto shape = compute_boundary_matrices(stree, 1, boundary2, boundary1);
   int n = std::get<1>(shape);
   // init the laplacian with the appropriate size
   // std::cout << "Ok boundaries !";
   // fmpz_mat_print_pretty(boundary2);
   set_laplaciank(boundary1, boundary2, shape, laplacian);
   fmpz_mat_init(null_basis, n, n); // enumeration over the k-th simplex
   fmpz_mat_nullspace(null_basis, laplacian);
   fmpz_mat_clear(boundary1);
   fmpz_mat_clear(boundary2);
   fmpz_mat_clear(laplacian);
   // return n; // the size of the laplacian matrix
}