#ifndef LAPLACIAN_H
#define LAPLACIAN_H


#include <gudhi/Simplex_tree.h>
using ST = Gudhi::Simplex_tree<>;

#include <utility>      // std::pair, std::make_pair
#include <cstdio>
#include "fmpz_matxx.h"
#include "fmpzxx.h"
#include "arithxx.h"


/**
 * @brief Fill the initialized boundary matrix with the coefficients given all the edges
 * 
 * @param stree the simplex tree from which we extract the boundaries
 * @param k the dimension of the boundary matrix
 * @param boundary the matrix to build
 * @return std::pair<int,int> the shape of the boundary matrix
 */
std::pair<int,int> compute_boundary_matrix(ST stree, const int k, fmpz_mat_t boundary);

/**
 * @brief 
 * 
 * @param stree the simplex tree from which we extract the boundaries
 * @param k the dimension of the boundary matrix
 * @param boundarykplus the boundary of order k+1 to build
 * @param boundaryk the boundary of order k to build
 * @return int 
 */
int compute_boundary_matrices(ST stree, const int k, fmpz_mat_t boundarykplus, fmpz_mat_t boundaryk);

/**
 * @brief Set the Laplacian of order k given the boundary matrices by matrix multiplication
 * 
 * @param boundaryk the computed boundary matrix of order k
 * @param boundarykplus the computed boundary matrix of order k+1
 * @param laplacian the matrix to store the Laplacian
 */
void set_laplaciank(fmpz_mat_t boundaryk, fmpz_mat_t boundarykplus, fmpz_mat_t laplacian);


/**
 * @brief Get the null space of laplacian of order 1
 * 
 * @param stree the simplex tree from which we get the boundaries
 * @param null_basis a matrix to store the null_space of the Laplacian
 */
void get_null_space_laplacian1(ST stree, fmpz_mat_t null_basis);

#endif //LAPLACIAN_H