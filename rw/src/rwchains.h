#ifndef RWCHAINS_H
#define RWCHAINS_H

#include <iostream>
#include <fstream>
#include <string>
// #include <limits>  // for numeric_limits
#include <algorithm>
#include <random>
#include <vector>
#include <list>
#include <iterator>
// #include <unordered_map>
#include <utility>
// #include <boost/program_options.hpp>
// #include <CGAL/Epick_d.h>
#include <boost/program_options.hpp>
#include <gudhi/Simplex_tree.h>

// struct MyOptions : Gudhi::Simplex_tree_options_full_featured
// {
//     // Not doing persistence, so we don't need those
//     static const bool store_key = false;
//     static const bool store_filtration = false;
//     // typedef short Vertex_handle; // short for now if we have few vertices, we'll see
// };

using ST = Gudhi::Simplex_tree<>;

class RWZ2chains
/**
   Class for the \f$ \mathbb{Z}/2\mathbb{Z} \f$ random walks that writes in a target file the visited chains
*/
{
private:
    ST stree;
    std::set<ST::Simplex_handle> chain;
    std::map<ST::Simplex_handle, std::vector<ST::Simplex_handle >> coface_dictionary; // map a simplex to the next possible simplex in the chain

public:
    
    /**
     * @brief Construct a new rw Z2chains object (use ST::Simplex_handle e = st.find(edge) to initialize of the chain);
     * 
     * @param stree the simplex tree associated to the simplicial complex
     * @param c0 the initial chain
     */
    RWZ2chains(ST st, std::set<ST::Simplex_handle> c0);
    
    /**
     * @brief Computes the transition to another chain.
     * @return the next chain boundary of a coface drawn uniformly
     */
    std::set<ST::Simplex_handle> makeTransition();
    
    /**
     * @brief Writes the chain in a text file to decode and plot in Python.
     */
    void writeChain(std::set<ST::Simplex_handle> c, std::ofstream &file);
    
    /**
     * @brief Runs the random walk given the initialization by the constructor and writes the progress in a file
     * 
     */
    void run(int n_steps, std::string);
    
    /**
     * @brief Runs a simulated annealing based on the random walk.
     * @param T0 the initial temperature
     * @param alpha the decay of the temperature at each step
     * @param name_file the name of the file where the chains are written at each step
     */
    void runSA(float T0, float alpha, std::string name_file);
    
    /**
     * @brief Computes the energy gap between chains.
     * @param c the chain to examine (for SA)
     * @return the difference of length between the current chain and the parameter c
     */
    int diffLength(std::set<ST::Simplex_handle>  c);
};

// /**
//  * @brief A structure for finite-length chains with integer weights
//  * 
//  */
// struct Chain 
// {
//     std::list<ST::Simplex_handle> simplices; // list of simplices in the chain
//     std::list<int> weights; // list of weights
// };



class RWZchains 
/**
   Class for the \f$ \mathbb{Z} \f$ random walks that writes in a target file the weights of each simplex in the chain.
   We avoid to store the weights of simplices as it is not necessary (only the non-zero weights) since in real-life applications we do not know the triangles of the simplicial complexes nor the edges, but we discover them with the random walk. It may not be appropriate since lists are not as easy to manage as in Python.
*/
{
private:
    ST stree;
    // Chain chain; 
    std::map<ST::Simplex_handle, long> chain;
    std::map<ST::Simplex_handle, long>  c0; // initial chain of simplices
    std::map<ST::Simplex_handle, std::vector<ST::Simplex_handle >> coface_dictionary; // map a simplex to the next possible simplex in the chain

public:
    
    RWZchains(ST st, std::map<ST::Simplex_handle, long> c0);

    /**
     * @brief Computes the energy of the chain
     * 
     * @param c 
     * @return int the energy
     */
    int computeEnergy(std::map<ST::Simplex_handle, long> c);
    
    /**
     * @brief Makes a transition in the simplex tree to another chain.
     * @return The possible chain to transition to
     */
    std::map<ST::Simplex_handle, long>  updateChain();
    
    /**
     * @brief Writes the chain in a text file to decode and plot in Python.
     */
    void writeChain(std::map<ST::Simplex_handle, long> chain, std::ofstream &file);
    
    /**
     * @brief Runs the random walk given the initialization by the constructor and writes the progress in a file
     * 
     */
    void run(int n_steps, std::string);
    
    /**
     * @brief Runs a simulated annealing based on the random walk.
     * @param T0 the initial temperature
     * @param alpha the decay of the temperature at each step
     * @param name_file the name of the file where the chains are written at each step
     */
    void runSA(float T0, float alpha, std::string name_file);
    
    /**
     * @brief Computes the energy gap between chains.
     * @param c the chain to examine (for SA)
     * @return the difference of length between the current chain and the parameter c
     */
    long diffWeight(std::map<ST::Simplex_handle, long>  c);
};

#endif //RWCHAINS_H