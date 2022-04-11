#include "rwchains.h"
// for absolute value of integer
#include <cstdlib>

// Store the chains and weights in a std::map

RWZchains::RWZchains(ST st, std::map<ST::Simplex_handle,long> c0)
{
    stree = st;
    chain = c0;
    coface_dictionary = {};
}

std::map<ST::Simplex_handle, long> RWZchains::updateChain()
{
    std::list<ST::Simplex_handle> candidates; // only need a forward list anyway
    candidates = {};
    // find the candidate simplexes, may be duplicates that put some weights on cofaces
    for (auto c : this->chain)
    {
        std::vector<ST::Simplex_handle> cofaces_e;

        // beware lot of recopies of lists
        if (this->coface_dictionary.count(c.first) > 0)
        {
            std::copy(this->coface_dictionary[c.first].begin(), this->coface_dictionary[c.first].end(), std::back_inserter(cofaces_e));
        }
        else
        {
            // look for all the cofaces
            auto cofaces = this->stree.cofaces_simplex_range(c.first, 1);
            this->coface_dictionary.insert(std::make_pair(c.first, cofaces));
            std::copy(cofaces.begin(), cofaces.end(), std::back_inserter(cofaces_e));
        }
        candidates.insert(candidates.begin(), cofaces_e.begin(), cofaces_e.end());
    }
    // std::cout << candidates.size() << "\n";
    // construct the next chain
    std::map<ST::Simplex_handle, long> next_chain = this->chain;

    // std::cout << next_chain.simplices.size() << "\n";

    if (!candidates.empty()) 
    {
        // coefficient for the magnitude of the update of the chain
        auto it_lbdamin = std::min_element(std::begin(this->chain), std::end(this->chain), [](const auto& l, const auto& r) { return l.second < r.second; }); // minimum of of a map second argument, requires C++ 14
        long lbdamin = it_lbdamin->second; // the map is not empty
        
        std::random_device seed; // seed for generator (engine)
        // Standard mersenne twister_engine seeded with rd()
        std::mt19937 engine(seed());
        // number distribution uniform among the candidates
        std::uniform_int_distribution<> distrib(0, candidates.size() - 1);
        auto it = candidates.begin();
        int rn = distrib(engine);
        std::advance(it, rn);
        

        bool orient = false; // coefficient to alternate (not included)
        // bool isInChain;
        for (auto& b : this->stree.boundary_simplex_range(*it))
        {
            orient = !orient; // alternate in the sum
            // isInChain = false;
            // stree.assign_key(b, index_b);
            // Modify the weight by 1 or -1 depending on the orientation in the boundary
            if (next_chain.count(b) > 0)
            {
                next_chain.at(b) -= (2 * orient - 1) * lbdamin;
            }
            else 
            {
                next_chain.insert({b, (1 - 2 * orient) * lbdamin});
            }
        }
    }
    // the chain contains all the visited simplices with their weights (can be 0)
    return next_chain;
}

int RWZchains::computeEnergy(std::map<ST::Simplex_handle, long>  c)
{
    int energy = 0;
    for (auto w : c)
    {
        energy = energy + abs(w.second);
    }
    return energy;
}

void RWZchains::writeChain(std::map<ST::Simplex_handle, long>  c, std::ofstream &file)
{
    file << "[";
    for (auto s : c)
    {
        file << "( ";
        for (auto vertex : this->stree.simplex_vertex_range(s.first))
        {
            file << vertex << " "; // space as separator
        }
        file << ")";
    }

    file << "] | ";
    file << "[";
    for (auto w : c)
    {
        file << w.second << " ";
    }
    file << "]" << std::endl;
}

void RWZchains::run(int n_steps, std::string name_file)
{
    std::ofstream file(name_file);
    writeChain(this->chain, file);
    for (int k = 0; k < n_steps; k++)
    {
        std::map<ST::Simplex_handle, long>  next_chain = updateChain();
        this->chain = next_chain;
        writeChain(this->chain, file);
    }
    file.close();
}

// Simulated annealing random walk
void RWZchains::runSA(float T0, float alpha, std::string name_file)
{
    std::ofstream file(name_file);
    writeChain(this->chain, file);
    float T = T0;
    float prob;
    while (T > 1)
    {
        std::map<ST::Simplex_handle, long>  next_chain = updateChain();
        // symmetric difference of the transition
        int delta_U = diffWeight(next_chain); // can be another type of energy gap
        // // std::cout << "DeltaU = " << delta_U << "\n";
        prob = 1;
        if (delta_U >= 0)
        {
            prob = std::exp(-delta_U / T);
            // std::cout << "p = " << prob << "\n";
        }
        float u = (float)std::rand() / (float)RAND_MAX;
        if (u < prob)
        {
            this->chain = next_chain;
            writeChain(this->chain, file);
        }
        T *= alpha; // T = T0 / std::log(m)
    }
    file.close();
}

// A way to see the gap of energy
long RWZchains::diffWeight(std::map<ST::Simplex_handle, long>  c)
{
    // compute the difference of weights
    return computeEnergy(c) - computeEnergy(this->chain); // integer energy exceptionally
}
