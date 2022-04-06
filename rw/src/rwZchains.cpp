#include "rwchains.h"
// for absolute value of integer
#include <cstdlib>

// Store the chains and weights in a std::map

RWZchains::RWZchains(ST st, Chain c0)
{
    stree = st;
    chain.simplices = c0.simplices;
    chain.weights = c0.weights;
    coface_dictionary = {};
}

Chain RWZchains::updateChain()
{
    std::list<ST::Simplex_handle> candidates; // only need a forward list anyway
    candidates = {};
    // find the candidate simplexes, may be duplicates that put some weights on cofaces
    for (auto e : (this->chain).simplices)
    {
        std::list<ST::Simplex_handle> cofaces_e = {};

        // beware lot of recopies of lists
        if (this->coface_dictionary.count(e))
        {
            std::copy(this->coface_dictionary[e].begin(), this->coface_dictionary[e].end(), std::back_inserter(cofaces_e));
        }
        else
        {
            // look for all the cofaces
            auto cofaces = this->stree.cofaces_simplex_range(e, 1);
            this->coface_dictionary.insert(std::make_pair(e, cofaces));
            std::copy(cofaces.begin(), cofaces.end(), std::back_inserter(cofaces_e));
        }
        candidates.insert(candidates.begin(), cofaces_e.begin(), cofaces_e.end());
    }
    // std::cout << candidates.size() << "\n";
    // construct the next chain
    Chain next_chain;
    next_chain.simplices.assign(this->chain.simplices.begin(), this->chain.simplices.end());
    next_chain.weights.assign(this->chain.weights.begin(), this->chain.weights.end());

    // std::cout << next_chain.simplices.size() << "\n";

    if (!candidates.empty()) 
    {
        std::random_device seed; // seed for generator (engine)
        // Standard mersenne twister_engine seeded with rd()
        std::mt19937 engine(seed());
        // number distribution uniform among the candidates
        std::uniform_int_distribution<> distrib(0, candidates.size() - 1);
        auto it = candidates.begin();
        std::next(it, distrib(engine));
        // coefficient for the magnitude of the update of the chain
        int lbdamin = *std::min_element(this->chain.weights.begin(), this->chain.weights.end());

        bool orient = false; // coefficient to alternate (not included)
        std::size_t index_b = 0;
        int j;
        bool isInChain;
        for (const auto &b : this->stree.boundary_simplex_range(*it))
        {
            orient = !orient; // alternate in the sum
            isInChain = false;
            // stree.assign_key(b, index_b);
            // Modify the weight by 1 or -1 depending on the orientation in the boundary
            j = 0;
            for (auto sh : this->chain.simplices)
            {
                // if (stree.key(sh) == index_b)
                if (sh == b)
                {
                    std::cout << "In chain \n";
                    isInChain = true;
                    *std::next(next_chain.weights.begin(), j)  -= (2 * orient - 1) * lbdamin;
                    break;
                }
                j++;
            }
            // Add the simplex and weight to the list otherwise
            if (!isInChain)
            {
                next_chain.simplices.push_back(b);
                next_chain.weights.push_back((1 - 2 * orient) * lbdamin);
            }
            index_b++;
        }
    }
    // the chain contains all the visited simplices with their weights (can be 0)
    return next_chain;
}

int RWZchains::computeEnergy(Chain c)
{
    int energy = 0;
    for (const auto &w : c.weights)
    {
        energy = energy + abs(w);
    }
    return energy;
}

void RWZchains::writeChain(Chain c, std::ofstream &file)
{
    file << "[";
    for (auto s : c.simplices)
    {
        file << "( ";
        for (auto vertex : this->stree.simplex_vertex_range(s))
        {
            file << vertex << " "; // space as separator
        }
        file << ")";
    }

    file << "] | ";
    file << "[";
    for (auto w : c.weights)
    {
        file << w << " ";
    }
    file << "]" << std::endl;
}

void RWZchains::run(int n_steps, std::string name_file)
{
    std::ofstream file(name_file);
    writeChain(this->chain, file);
    for (int k = 0; k < n_steps; k++)
    {
        Chain next_chain = updateChain();
        this->chain.simplices.assign(next_chain.simplices.begin(), next_chain.simplices.end());
        this->chain.weights.assign(next_chain.weights.begin(), next_chain.weights.end());
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
    while (T > 1)
    {
        Chain next_chain = updateChain();
        // symmetric difference of the transition
        int delta_U = diffWeight(next_chain); // can be another type of energy gap
        // // std::cout << "DeltaU = " << delta_U << "\n";
        float prob = 1;
        if (delta_U >= 0)
        {
            prob = std::exp(-delta_U / T);
        }
        float u = (float)std::rand() / (float)RAND_MAX;
        if (prob < u)
        {
            this->chain.simplices = next_chain.simplices;
            this->chain.weights = next_chain.weights;
            writeChain(this->chain, file);
        }
        T = T * alpha; // T = T0 / std::log(m)
    }
    file.close();
}

// A way to see the gap of energy
int RWZchains::diffWeight(Chain c)
{
    // compute the difference of weights
    return computeEnergy(c) - computeEnergy(this->chain); // integer energy exceptionally
}
