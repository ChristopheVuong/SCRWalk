// #include <boost/program_options.hpp>
#include "rwchains.h"

RWZ2chains::RWZ2chains(ST st, std::set<ST::Simplex_handle> c0)
{
    stree = st;
    chain = c0;
    coface_dictionary = {};
}

std::set<ST::Simplex_handle> RWZ2chains::makeTransition()
{
    std::list<ST::Simplex_handle> candidates; // only need a forward list anyway
    candidates = {};
    // find the candidate simplexes, may be duplicates that put some weights on cofaces
    std::vector<ST::Simplex_handle> cofaces_e;
    for (auto e : this->chain)
    {

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
        cofaces_e.clear();
    }
    if (!candidates.empty()) {
        // find occurrences in candidates
        std::random_device seed; // seed for generator (engine)
        // Standard mersenne twister_engine seeded with rd()
        std::mt19937 engine(seed());
        // number distribution uniform among the candidates
        std::uniform_int_distribution<> distrib(0, candidates.size() - 1);
        auto it = candidates.begin();
        int rn = distrib(engine);
        std::advance(it, rn);

        // for (auto v : stree.simplex_vertex_range(*it))
        // {
        //     std::clog << "(" << v << ")";
        //     std::clog << std::endl;
        // }
        std::set<ST::Simplex_handle> boundaries_set = std::set<ST::Simplex_handle>();
        auto boundaries = this->stree.boundary_simplex_range(*it);
        boundaries_set.insert(boundaries.begin(), boundaries.end());
        // do not care abound the orientation of the simplex in the boundary in Z/2Z (presence, absence)
        return boundaries_set;
    }

    else
    {
        return std::set<ST::Simplex_handle>();
    }
    
}

void RWZ2chains::writeChain(std::set<ST::Simplex_handle> c, std::ofstream &file)
{
    file << "[";
    for (auto s : c)
    {
        file << "( ";
        for (auto vertex : this->stree.simplex_vertex_range(s))
        {
            file << vertex << " "; // space between the vertices
        }
        file << ")";
    }

    file << "]" << std::endl;

}

void RWZ2chains::run(int n_steps, std::string name_file)
{
    std::ofstream file(name_file);
    writeChain(this->chain, file);
    for (int k = 0; k < n_steps; k++)
    {
        std::set<ST::Simplex_handle> tau_prime = makeTransition();
        // symmetric difference of the transition
        std::set<ST::Simplex_handle> out;
        // std::set_symmetric_difference(this->chain.begin(), this->chain.end(), next_chain.begin(), next_chain.end(), std::back_inserter(out)); // not working for set iterator inputs
        std::set_symmetric_difference(this->chain.begin(), this->chain.end(), tau_prime.begin(), tau_prime.end(), std::inserter(out, out.end()));
        this->chain = out;
        writeChain(this->chain, file);
        // std::cout << "Step " << k + 1 << "\n";
    }
    file.close();
}

// Simulated annealing random walk
void RWZ2chains::runSA(float T0, float alpha, std::string name_file)
{
    std::ofstream file(name_file);
    writeChain(this->chain, file);
    float T = T0;
    float prob;
    while (T > 1)
    {
        std::set<ST::Simplex_handle> out;
        std::set<ST::Simplex_handle> tau_prime = makeTransition();
        // std::set_symmetric_difference(this->chain.begin(), this->chain.end(), next_chain.begin(), next_chain.end(), std::back_inserter(out)); // not the right symmetric difference
        std::set_symmetric_difference(this->chain.begin(), this->chain.end(), tau_prime.begin(), tau_prime.end(), std::inserter(out, out.end()));
        int delta_U = diffLength(out); // can be another type of energy gap
        prob = 1;
        if (delta_U >= 0)
        {
            prob = std::exp(-delta_U / T);
        }
        float u = (float)std::rand() / (float)RAND_MAX;
        if (u < prob)
        {
            // symmetric difference between the current chain and the boundary of the chosen coface
            this->chain = out;
            
        }
        writeChain(this->chain, file); // write chain even if it remains the same, takes some time actually
        T *= alpha; // T = T0 / std::log(m) // cooling rate for probabilty convergence in SA
    }
    file.close();
}

// A way to see the gap of energy
int RWZ2chains::diffLength(std::set<ST::Simplex_handle> c)
{
    return c.size() - this->chain.size(); // integer energy exceptionally
}
