"""Plot the results of simulated annealing random walk from the C++ output files."""
# Author: Christophe Vuong <christophe.vuong@telecom-paris.fr>
# License: BSD 3 clause

import os, sys

p = os.path.abspath(__file__)
# print(p)
parent = os.path.dirname(p)
grandparent = os.path.dirname(parent)
# print(grandparent)
sys.path.insert(1, grandparent)

# Suppose that the C++ code has been compiled with CMake
path_to_experiments = grandparent + "/build/examples/"

# import gudhi
import matplotlib.pyplot as plt
from rw.python.utils import *


epsilon = 0.6

points = read_2Dpoints(path_to_experiments + "Points-Rips-SA")

chains = read_chains(path_to_experiments + "Rips-SA")
print(chains)

# R = gudhi.RipsComplex(points=points, max_edge_length=epsilon)
# s_tree = R.create_simplex_tree(max_dimension=2)


# plt.figure()
# plt.plot(rw0.get_lengths())
# plt.xlabel('Steps')
# plt.ylabel('Length of RW')
# plt.draw()

# chains0_to_display = tuple_to_list_paths(chains0)


# weights = rw0.get_weights(list_to_tuple(edges))
# plt.figure()
# plot_2dim_complex(s_tree1, np.array(points), [pathlist_init ], title=r'$\epsilon$ = %g' % epsilon)
# plt.draw()