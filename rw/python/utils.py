"""Some snippets for generating Python lists of result files generated in C++."""
# Author: Christophe Vuong <christophe.vuong@telecom-paris.fr>
# License: BSD 3 clause

import re
import pandas as pd
import io


def read_chains(filename):
    """Reads a text file and transpose the chains into a list of lists of simplices"""
    chains = []
    with open(filename, "r") as f:
        lines = f.readlines()
    for line in lines:
        c = []
        line.strip("[")
        line.strip("]")
        lchain = re.split("( | )", line)
        for s in lchain:
            verts = s.split(" ")
            c.append([int(v) for v in verts])
        chains.append(c)
    return chains


# use the points to regenerate the complex from the Python side (less heavier in terms of text files)?


def points_to_np_array(filename):
    """Reads a .xyz file (exported with CGAL) and put the point coordinates in a numpy array"""
    # read xyz file and convert to a dataframe
    xyz_df = pd.read_csv(
        io.StringIO(filename),
        delim_whitespace=True,
        skiprows=2,
        names=["element", "x", "y", "z"],
    )
    # convert it to numpy array
    return xyz_df.to_numpy() # xyz_df.values
