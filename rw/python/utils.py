"""Some snippets for generating Python lists of result files generated in C++."""
# Author: Christophe Vuong <christophe.vuong@telecom-paris.fr>
# License: BSD 3 clause

import re
# import pandas as pd
import io
import numpy as np


def read_2Dpoints(filename):
    """Reads the 2D coordinates of points in a file"""
    with open(filename, "r") as f:
        lines = f.readlines()
    points = np.zeros((len(lines), 2))
    idx = 0
    for line in lines:
        point = re.split(" ", line)
        points[idx,0] = point[0]
        points[idx,1] = point[1]
        idx += 1
    return points

def read_edges(filename):
    """Reads a text file that contains the edges of a complex"""
    edges = []
    with open(filename, "r") as  f:
        lines = f.readlines()
    for line in lines:
        verts = [u for u in line[1:-2].split() if u!=''] # strip the trailing space and last character
        edges.append([int(v) for v in verts])
    return edges

def read_triangles(filename):
    """Reads a text file that contains the edges of a complex"""
    triangles = []
    with open(filename, "r") as  f:
        lines = f.readlines()
    for line in lines:
        verts = [u for u in line[1:-2].split() if u!=''] # strip the trailing space and last character
        triangles.append([int(v) for v in verts])
    return triangles

def read_chains(filename):
    """Reads a text file and transpose the chains into a list of lists of simplices"""
    chains = []
    with open(filename, "r") as f:
        lines = f.readlines()
    for line in lines:
        c = []
        line_strip = line[1:-2] # strip the trailing space and last character
        # print("My line after stripping: %s"  % line_strip)
        lchain = re.split(r"\(|\)", line_strip)
        lchain = [u for u in lchain if u!= '']
        # print(lchain)
        for s in lchain:
            verts = s.split()
            c.append([int(v) for v in verts])
        chains.append(c)
    return chains

def read_weighted_chains(filename):
    """Reads a text file and restitute the chains along with the weights of each edges in the chain"""
    chains = []
    weights = []
    with open(filename, "r") as f:
        lines = f.readlines()
    for line in lines:
        c = []
        w = []
        line_strip = line[1:-2] # strip the trailing space and last character
        chainsAndWeights = re.split("|")
        lchains = chainsAndWeights[0][:-2]
        lweights = chainsAndWeights[1][:-2]
        for s in lchains:
            verts = s.split()
            c.append([int(v) for v in verts])
        for u in lweights:
            wei = u.split()
            w.append([int(we) for we in wei])
        chains.append(c)
        weights.append(w)
    return chains, weights


# use the points to regenerate the complex from the Python side (less heavier in terms of text files)?


# def points_to_np_array(filename):
#     """Reads a .xyz file (exported with CGAL) and put the point coordinates in a numpy array"""
#     # read xyz file and convert to a dataframe
#     xyz_df = pd.read_csv(
#         io.StringIO(filename),
#         delim_whitespace=True,
#         skiprows=2,
#         names=["element", "x", "y", "z"],
#     )
#     # convert it to numpy array
#     return xyz_df.to_numpy() # xyz_df.values
