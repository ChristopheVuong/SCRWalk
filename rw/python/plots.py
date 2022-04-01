"""Some snippets of code for plotting simplicial complexes and animation of random walks."""
# Author: Christophe Vuong <christophe.vuong@telecom-paris.fr>
# License: BSD 3 clause
import matplotlib.collections as clt
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

import networkx as nx


def plot_2dim_tree(
    tree,
    pts_array,
    l_paths=None,
    labelling=False,
    style="fill",
    title="",
    verbose=False,
    margin=1.0,
):
    """Plot 2-simplicial complex based on Gudhi simplex tree (can be easily extended with 3DCollection) along
    with animation of paths if indicated
    Use Matplotlib that may be slow for this task"""
    # We need to set the plot limits.
    fig, ax = plt.subplots()
    ax.set_xlim(pts_array[:, 0].min() - margin, pts_array[:, 0].max() + margin)
    ax.set_ylim(pts_array[:, 1].min() - margin, pts_array[:, 1].max() + margin)
    plt.axis("equal")
    plt.axis("off")
    # list of triangles (point indices)
    triangles = np.array([s[0] for s in tree.get_skeleton(2) if len(s[0]) == 3])
    # list of edges (point coordinates)
    edges = np.array([s[0] for s in tree.get_skeleton(1) if len(s[0]) == 2])

    plt.plot(*pts_array.T, "k.", markersize=10)
    if title != "":
        plt.title(title)
    # Plot graph
    if labelling:
        num_seq = range(pts_array.shape[0])
        for num in num_seq:
            label = f"{num}"
            plt.annotate(
                label,
                (
                    pts_array[num, 0] - 0.5 * pts_array[:, 0].min(),
                    pts_array[num, 1] + 0.8 * pts_array[:, 1].min(),
                ),
            )

    ax.add_collection(clt.LineCollection(pts_array[edges], linewidths=0.5, colors="k"))
    if style == "fill":
        ax.add_collection(
            clt.PolyCollection(pts_array[triangles], facecolors="slateblue", alpha=0.4)
        )
    else:
        ax.add_collection(
            clt.PolyCollection(
                pts_array[triangles],
                linewidths=1.0,
                edgecolors=style,
                facecolors="white",
            )
        )
    if l_paths is not None:
        path = l_paths[0]
        collection = clt.LineCollection(
            pts_array[np.array(path)], linewidths=0.9, colors="r"
        )
        ax.add_collection(collection)
        len_path = len(l_paths)
        if len_path > 1:

            def animate(frame):
                collection.set_paths(pts_array[np.array(l_paths[frame])])
                if verbose:
                    ax.set_title("Step n°%d" % frame)
                return (collection,)

            anim = animation.FuncAnimation(
                fig, animate, frames=len_path, interval=800, blit=(not verbose)
            )
            return anim


def plot_2dim_complex(
    edges,
    triangles,
    pts_array,
    l_paths=None,
    labelling=False,
    style="fill",
    title="",
    verbose=False,
    margin=1.0,
):
    """Plot a 2-simplicial complex (can be easily extended with 3DCollection) along
    with animation of paths if indicated
    Use Matplotlib that may be slow for this task"""
    # We need to set the plot limits.
    fig, ax = plt.subplots()
    ax.set_xlim(pts_array[:, 0].min() - margin, pts_array[:, 0].max() + margin)
    ax.set_ylim(pts_array[:, 1].min() - margin, pts_array[:, 1].max() + margin)
    plt.axis("equal")
    plt.axis("off")

    plt.plot(*pts_array.T, "k.", markersize=10)
    if title != "":
        plt.title(title)
    # Plot graph
    if labelling:
        num_seq = range(pts_array.shape[0])
        for num in num_seq:
            label = f"{num}"
            plt.annotate(
                label,
                (
                    pts_array[num, 0] - 0.5 * pts_array[:, 0].min(),
                    pts_array[num, 1] + 0.8 * pts_array[:, 1].min(),
                ),
            )

    ax.add_collection(clt.LineCollection(pts_array[edges], linewidths=0.5, colors="k"))
    if style == "fill":
        ax.add_collection(
            clt.PolyCollection(pts_array[triangles], facecolors="slateblue", alpha=0.4)
        )
    else:
        ax.add_collection(
            clt.PolyCollection(
                pts_array[triangles],
                linewidths=1.0,
                edgecolors=style,
                facecolors="white",
            )
        )
    if l_paths is not None:
        path = l_paths[0]
        collection = clt.LineCollection(
            pts_array[np.array(path)], linewidths=0.9, colors="r"
        )
        ax.add_collection(collection)
        len_path = len(l_paths)
        if len_path > 1:

            def animate(frame):
                collection.set_paths(pts_array[np.array(l_paths[frame])])
                if verbose:
                    ax.set_title("Step n°%d" % frame)
                return (collection,)

            anim = animation.FuncAnimation(
                fig, animate, frames=len_path, interval=800, blit=(not verbose)
            )
            return anim


def plot_2dim_trees(
    l_trees,
    l_pts_array,
    l_paths=None,
    l_titles_txt=[],
    labelling=False,
    style="fill",
    title="",
    verbose=False,
    margin=1.0,
):
    """Plot several 2-simplicial complex based on Gudhi simplex tree along
    with animation of paths if indicated
    Use Matplotlib that may be slow for this task"""
    # We need to set the plot limits.
    fig, axs = plt.subplots(nrows=1, ncols=len(l_trees))
    for i in range(len(l_trees)):
        axs[i].set_xlim(
            l_pts_array[i][:, 0].min() - margin, l_pts_array[i][:, 0].max() + margin
        )
        axs[i].set_ylim(
            l_pts_array[i][:, 1].min() - margin, l_pts_array[i][:, 1].max() + margin
        )
        axs[i].axis("equal")
        axs[i].axis("off")
        # list of triangles (point indices)
        triangles = np.array(
            [s[0] for s in l_trees[i].get_skeleton(2) if len(s[0]) == 3]
        )
        # list of edges (point coordinates)
        edges = np.array([s[0] for s in l_trees[i].get_skeleton(1) if len(s[0]) == 2])

        axs[i].plot(*l_pts_array[i].T, "k.", markersize=10)
        # Plot graph
        if labelling:
            num_seq = range(l_pts_array[i].shape[0])
            for num in num_seq:
                label = f"{num}"
                axs[i].annotate(
                    label,
                    (
                        l_pts_array[i][num, 0] - 0.5 * l_pts_array[i][:, 0].min(),
                        l_pts_array[i][num, 1] + 0.8 * l_pts_array[i][:, 1].min(),
                    ),
                )

        axs[i].add_collection(
            clt.LineCollection(l_pts_array[i][edges], linewidths=0.5, colors="k")
        )
        if style == "fill":
            axs[i].add_collection(
                clt.PolyCollection(
                    l_pts_array[i][triangles], facecolors="slateblue", alpha=0.4
                )
            )
        else:
            axs[i].add_collection(
                clt.PolyCollection(
                    l_pts_array[i][triangles],
                    linewidths=1.0,
                    edgecolors=style,
                    facecolors="white",
                )
            )
        if l_titles_txt:
            axs[i].title.set_text(l_titles_txt[i])
    if title != "":
        plt.suptitle(title)
    if l_paths is not None:
        l_collections = []
        for i in range(len(l_paths)):
            path = l_paths[i][0]
            collection = clt.LineCollection(
                l_pts_array[i][np.array(path)], linewidths=0.9, colors="r"
            )
            l_collections.append(collection)
            axs[i].add_collection(collection)
        len_path = len(l_paths[0])

        if len_path > 1:

            def animate(frame):
                for i in range(len(l_paths)):
                    l_collections[i].set_paths(
                        l_pts_array[i][np.array(l_paths[i][frame])]
                    )
                    if verbose:
                        axs[i].set_title("Step n°%d" % frame)
                return l_collections

            anim = animation.FuncAnimation(
                fig, animate, frames=len_path, interval=800, blit=(not verbose)
            )
            return anim

def plot_2dim_complexes(
    l_edges,
    l_triangles,
    l_pts_array,
    l_paths=None,
    l_titles_txt=[],
    labelling=False,
    style="fill",
    title="",
    verbose=False,
    margin=1.0,
):
    """Plot several 2-simplicial complex along
    with animation of paths if indicated
    Use Matplotlib that may be slow for this task"""
    # We need to set the plot limits.
    fig, axs = plt.subplots(nrows=1, ncols=len(l_edges))
    for i in range(len(l_edges)):
        axs[i].set_xlim(
            l_pts_array[i][:, 0].min() - margin, l_pts_array[i][:, 0].max() + margin
        )
        axs[i].set_ylim(
            l_pts_array[i][:, 1].min() - margin, l_pts_array[i][:, 1].max() + margin
        )
        axs[i].axis("equal")
        axs[i].axis("off")
        # list of triangles (point indices)
        triangles = np.array(l_triangles[i])
        # list of edges (point coordinates)
        edges = np.array(l_edges[i])

        axs[i].plot(*l_pts_array[i].T, "k.", markersize=10)
        # Plot graph
        if labelling:
            num_seq = range(l_pts_array[i].shape[0])
            for num in num_seq:
                label = f"{num}"
                axs[i].annotate(
                    label,
                    (
                        l_pts_array[i][num, 0] - 0.5 * l_pts_array[i][:, 0].min(),
                        l_pts_array[i][num, 1] + 0.8 * l_pts_array[i][:, 1].min(),
                    ),
                )

        axs[i].add_collection(
            clt.LineCollection(l_pts_array[i][edges], linewidths=0.5, colors="k")
        )
        if style == "fill":
            axs[i].add_collection(
                clt.PolyCollection(
                    l_pts_array[i][triangles], facecolors="slateblue", alpha=0.4
                )
            )
        else:
            axs[i].add_collection(
                clt.PolyCollection(
                    l_pts_array[i][triangles],
                    linewidths=1.0,
                    edgecolors=style,
                    facecolors="white",
                )
            )
        if l_titles_txt:
            axs[i].title.set_text(l_titles_txt[i])
    if title != "":
        plt.suptitle(title)
    if l_paths is not None:
        l_collections = []
        for i in range(len(l_paths)):
            path = l_paths[i][0]
            collection = clt.LineCollection(
                l_pts_array[i][np.array(path)], linewidths=0.9, colors="r"
            )
            l_collections.append(collection)
            axs[i].add_collection(collection)
        len_path = len(l_paths[0])

        if len_path > 1:

            def animate(frame):
                for i in range(len(l_paths)):
                    l_collections[i].set_paths(
                        l_pts_array[i][np.array(l_paths[i][frame])]
                    )
                    if verbose:
                        axs[i].set_title("Step n°%d" % frame)
                return l_collections

            anim = animation.FuncAnimation(
                fig, animate, frames=len_path, interval=800, blit=(not verbose)
            )
            return anim



def save_animation(my_anim, name="random-walk-chains", writer="imagemagick", fps=30):
    """Save an animation as a GIF"""
    my_anim.save(name + ".gif", writer=writer, fps=fps)



def plot_2dim_weights_complex(
    points, aretes, weights, title="", name_="figure.png", save=False
):
    """Procedure for Zhihan's code for weight coloring of edges"""
    G = nx.Graph()
    for i in range(len(points)):
        G.add_node(i, pos=(points[i][0], points[i][1]))
    for j in range(len(weights)):
        G.add_edge(aretes[j][0], aretes[j][1], weight=weights[j])
    pos = nx.get_node_attributes(G, "pos")
    edges12, weights12 = zip(*nx.get_edge_attributes(G, "weight").items())
    nodes32 = nx.draw_networkx_nodes(
        G, pos, node_color="k", node_size=1, with_labels=False
    )
    #    Widths = tuple(10*(xyy-min(weights12))/(max(weights12)-min(weights12)) for xyy in weights12)
    Widths = tuple(10 * (xyy - 0) / (max(weights12) - 0) for xyy in weights12)
    cmap = plt.cm.Blues
    edges32 = nx.draw_networkx_edges(
        G, pos, edge_color=weights12, width=Widths, edge_cmap=cmap
    )  # ,edge_cmap=plt.cm.Greys
    # edges32 = nx.draw_networkx_edges(G,pos,edge_color=weights12,width=10*(weights12-min(weights12))/(max(weights12)-min(weights12)),edge_cmap=plt.cm.Greys)
    #    edges32 = nx.draw_networkx_edges(G,pos,edge_color=weights12,width=10*(max(weights12)-weights12)/(max(weights12)-min(weights12)),edge_cmap=plt.cm.Greys)

    # plt.colorbar(edges32)
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max(weights12))
    )
    sm._A = []
    plt.axis("equal")
    if title != "":
        plt.suptitle(title)
    plt.colorbar(sm)
    if save:
        plt.savefig(name_ + ".pdf", dpi=400)


def plot_2dim_weights_complexes(
    l_points,
    l_aretes,
    l_weights,
    l_titles_txt=[],
    title="",
    name_="figure.png",
    save=False,
):
    """Procedure for Zhihan's code for weight coloring of edges"""
    fig, axs = plt.subplots(nrows=1, ncols=len(l_aretes))
    for k in range(len(l_aretes)):
        G = nx.Graph()
        for i in range(len(l_points[k])):
            G.add_node(i, pos=(l_points[k][i][0], l_points[k][i][1]))
        for j in range(len(l_weights[k])):
            G.add_edge(l_aretes[k][j][0], l_aretes[k][j][1], weight=l_weights[k][j])
        pos = nx.get_node_attributes(G, "pos")
        edges12, weights12 = zip(*nx.get_edge_attributes(G, "weight").items())
        nodes32 = nx.draw_networkx_nodes(
            G, pos, node_color="k", node_size=1, with_labels=False, ax=axs[k]
        )
        #    Widths = tuple(10*(xyy-min(weights12))/(max(weights12)-min(weights12)) for xyy in weights12)
        Widths = tuple(10 * (xyy - 0) / (max(weights12) - 0) for xyy in weights12)
        cmap = plt.cm.Blues
        edges32 = nx.draw_networkx_edges(
            G, pos, edge_color=weights12, width=Widths, edge_cmap=cmap, ax=axs[k]
        )  # ,edge_cmap=plt.cm.Greys
        # edges32 = nx.draw_networkx_edges(G,pos,edge_color=weights12,width=10*(weights12-min(weights12))/(max(weights12)-min(weights12)),edge_cmap=plt.cm.Greys)
        #    edges32 = nx.draw_networkx_edges(G,pos,edge_color=weights12,width=10*(max(weights12)-weights12)/(max(weights12)-min(weights12)),edge_cmap=plt.cm.Greys)

        # plt.colorbar(edges32)
        sm = plt.cm.ScalarMappable(
            cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max(weights12))
        )
        sm._A = []
        axs[k].axis("equal")
        axs[k].axis("off")
        fig.colorbar(sm, ax=axs[k])
        axs[k].title.set_text(l_titles_txt[k])
    if title != "":
        plt.suptitle(title)
    if save:
        fig.savefig(name_ + ".pdf", dpi=400)
