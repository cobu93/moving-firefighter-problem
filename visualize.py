import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_tree(tree):
    nodes = tree.nodes_positions

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2])

    for e in tree.edges:
        o_x, o_y, o_z = nodes[e[0]]
        d_x, d_y, d_z = nodes[e[1]]
        ax.plot3D([o_x, d_x], [o_y, d_y], [o_z, d_z])

    plt.show()