from ctypes import *
import numpy as np

class Tree:
    def __init__(self, nodes, edges, nodes_positions=None):
        self.nodes = nodes
        self.edges = edges
        self.nodes_positions = np.array([nodes_positions[i] for i in self.nodes])


class TREE(Structure):
    _fields_ = [
        ("n_nodes", c_int),
        ("nodes", POINTER(c_int)),
        ("nodes_x", POINTER(c_float)),
        ("nodes_y", POINTER(c_float)),
        ("nodes_z", POINTER(c_float)),

        ("n_egdes", c_int),
        ("edges_o", POINTER(c_int)),
        ("edges_d", POINTER(c_int))
    ]


def tree_to_structure(tree):
    n_nodes = tree.nodes.shape[0]
    n_edges = tree.edges.shape[0]

    return TREE(
        n_nodes,
        (c_int * n_nodes)(*tree.nodes.tolist()),
        (c_float * n_nodes)(*tree.nodes_positions[:, 0].tolist()),
        (c_float * n_nodes)(*tree.nodes_positions[:, 1].tolist()),
        (c_float * n_nodes)(*tree.nodes_positions[:, 2].tolist()),

        n_edges, 
        (c_int * n_edges)(*tree.edges[:, 0].tolist()),
        (c_int * n_edges)(*tree.edges[:, 1].tolist())
        )
        

