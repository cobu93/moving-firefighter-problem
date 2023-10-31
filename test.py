from generators import generate_random_tree
from visualize import plot_tree
from utils import tree_to_structure, TREE
from ctypes import *

import numpy as np

tree = generate_random_tree(20)
#plot_tree(tree)

so_file = "./mfp.so"
mfp_functions = CDLL(so_file)

print("Python", tree.nodes.shape, tree.edges.shape)


root = np.random.choice(tree.nodes)
current_node = np.random.choice(tree.nodes)
tree_structure = tree_to_structure(tree)

print(root, current_node)

mfp_functions.mfp_dp_solver(
        tree_structure, 
        (c_int)(root),
        (c_int)(current_node)
    )



