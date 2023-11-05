from utils import tree_to_structure
from ctypes import *
import numpy as np
import time

class DynamicProgramming:
    
    def __init__(self, so_file="./dp/mfp.so"):
        self.mfp_functions = CDLL(so_file)
        
    def run(self, i_tree, root, initial_ff_position, t_propagation):
        i_tree.add_firefighter_position(initial_ff_position)
        tree, _ = i_tree.to_directed(root)
        n_leaves = (int) (np.argwhere(tree.edges.sum(axis=-1) == 0).flatten().shape[0])
        tree_structure = tree_to_structure(tree)

        optimal = c_int()
        c_optimal_path = (c_int * n_leaves)()

        start = time.time()
        self.mfp_functions.mfp_dp_solver(
            tree_structure, 
            (c_int)(root),
            (c_int)(len(tree.nodes)),
            (c_float)(t_propagation),
            byref(optimal),
            c_optimal_path
        )

        optimal_path = [-1]

        for i in range(n_leaves):
            if c_optimal_path[i] != -1:
                optimal_path += [c_optimal_path[i]]

        end = time.time()

        return optimal.value, optimal_path, end - start
