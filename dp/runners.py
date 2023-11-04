from utils import tree_to_structure
from ctypes import *
import time

class DynamicProgramming:
    
    def __init__(self, so_file="./dp/mfp.so"):
        self.mfp_functions = CDLL(so_file)
        
    def run(self, tree, root, initial_ff_position, t_propagation):
        tree.add_firefighter_position(initial_ff_position)
        tree_structure = tree_to_structure(tree)
        optimal = c_int()

        start = time.time()
        self.mfp_functions.mfp_dp_solver(
            tree_structure, 
            (c_int)(root),
            (c_int)(len(tree.nodes)),
            (c_float)(t_propagation),
            byref(optimal)
        )
        end = time.time()

        return optimal.value, end - start
