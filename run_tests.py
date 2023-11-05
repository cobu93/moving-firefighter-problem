from generators import generate_random_tree
from dp.runners import DynamicProgramming
from iqcp.runners import IQCP
from ilp.runners import ILP
import numpy as np
from config import MAX_NODES, MIN_NODES, N_SAMPLES, RESULTS_DIR, RESULTS_FILE, EXPERIMENTS_FILE
import json 
import os

runners = {
    # Restricted sequence length to number of leaves
    # Pre-filled hops and subtrees
    # Accepts propagation time
    "dynamic_programming": DynamicProgramming(),
    # Unrestricted sequence length to number of leaves
    # Pre-filled hops and subtrees
    # Accepts propagation time
    "dynamic_programming_unrestricted_length": DynamicProgramming(so_file="./dp/mfp_unrestricted_length.so"),
    # Dynamic computing of hops and subtrees
    # Default propagation time to 1
    "dynamic_programming_dynamic_computing": DynamicProgramming(so_file="./dp/mfp_dynamic_computing.so"),
    "iqcp": IQCP(),
    "ilp": ILP(),
}

n_nodes = np.arange(MIN_NODES, MAX_NODES + 1).astype(int)
n_samples = N_SAMPLES


results = {}
experiments = []

propagation_time = 1.

for n_i, n in enumerate(n_nodes):
    for s in range(n_samples):

        initial_ff_position = np.random.rand(3)
        tree, sequence = generate_random_tree(n)
        root = np.random.choice(tree.nodes)

        experiment_id = len(experiments) + 1

        experiments.append({
            "id": experiment_id,
            "n_nodes": tree.nodes.shape[0],
            "sequence": sequence,
            "root": int(root),
            "initial_firefighter_position": initial_ff_position.tolist(),
            "propagation_time": propagation_time
        })

        optimals = np.zeros(len(runners))

        for r_i, r in enumerate(runners):
            optimal, duration = runners[r].run(tree, root, initial_ff_position, propagation_time)

            if r not in results:
                results[r] = []

            results[r].append({
                "experiment": experiment_id,
                "duration": duration,
                "optimal": optimal
            })

            optimals[r_i] = optimal

            if not np.all(optimals[:r_i] == optimal):
                raise ValueError(f"Runner {r} gaves an inconsistent result:", optimals)


if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

with open(os.path.join(RESULTS_DIR, RESULTS_FILE), "w") as outfile:
    outfile.write(json.dumps(results, indent=4))

with open(os.path.join(RESULTS_DIR, EXPERIMENTS_FILE), "w") as outfile:
    outfile.write(json.dumps(experiments, indent=4))
            
        



