from generators import generate_random_tree
from dp.runners import DynamicProgramming
from iqcp.runners import IQCP
from ilp.runners import ILP
import numpy as np
from config import N_NODES, N_DEGREES, N_SAMPLES, RESULTS_DIR, RESULTS_FILE, EXPERIMENTS_FILE
import json 
import os
import time

runners = {
    # Restricted sequence length to number of leaves
    # Pre-filled hops and subtrees
    # Accepts propagation time
    # Recovers the solution
    "dynamic_programming": {
            "runner": DynamicProgramming(),
            "max_nodes": 110
    },
    # Unrestricted sequence length to number of leaves
    # Pre-filled hops and subtrees
    # Accepts propagation time
    #"dynamic_programming_unrestricted_length": {
    #        "runner": DynamicProgramming(so_file="./dp/mfp_unrestricted_length.so"),
    #        "max_nodes": 15
    #},
    # Dynamic computing of hops and subtrees
    # Default propagation time to 1
    #"dynamic_programming_dynamic_computing": {
    #        "runner": DynamicProgramming(so_file="./dp/mfp_dynamic_computing.so"),
    #        "max_nodes": 15
    #},
    "iqcp": {
            "runner": IQCP(),
            "max_nodes": 110
    },
    "ilp": {
            "runner": ILP(),
            "max_nodes": 110
    },
}

n_nodes = np.array(N_NODES).astype(int)
n_degrees = np.array(N_DEGREES).astype(int)
n_samples = N_SAMPLES


results = {}
experiments = []


propagation_time = 1.

for n_i, n in enumerate(n_nodes):
    for d_i, d in enumerate(n_degrees):
        skipped = np.zeros(len(runners), dtype=bool)
        for s in range(n_samples):

            initial_ff_position = np.random.rand(3)
            tree, sequence = generate_random_tree(n, d)
            root = np.random.choice(tree.nodes)

            experiment_id = len(experiments) + 1

            experiments.append({
                "id": experiment_id,
                "n_nodes": tree.nodes.shape[0],
                "degree": int(d),
                "n_degrees": tree.edges.sum(axis=-1).tolist(),
                "sequence": sequence,
                "nodes_positions": tree.nodes_positions.tolist(),
                "root": int(root),
                "initial_firefighter_position": initial_ff_position.tolist(),
                "propagation_time": propagation_time
            })

            optimals = np.zeros(len(runners))
            
            for r_i, r in enumerate(runners):
                if n > runners[r]["max_nodes"]:
                    print("Runner execution skipped due to the number of nodes")
                    skipped[r_i] = True
                    continue

                skipped[r_i] = False
                start = time.time()
                optimal, solution = runners[r]["runner"].run(tree, root, initial_ff_position, propagation_time)
                end = time.time()

                if r not in results:
                    results[r] = []

                results[r].append({
                    "experiment": experiment_id,
                    "duration": end - start,
                    "solution": solution,
                    "optimal": optimal
                })

                optimals[r_i] = optimal

                non_skipped = np.argwhere(skipped[:r_i] == False).flatten()

                if not np.all(optimals[non_skipped] == optimal):

                    if not os.path.exists(RESULTS_DIR):
                        os.makedirs(RESULTS_DIR)

                    with open(os.path.join(RESULTS_DIR, RESULTS_FILE), "w") as outfile:
                        outfile.write(json.dumps(results, indent=4))

                    with open(os.path.join(RESULTS_DIR, EXPERIMENTS_FILE), "w") as outfile:
                        outfile.write(json.dumps(experiments, indent=4))
                                
                    raise ValueError(f"Runner {r} gaves an inconsistent result:", optimals)


if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

with open(os.path.join(RESULTS_DIR, RESULTS_FILE), "w") as outfile:
    outfile.write(json.dumps(results, indent=4))

with open(os.path.join(RESULTS_DIR, EXPERIMENTS_FILE), "w") as outfile:
    outfile.write(json.dumps(experiments, indent=4))
            
        



