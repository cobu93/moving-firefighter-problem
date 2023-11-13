import json
import os
import numpy as np
from config import RESULTS_DIR, EXPERIMENTS_FILE, RESULTS_FILE
from generators import tree_from_sequence

experiment_id = 1
experiments_file = os.path.join(RESULTS_DIR, EXPERIMENTS_FILE)

if not os.path.exists(experiments_file):
    raise ValueError("File {experiments_file} doesn't exists")

with open(experiments_file, "r") as f:
    experiments = json.load(f)

results_file = os.path.join(RESULTS_DIR, RESULTS_FILE)

if not os.path.exists(results_file):
    raise ValueError("File {experiments_file} doesn't exists")

with open(results_file, "r") as f:
    results = json.load(f)

experiment = None
for e in experiments:
    if e["id"] == experiment_id:
        experiment = e
        break

if experiment is None:
    raise ValueError("Experiment not found")

experiment_results = {}

for m in results:
    for e in results[m]:
        if e["experiment"] == experiment_id:
            experiment_results[m] = e
            break


tree = tree_from_sequence(np.array(experiment["sequence"]), add_positions=False)
tree.nodes_positions = np.array(np.concatenate((
                    experiment["nodes_positions"], 
                    [experiment["initial_firefighter_position"]]
                ), axis=0))

root = experiment["root"]
n_nodes = experiment["n_nodes"]

tree, _ = tree.to_directed(root)
propagation_time = experiment["propagation_time"]
path_length = np.argwhere(tree.edges.sum(axis=-1) == 0).shape[0]

print("Max path length:", path_length)

for m in experiment_results:
    print(f"============================== {m}")
    path = experiment_results[m]["solution"]
    forest = np.ones(n_nodes).astype(bool)
    forest[root] = False
    cum_dist = 0
    last_dist = 0
    saved_nodes = []
    saving_times = []

    print("Path:", path)
    
    for ff_t in range(path_length):

        if ff_t + 1 > len(path) - 1:
            print("Left movements: {}".format(path_length - ff_t))
            break

        path_to_root = tree.get_path_to_root(path[ff_t + 1])

        time_limit = path_to_root.shape[0] * propagation_time
        cum_dist += np.linalg.norm(tree.nodes_positions[path[ff_t]] - tree.nodes_positions[path[ff_t + 1]])

        in_forest = forest[path[ff_t + 1]]
        on_time = cum_dist < time_limit

        valid = in_forest and on_time

        subtree_nodes = tree.get_subtree_nodes(path[ff_t + 1])
        saved_nodes.append(len(subtree_nodes))
        saving_times.append(cum_dist)

        unavailable_nodes = path_to_root.tolist() + subtree_nodes.tolist()
        forest[unavailable_nodes] = False

        print(path[ff_t], "--", "✓" if valid else "x", "-->", path[ff_t + 1])

    print("Saving times:", saving_times)
    print("Saved nodes:", saved_nodes)
    print("Total saved nodes:", np.sum(saved_nodes))
    

