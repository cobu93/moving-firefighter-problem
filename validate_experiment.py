import json
import os
import numpy as np
from config import RESULTS_DIR, EXPERIMENTS_FILE, RESULTS_FILE
from generators import tree_from_sequence

experiment_id = 2

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


print(experiment)
print(experiment_results)

tree = tree_from_sequence(np.array(experiment["sequence"]), add_positions=False)
tree.nodes_positions = np.array(np.concatenate((
                    experiment["nodes_positions"], 
                    [experiment["initial_firefighter_position"]]
                ), axis=0))
tree, _ = tree.to_directed(experiment["root"])
propagation_time = experiment["propagation_time"]
path_length = np.argwhere(tree.edges.sum(axis=-1) == 0).shape[0] + 1

print(path_length)

for m in experiment_results:
    print(f"============= {m}")
    path = experiment_results[m]["solution"]
    cum_dist = 0
    last_dist = 0

    print(path)

    t = 0
    for _ in range(len(path) - 1):
        current_time = (t + 1) * propagation_time
        print(f"........ Propagation {t}")

        cum_dist = last_dist

        while cum_dist < current_time:

            if(t > len(path) - 2):
                break

            dist = np.linalg.norm(tree.nodes_positions[path[t]] - tree.nodes_positions[path[t + 1]])
            cum_dist += dist

            if cum_dist < current_time:
                print(f"Distance from {path[t]} to {path[t+1]} id {dist}. Total distance: {cum_dist}")
                t += 1
            else:
                last_dist = cum_dist - dist
                break
            

