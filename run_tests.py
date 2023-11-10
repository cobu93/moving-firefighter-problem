from generators import generate_random_tree, tree_from_sequence
from dp.runners import DynamicProgramming
from greedy.runners import Greedy
from iqcp.runners import IQCP
from ilp.runners import ILP
import numpy as np
from config import N_NODES, ROOT_DEGREE, N_SAMPLES, RESULTS_DIR, RESULTS_FILE, EXPERIMENTS_FILE, RUNNER_TIMEOUT_SEC
import json 
import os
import time
import multiprocessing

runners = {
    "dynamic_programming": {
            "runner": DynamicProgramming(),
            "max_nodes": 110,
            "validate_optimal": True
    },
    "greedy": {
            "runner": Greedy(),
            "max_nodes": 110,
            "validate_optimal": False
    },
    "iqcp": {
            "runner": IQCP(),
            "max_nodes": 110,
            "validate_optimal": True
    },
    "ilp": {
            "runner": ILP(),
            "max_nodes": 110,
            "validate_optimal": True
    },
}

queue = multiprocessing.Queue()

def runner_wrapper(func, tree, root, initial_ff_position, propagation_time):
    optimal, solution  = func(tree, root, initial_ff_position, propagation_time)
    queue.put((optimal, solution))



n_nodes = np.array(N_NODES).astype(int)
root_degree = ROOT_DEGREE
type_root_degree = None

if isinstance(root_degree, list):
    root_degree = np.array(root_degree)
    type_root_degree = "exact"
elif isinstance(root_degree, int):
    root_degree = np.array([root_degree])
    type_root_degree = "min"


n_samples = N_SAMPLES


experiments = []
propagation_time = 1.


experiments_file = os.path.join(RESULTS_DIR, EXPERIMENTS_FILE)

# Generate experiments
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

existing_experiments = None
considered_experiments = {}
current_id = 0

if os.path.exists(experiments_file):
    with open(experiments_file, "r") as f:
        existing_experiments = json.load(f)
    
    for e in existing_experiments:
        considered_experiments[e["n_nodes"]] = considered_experiments.get(e["n_nodes"], {})
        considered_experiments[e["n_nodes"]][e["degree"]] = considered_experiments[e["n_nodes"]].get(e["degree"], [])
        considered_experiments[e["n_nodes"]][e["degree"]].append(e)

        if e["id"] > current_id:
            current_id = e["id"]

current_id += 1

for n in n_nodes:
    for d in root_degree:
        for s in range(n_samples):
            sample = True
            exp = None

            if n in considered_experiments and d in considered_experiments[n] and considered_experiments[n][d]:
                if len(considered_experiments[n][d]) > 0:
                    exp = considered_experiments[n][d].pop()
                    experiments.append(exp)
                    sample = False

            if not sample:
                print(f"Found experiment for Nodes:{n} Root degree:{d} [Id: {exp['id']}]")
                continue

            print(f"Sampling experiment for Nodes:{n} Root degree:{d}")

            initial_ff_position = np.random.rand(3)
            tree, sequence, root = generate_random_tree(n, d, type_root_degree)
            
            experiment_id = current_id
            current_id += 1

            experiments.append({
                "id": experiment_id,
                "n_nodes": tree.nodes.shape[0],
                "type_root_degree": type_root_degree,
                "expected_root_degree": int(d),
                "root_degree": int(tree.edges[root].sum()),
                "sequence": sequence,
                "nodes_positions": tree.nodes_positions.tolist(),
                "root": int(root),
                "initial_firefighter_position": initial_ff_position.tolist(),
                "propagation_time": propagation_time
            })


with open(experiments_file, "w") as outfile:
    outfile.write(json.dumps(experiments, indent=4))

results = {}
results_file = os.path.join(RESULTS_DIR, RESULTS_FILE)

if os.path.exists(results_file):
    with open(results_file, "r") as f:
        results = json.load(f)


consistent_experiments = np.array([True] * len(experiments))

for e_i, e in enumerate(experiments):
    optimals = np.zeros(len(runners)) - 1
    optimals_validate = np.array([False] * len(runners))
    optimals_execute = np.array([False] * len(runners))

    print("============= Experiment", e["id"])

    for m_i, m in enumerate(runners):
        optimals_execute[m_i] = (n <= runners[m]["max_nodes"])
        if optimals_execute[m_i]:
            optimals_validate[m_i] = runners[m]["validate_optimal"]
        else:
            optimals_validate[m_i] = False

        results[m] = results.get(m, [])

        for r in results[m]:
            if r["experiment"] == e["id"]:
                optimals[m_i] = r["optimal"]
                break

    all_experiments_executed = np.all(optimals[optimals_execute] >= 0)
    all_experiments_consistent = np.all(optimals[optimals_validate] == optimals[optimals_validate][0])

    if all_experiments_executed and all_experiments_consistent:
        print(f"Every method was executed and all of them are consistent :D")

    if not all_experiments_executed:

        tree = tree_from_sequence(np.array(e["sequence"]), add_positions=False)
        tree.nodes_positions = np.array(e["nodes_positions"])
        root = e["root"]
        initial_ff_position = np.array(e["initial_firefighter_position"])
        propagation_time = e["propagation_time"]


        for m_i, m in enumerate(runners):
            if optimals_execute[m_i] and optimals[m_i] < 0:

                print(f"Executing experiment for runner '{m}'")
    
                p = multiprocessing.Process(
                    target=runner_wrapper, 
                    args=(
                        runners[m]["runner"].run, 
                        tree, 
                        root, 
                        initial_ff_position, 
                        propagation_time
                    ),
                    name=m
                )

                try:
                    start = time.time()
                    #optimal, solution = runners[m]["runner"].run(tree, root, initial_ff_position, propagation_time)
                    p.start()
                    p.join(RUNNER_TIMEOUT_SEC)
                    end = time.time()

                    if p.is_alive():
                        p.terminate()
                        p.join()
                        raise TimeoutError("Runner exceded the timeout limit")
                    
                    optimal, solution = queue.get()

                    message = "Done!"
                except Exception as ex:
                    end = time.time()
                    optimal = None
                    solution = None
                    message = "ERROR: " + str(ex)
                    print(message)

                results[m].append({
                    "experiment": e["id"],
                    "duration": end - start,
                    "solution": solution,
                    "optimal": optimal,
                    "message": message
                })

                
                print(f"Saving results....")
                with open(results_file, "w") as outfile:
                    outfile.write(json.dumps(results, indent=4))
                print(f"Done")
                
                optimals[m_i] = optimal

        all_experiments_consistent = np.all(optimals[optimals_validate] == optimals[optimals_validate][0])

    if not all_experiments_consistent:
        print("WARNING: Not all experiments' results are consistent, please verify! D:")
        for m_i, m in enumerate(runners):
            print(f"Method: {m}   Optimal: {optimals[m_i]}   Execute?: {optimals_execute[m_i]}   Validate?: {optimals_validate[m_i]}")

        consistent_experiments[e_i] = False
        
print("============= Executions Finished")
if not np.all(consistent_experiments):
    print("WARNING: Not all experiments' results are consistent, please verify! D:")
    for e_i, e in enumerate(experiments):
        if not consistent_experiments[e_i]:
            print(f"Experiment {e['id']} inconsistent")
else:
    print("Everything was excellent! Don't worry :D")

