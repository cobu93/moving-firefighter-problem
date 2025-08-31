# First set of experiments
# N_NODES = [20, 25, 40, 50, 75, 100]
# ROOT_DEGREE = 3
# RESULTS_FILE = "results_moving_nodes.json"
# EXPERIMENTS_FILE = "experiments_moving_nodes.json"

# Second set of experiments
# N_NODES = [40]
# ROOT_DEGREE = [3, 4, 5, 6, 7]
# RESULTS_FILE = "results_moving_roots.json"
# EXPERIMENTS_FILE = "experiments_moving_roots.json"

# N_SAMPLES = 10
# RESULTS_DIR = "results"
# 24 hours
# RUNNER_TIMEOUT_SEC = 86400

EXPERIMENTS = [
    # First set of experiments
    {
        "n_nodes": [20, 25, 40, 50, 75, 100],
        "root_degree": 3,
        "results_file": "results_moving_nodes.json",
        "experiments_file": "experiments_moving_nodes.json",
        "n_samples": 10, 
        "results_dir": "results",
        "runner_timeout_sec": 86400
    },
    # Second set of experiments
    {

        "n_nodes": [40],
        "root_degree": [3, 4, 5, 6, 7],
        "results_file": "results_moving_roots.json",
        "experiments_file": "experiments_moving_roots.json",
        "n_samples": 10, 
        "results_dir": "results",
        "runner_timeout_sec": 86400
    }
]