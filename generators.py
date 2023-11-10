from utils import Tree
import networkx as nx
import numpy as np

def generate_prufer_sequence(n_nodes):
    return np.random.randint(0, n_nodes, size=n_nodes - 2)

def tree_from_sequence(sequence, add_positions=True):
    edges = []
    n_nodes = sequence.shape[0] + 2
    degrees = np.ones(n_nodes)


    counts = np.bincount(sequence, minlength=n_nodes)
    degrees += counts

    for a in sequence:
        origin = np.argwhere(degrees == 1)[0, 0]
        edges.append((a, origin))
        degrees[origin] -= 1
        degrees[a] -= 1

    remaining = np.argwhere(degrees == 1)[:, 0]
    assert remaining.shape[0] == 2, "There are more than 2 remaining degrees = 1"

    edges.append((remaining[1], remaining[0]))
    positions = None

    if add_positions:
        tree = nx.Graph()
        tree.add_nodes_from(np.arange(n_nodes))
        tree.add_edges_from(edges)
        positions = nx.drawing.layout.fruchterman_reingold_layout(tree, dim=3, scale=1.)
        
    return Tree(np.arange(n_nodes), np.array(edges), positions)

def generate_random_tree(n_nodes, root_degree, type_root_degree, add_positions=True, max_trials=1000):    
    
    sequence = None
    found_tree = False

    for i in range(max_trials):
        
        sequence = generate_prufer_sequence(n_nodes)
        counts = np.bincount(sequence, minlength=n_nodes)

        candidate_roots = None

        if type_root_degree == "exact":
            candidate_roots = np.argwhere(counts == root_degree - 1).flatten()
        elif type_root_degree == "min":
            candidate_roots = np.argwhere(counts >= root_degree - 1).flatten()
        else:
            raise ValueError(f"Root degree type {type_root_degree} not recongized!")

        if len(candidate_roots) > 0:
            found_tree = True
            break

    if not found_tree:
        raise ValueError(f"Can't find a tree of {n_nodes} nodes and a root of degree {root_degree} in {max_trials} trials.")


    root = np.random.choice(candidate_roots)
    return tree_from_sequence(sequence, add_positions), sequence.tolist(), int(root)

    