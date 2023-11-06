from utils import Tree
import networkx as nx
import numpy as np

def generate_prufer_sequence(n_nodes, mean_degree):
    sequence = np.repeat(np.arange(n_nodes), mean_degree - 1)
    return np.random.choice(sequence, size=n_nodes - 2)

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

def generate_random_tree(n_nodes, mean_degree, add_positions=True):    
    sequence = generate_prufer_sequence(n_nodes, mean_degree)
    return tree_from_sequence(sequence, add_positions), sequence.tolist()

    