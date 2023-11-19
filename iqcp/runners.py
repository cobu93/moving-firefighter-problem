import gurobipy as gp
from gurobipy import GRB

import numpy as np

class IQCP:
    
    def __init__(self):
        pass

    def __setup_problem__(self, i_tree, root, t_propagation, m):
        tree, _ = i_tree.to_directed(root)
        leaves = np.argwhere(tree.edges.sum(axis=-1) == 0).flatten()
        max_path_len = np.argwhere(tree.edges.sum(axis=-1) == 0).shape[0]
        assert max_path_len >= 1, "There aren't leaves in the tree"
        n_nodes = tree.nodes.shape[0]

        
        # Shapé: (nodes, max_path_len)
        X = m.addMVar(shape=(n_nodes, max_path_len), vtype=GRB.BINARY, name="X")
        m.addConstr(X[root].sum() == 0, f"(5)")
        m.addConstr(X.sum(axis=-1) <= 1, "(6)")
        m.addConstr(X.sum(axis=0) <= 1, "(7)")
        m.addConstr(X.sum(axis=0)[1:] <= X.sum(axis=0)[:-1], "(8)")
        
        paths_to_root = []
        for i in range(n_nodes):
            paths_to_root.append(tree.get_path_to_root(i))

        for i in leaves:
            m.addConstr(X[paths_to_root[i]].sum() <= 1, f"(9.{i})")

        
        # Compute distances from initial to each node
        i_distances = np.linalg.norm(tree.nodes_positions[:-1] - tree.nodes_positions[-1], axis=-1)
        distances = np.zeros((n_nodes, n_nodes))

        for i in range (n_nodes):
            for j in range (n_nodes):
                distances[i, j] = np.linalg.norm(tree.nodes_positions[i] - tree.nodes_positions[j])

        C = (max_path_len - 1) * np.max(distances) + np.max(i_distances)
        hops = np.array([paths_to_root[j].shape[0] - 1 for j in range(n_nodes)]) * t_propagation
                
        for i in range(max_path_len):

            feasible_cstr = gp.LinExpr()
            feasible_cstr += (i_distances * X[:, 0]).sum()

            for j in range(n_nodes):
                for k in range(n_nodes):
                    for l in range(i):
                        feasible_cstr += distances[j, k] * X[j, l] * X[k, l + 1]

            m.addConstr(feasible_cstr <= (hops * X[:, i]).sum() + C * (1 - X[:, i].sum()), f"(10.{i})")
        
        
        subtrees_cardinalities = np.array([tree.get_subtree_nodes(i).shape[0] for i in range(n_nodes)])
        m.setObjective((X.sum(axis=-1) * subtrees_cardinalities).sum(), GRB.MAXIMIZE)
        
        
    def run(self, tree, root, initial_ff_position, t_propagation):
        tree.add_firefighter_position(initial_ff_position)
        
        optimal = 0


        with gp.Env() as env, gp.Model(env=env) as m:
            m.setParam("MIPGap", 0)
            
            self.__setup_problem__(tree, root, t_propagation, m)
            m.optimize()
            
            optimal = m.ObjVal
            
            optimal_path = []
            for v in m.getVars():
                optimal_path.append(v.X)

            optimal_path = np.array(optimal_path).reshape(tree.nodes.shape[0], -1)
            optimal_path = [-1] + np.argwhere(optimal_path.T == 1)[:, 1].tolist()
            
        return optimal, optimal_path
