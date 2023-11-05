import gurobipy as gp
from gurobipy import GRB

import numpy as np

class ILP:
    
    def __init__(self):
        pass

    def __setup_problem__(self, i_tree, root, t_propagation, m):
        tree, _ = i_tree.to_directed(root)
        leaves = np.argwhere(tree.edges.sum(axis=-1) == 0).flatten()
        max_path_len = np.argwhere(tree.edges.sum(axis=-1) == 0).shape[0]
        assert max_path_len >= 1, "There aren't leaves in the tree"
        n_nodes = tree.nodes.shape[0]
        
        # Shape: (nodes + 1, nodes, max_path_len)
        # The extra value in the first dimension corresponds to a
        X = m.addMVar(shape=(n_nodes, n_nodes, max_path_len - 1), vtype=GRB.BINARY, name="X")
        A = m.addMVar(shape=n_nodes, vtype=GRB.BINARY, name="A")
        
        if max_path_len > 1:
            m.addConstr(A + X.sum(axis=(0, -1)) <= 1, f"(13)")
        else:
            m.addConstr(A <= 1, f"(13)")

        if max_path_len > 1:
            m.addConstr(X[:, :, 0].sum(axis=(0, 1)) <= A.sum(), "(14.0)")
            if max_path_len > 2:
                m.addConstr(X[:, :, 1:].sum(axis=(0, 1)) <= X[:, :, :-1].sum(axis=(0, 1)), "(14)")
            

        
        paths_to_root = []
        for i in range(n_nodes):
            paths_to_root.append(tree.get_path_to_root(i))

        for i in leaves:
            if max_path_len > 1:
                m.addConstr(
                    (A[paths_to_root[i]] + X[:, paths_to_root[i], :].sum(axis=(0, -1))).sum() <= 1, 
                    f"(15.{i})"
                )
            else:
                m.addConstr(
                    A[paths_to_root[i]].sum() <= 1, 
                    f"(15.{i})"
                )


        # Compute distances from initial to each node
        i_distances = np.linalg.norm(tree.nodes_positions[:-1] - tree.nodes_positions[-1], axis=-1)
        distances = np.zeros((n_nodes, n_nodes))

        for i in range (n_nodes):
            for j in range (n_nodes):
                distances[i, j] = np.linalg.norm(tree.nodes_positions[i] - tree.nodes_positions[j])

        C = (n_nodes + 1) * np.max(distances) + np.max(i_distances)
        hops = np.array([paths_to_root[j].shape[0] - 1 for j in range(n_nodes)]) * t_propagation
        
        
        ###################################################################
        m.addConstr((i_distances * A).sum() <= (hops * A).sum() + C * (1 - A.sum()), f"(16.base)")
        
        for i in range(max_path_len - 1):

            feasible_cstr = gp.LinExpr()
            feasible_cstr += (i_distances * A).sum()

            for j in range(n_nodes):
                for k in range(n_nodes):
                    for l in range(i + 1): 
                        feasible_cstr += distances[j, k] * X[j, k, l]

            comp_feasible_cstr = gp.LinExpr()

            for j in range(n_nodes):
                comp_feasible_cstr += (hops * X[j, :, i]).sum()
            
            comp_feasible_cstr += C * (1 - X[:, :, i].sum())#

            m.addConstr(feasible_cstr <= comp_feasible_cstr, f"(16.{i})")
        ###################################################################3

        
        for i in range(n_nodes):
            path_cstr = gp.LinExpr()
            path_cstr += A[i]

            if max_path_len > 1:
                for j in range(n_nodes):
                    if (i != j):
                        path_cstr += X[j, :, 0].sum()
                
            m.addConstr(path_cstr <= 1, f"(17.{i})")
        
        
        for i in range(max_path_len - 2):
            for j in range(n_nodes):
                path_cstr = gp.LinExpr()
                path_cstr += X[:, j, i]

                for k in range(n_nodes):
                    if (k != j):
                        path_cstr += X[k, :, i + 1].sum()
                
                m.addConstr(path_cstr <= 1, f"(18.{i})")
        
        m.addConstr(A[root] == 0, f"(No root available ini)")
        m.addConstr(A.sum() <= 1, f"(One initial node)")
        if max_path_len > 1:
            m.addConstr(X[:, root].sum() == 0, f"(No root available)")
            
        subtrees_cardinalities = np.array([tree.get_subtree_nodes(i).shape[0] for i in range(n_nodes)])
        
        if max_path_len > 1:
            m.setObjective(
                (A * subtrees_cardinalities).sum() + (X[:, :, :].sum(axis=(0, -1)) * subtrees_cardinalities).sum(), 
                GRB.MAXIMIZE
                )
        else:
            m.setObjective(
                (A * subtrees_cardinalities).sum(), 
                GRB.MAXIMIZE
                )
        
        
    def run(self, tree, root, initial_ff_position, t_propagation):
        tree.add_firefighter_position(initial_ff_position)
        
        optimal = 0


        with gp.Env() as env, gp.Model(env=env) as m:
            m.setParam("MIPGap", 0)

            self.__setup_problem__(tree, root, t_propagation, m)
            m.optimize()
            
            optimal = m.ObjVal
            result = []
            first = []
            for v in m.getVars():
                if v.varName.startswith("A"):
                    first.append(v.X)
                else:
                    result.append(v.X)

            first = np.array(first)
            result = np.array(result).reshape(tree.nodes.shape[0], tree.nodes.shape[0], -1)

            optimal_path = [-1] + [int(np.argwhere(first == 1)[0, 0])]
            for i in range(result.shape[-1]):
                next_node = np.argwhere(result[:, :, i] == 1)
                if next_node.shape[0] > 0:
                    optimal_path += [int(np.argwhere(result[:, :, i] == 1)[0, 1])]

        return optimal, optimal_path
