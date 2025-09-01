import gurobipy as gp
from gurobipy import GRB

import numpy as np
import signal
import os
import joblib

class EILP:
    
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
        X = m.addMVar(shape=(n_nodes + 1, n_nodes), vtype=GRB.BINARY, name="X")
        X_path = m.addMVar(shape=n_nodes, vtype=GRB.BINARY, name="X_path")
        t_v = m.addMVar(shape=n_nodes + 1, vtype=GRB.CONTINUOUS, name="t_v")

        m.addConstr(X[0].sum() == 1, f"(24)")
        m.addConstr(X[1:].sum(axis=-1) <= X_path, f"(25)")
        m.addConstr(X.sum(axis=0) == X_path, f"(26)")

        diag_cstr = gp.LinExpr()
        for i in range(n_nodes):
            diag_cstr += X[1 + i, i]
            
        m.addConstr(diag_cstr <= 0, f"(27)")

        for n_leave, leave in enumerate(leaves):
            ancestors = tree.get_path_to_root(leave)
            m.addConstr(X_path[ancestors].sum() <= 1, f"(28.{n_leave})")

        m.addConstr(t_v[0] == 0, f"(29)")
        m.addConstr(t_v >= 0, f"(30)")

        # Compute distances from initial to each node
        i_distances = np.linalg.norm(tree.nodes_positions[:-1] - tree.nodes_positions[-1], axis=-1)
        distances = np.zeros((n_nodes, n_nodes))

        for i in range (n_nodes):
            for j in range (n_nodes):
                distances[i, j] = np.linalg.norm(tree.nodes_positions[i] - tree.nodes_positions[j])

        paths_to_root = []
        for i in range(n_nodes):
            paths_to_root.append(tree.get_path_to_root(i))

        C = (max_path_len - 1) * np.max(distances) + np.max(i_distances)
        hops = np.array([paths_to_root[j].shape[0] - 1 for j in range(n_nodes)]) * t_propagation


        for v_i, v in enumerate(range(1, n_nodes + 1)):
            m.addConstr(t_v[v] >= t_v[0] + i_distances[v - 1] - C * (1 - X[0, v - 1]), f"(31.base.{v_i})")
               
        for u_i, u in enumerate(range(1, n_nodes + 1)):
            for v_i, v in enumerate(range(1, n_nodes + 1)):
                 m.addConstr(t_v[v] >= t_v[u] + distances[u - 1, v - 1] - C * (1 - X[u, v - 1]), f"(31.{u_i}.{v_i})")

        for v_i, v in enumerate(range(1, n_nodes + 1)):
            m.addConstr(t_v[v] <= hops[v - 1] + C * (1 - X_path[v - 1]), f"(32.{v_i})")
            
        ########################
        subtrees_cardinalities = np.array([tree.get_subtree_nodes(i).shape[0] for i in range(n_nodes)])
        
        m.setObjective(
            (subtrees_cardinalities * X_path).sum(), 
            GRB.MAXIMIZE
            )
        
    def get_extra_info(self, checkpoint_prefix):
        
        if os.path.exists(f"{checkpoint_prefix}.jl"):
            return joblib.load(f"{checkpoint_prefix}.jl")
        
        return {}       
        
    def run(self, tree, root, initial_ff_position, t_propagation, checkpoint_prefix):
        tree.add_firefighter_position(initial_ff_position)
        optimal = 0


        with gp.Env() as env, gp.Model(env=env) as m:
            m.setParam("MIPGap", 0)
            m.setParam("NodefileStart", 25)
            # m.setParam("Presolve", 2)
            # m.setParam("PreSparsify", 1)
            

            self.__setup_problem__(tree, root, t_propagation, m)

            def handle_sigterm(signum, frame):
                m.terminate()

            signal.signal(signal.SIGTERM, handle_sigterm)

            m.optimize()
            m.write(f"{checkpoint_prefix}.lp")
            m.write(f"{checkpoint_prefix}.sol")

            joblib.dump({
                "primal_bound": m.ObjVal,
                "dual_bound": m.ObjBound,
                "gap": m.MIPGap
            }, f"{checkpoint_prefix}.jl")
            
            optimal = m.ObjVal

            result = []
            first = []
            for v in m.getVars():
                if v.varName.startswith("X[0"):
                    first.append(v.X)
                elif v.varName.startswith("X["):
                    result.append(v.X)

            first = np.array(first)
            result = np.array(result).reshape(tree.nodes.shape[0], tree.nodes.shape[0], -1)

            optimal_path = [-1] 
            
            next_node = np.argwhere(first == 1)
            if next_node.shape[0] > 0:
                optimal_path += [int(next_node[0, 0])]            

            for i in range(result.shape[-1]):
                next_node = np.argwhere(result[:, :, i] == 1)
                if next_node.shape[0] > 0:
                    optimal_path += [int(next_node[0, 1])]


        return optimal, optimal_path