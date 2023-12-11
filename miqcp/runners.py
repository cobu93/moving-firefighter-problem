import gurobipy as gp
from gurobipy import GRB

import numpy as np

class MIQCP:
    
    def __init__(self):
        self.B = 0
        self.D = 0

    def __setup_problem__(self, i_tree, root, t_propagation, m):
        tree, height = i_tree.to_directed(root)
        leaves = np.argwhere(tree.edges.sum(axis=-1) == 0).flatten()
        max_path_len = np.argwhere(tree.edges.sum(axis=-1) == 0).shape[0]
        assert max_path_len >= 1, "There aren't leaves in the tree"
        n_nodes = tree.nodes.shape[0]

        #def mfp_constraints(D, B, n, graph, time, firefighters = 1, return_matrices=False):
        ############# Create the graph ###############
        # Adjacency matrix
        adj = tree.edges
        # Distances matrix
        dists = np.zeros((n_nodes + 1, n_nodes + 1))
        for i, node_i in enumerate(tree.nodes_positions):
            for j, node_j in enumerate(tree.nodes_positions):
                dists[i, j] = np.linalg.norm(node_i - node_j)
        # Rounds definition
        B = int(height)
        D = int(leaves.shape[0])
        n = int(n_nodes + 1)

        self.B = B
        self.D = D
        
        ################ Create the model ##############
        # ------------------------------INPUT --------------------------------
        I = [root]

        b0 = []
        d0 = []
        for i in range(n):
            if i in I:
                b0.append(1)  # These are the vertices burned at time j=0
            else:
                b0.append(0)
        for i in range(n):
            if i == n-1:
                d0.append(1)  # The bulldozer begins at vertex n-1
            else:
                d0.append(0)

        # ---------------------------- VARIABLES -----------------------------------------------------------

        b = []
        for i in range(n):
            temp = []
            for j in range(B):
                temp.append(0)
            b.append(temp)
        for i in range(n):
            for j in range(B):
                b[i][j] = m.addVar(vtype=GRB.BINARY, name="b,%s" % str(i + 1) + "," + str(j + 1))

        d = []
        for i in range(n):
            temp = []
            for j in range(B):
                temp.append(0)
            d.append(temp)
        for i in range(n):
            for j in range(B):
                d[i][j] = m.addVar(vtype=GRB.BINARY, name="d,%s" % str(i + 1) + "," + str(j + 1))

        d_prime = []
        for j in range(B):
            temp_1 = []
            for i in range(n):
                temp_2 = []
                for l in range(D):
                    temp_2.append(0)
                temp_1.append(temp_2)
            d_prime.append(temp_1)
        for j in range(B):
            for i in range(n):
                for l in range(D):
                    d_prime[j][i][l] = m.addVar(vtype=GRB.BINARY,
                                                name="d_prime,%s" % str(j + 1) + "," + str(i + 1) + "," + str(l + 1))

        p = []
        for j in range(B):
            temp_1 = []
            for i in range(n):
                temp_2 = []
                for l in range(D):
                    temp_2.append(0)
                temp_1.append(temp_2)
            p.append(temp_1)
        for j in range(B):
            for i in range(n):
                for l in range(D):
                    p[j][i][l] = m.addVar(vtype=GRB.BINARY,
                                        name="p,%s" % str(j + 1) + "," + str(i + 1) + "," + str(l + 1))

        t = []
        for j in range(B):
            t.append(0)
        for j in range(B):
            t[j] = m.addVar(vtype=GRB.CONTINUOUS, name="t,%s" % str(j + 1))

        y = []
        for j in range(B):
            y.append(0)
        for j in range(B):
            y[j] = m.addVar(vtype=GRB.BINARY, name="y,%s" % str(j + 1))

        # ---------------------------- CONSTRAINTS ---------------------------------------------------------

        for i in range(n):  # ---------------------------( 2 )
            for j in range(B):
                if j == 0:
                    m.addConstr(b[i][j] >= b0[i])
                else:
                    m.addConstr(b[i][j] >= b[i][j - 1])

        for i in range(n):  # -------------------------( 3 )
            for j in range(B):
                if j == 0:
                    m.addConstr(d[i][j] >= d0[i])
                else:
                    m.addConstr(d[i][j] >= d[i][j - 1])

        for i in range(n):  # ---------------------------( 4 )
            for j in range(B):
                m.addConstr(b[i][j] + d[i][j] <= 1)

        for i in range(n-1):  # ---------------------------( 5 )
            for j in range(B):
                for k in range(n-1):
                    if j == 0:
                        m.addConstr(b[i][j] + d[i][j] >= b0[k] * adj[k, i])
                    else:
                        m.addConstr(b[i][j] + d[i][j] >= b[k][j - 1] * adj[k, i])
                # k == n such as d[n][0] = 1, remains defended for every t.

        for i in range(n):  # ---------------------------( 6 y 7 )
            d0[i] = 0
        d0[n-1] = 1

        for i in range(n):  # ---------------------------( 8 )
            for j in range(B):
                if j == 0:
                    m.addConstr(d_prime[j][i][0] >= d0[i])
                else:
                    m.addConstr(d_prime[j][i][0] >= d[i][j - 1])

        for i in range(n):  # ---------------------------( 9 )
            for j in range(B):
                m.addConstr(d_prime[j][i][D - 1] == d[i][j])

        for i in range(n):  # ---------------------------( 10 )
            for j in range(B):
                for k in range(1, D):
                    m.addConstr(d_prime[j][i][k] >= d_prime[j][i][k - 1])

        for j in range(B):  # ---------------------------( 11 y 12 )
            for i in range(n):
                for k in range(D):
                    if k == 0:
                        if j == 0:
                            m.addConstr(p[j][i][k] >= d_prime[j][i][k] - d0[i])
                        else:
                            m.addConstr(p[j][i][k] >= d_prime[j][i][k] - d[i][j - 1])
                    else:
                        m.addConstr(p[j][i][k] >= d_prime[j][i][k] - d_prime[j][i][k - 1])

        for j in range(B):  # ---------------------------( 13 )
            for k in range(D):
                sum_ = 0
                for i in range(n):
                    sum_ = sum_ + p[j][i][k]
                m.addConstr(sum_ == 1)

        for j in range(B):  # ---------------------------( 14, 15 y 16)
            for k in range(D):
                sum_ = 0
                for i in range(n):
                    if k == 0:
                        if j == 0:
                            sum_ = sum_ + d_prime[j][i][k] - d0[i]
                        else:
                            sum_ = sum_ + d_prime[j][i][k] - d[i][j - 1]
                    else:
                        sum_ = sum_ + d_prime[j][i][k] - d_prime[j][i][k - 1]
                for i in range(n):
                    if k == 0:
                        if j == 0:
                            m.addConstr(p[j][i][k] >= d0[i] * (1 - sum_))
                        else:
                            m.addConstr(p[j][i][k] >= p[j - 1][i][D - 1] * (1 - sum_))
                    else:
                        m.addConstr(p[j][i][k] >= p[j][i][k - 1] * (1 - sum_))

        for j in range(B):  # ---------------------------( 17 y 18 )
            sum_1 = 0
            for l in range(n):
                sum_1_a = 0
                for i in range(n):
                    sum_1_a = sum_1_a + p[j][i][0] * dists[i][l]
                if j == 0:
                    sum_1 = sum_1 + sum_1_a * d0[l]
                else:
                    sum_1 = sum_1 + sum_1_a * p[j - 1][l][D - 1]
            sum_2 = 0
            for k in range(1, D):
                sum_2_a = 0
                for l in range(n):
                    sum_2_b = 0
                    for i in range(n):
                        sum_2_b = sum_2_b + p[j][i][k] * dists[i][l]
                    sum_2_a = sum_2_a + sum_2_b * p[j][l][k - 1]
                sum_2 = sum_2 + sum_2_a
            sum_3 = sum_1 + sum_2
            if j == 0:
                m.addConstr(t[j] == sum_3 + 0)
            else:
                m.addConstr(t[j] == sum_3 + t[j - 1])

        for j in range(B):  # ---------------------------( 19 )
            sum_ = 0
            for i in range(n):
                for k in range(D):
                    if k == 0:
                        if j == 0:
                            m.addConstr(y[j] >= d_prime[j][i][k] - d0[i])
                            sum_ = sum_ + d_prime[j][i][k] - d0[i]
                        else:
                            m.addConstr(y[j] >= d_prime[j][i][k] - d[i][j - 1])
                            sum_ = sum_ + d_prime[j][i][k] - d[i][j - 1]
                    else:
                        m.addConstr(y[j] >= d_prime[j][i][k] - d_prime[j][i][k - 1])
                        sum_ = sum_ + d_prime[j][i][k] - d_prime[j][i][k - 1]
            m.addConstr(y[j] <= sum_)


        for j in range(B):  # ---------------------------( 20 y 21 )
            # m.addConstr(t[j] >= j  * y[j])
            # m.addConstr(t[j] <= (j+1) + (1 - y[j]) * M )
            m.addConstr(t[j] <= j + 1)

        # ---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        b_transpose = np.array(b).T.tolist()
        m.setObjective(sum(b_transpose[B - 1]), GRB.MINIMIZE)  # -----------------------------------------------(1)
        # ---------------------------- OPTIMIZATION -------------------------------------------------------

    def run(self, tree, root, initial_ff_position, t_propagation):
        tree.add_firefighter_position(initial_ff_position)
        
        optimal = 0

        with gp.Env() as env, gp.Model(env=env) as m:
            m.setParam("MIPGap", 0)
            m.setParam("NodefileStart", 25)
            #m.setParam("Presolve", 2)
            #m.setParam("PreSparsify", 1)
            
            
            self.__setup_problem__(tree, root, t_propagation, m)
            m.optimize()
            
            optimal = len(tree.nodes) - m.ObjVal
            
            optimal_path = [-1]
            last_defended = []

            for v in m.getVars():
                v_parts = v.varName.split(",")
                if v_parts[0] == "p":
                    last_defended.append(v.x)

            last_defended = np.array(last_defended).reshape(self.B, tree.edges.shape[0] + 1, self.D).transpose(0, 2, 1)

            for br in last_defended:
                for dr in br:
                    last_defended = int(np.argwhere(np.round(dr) == 1)[0, 0])
                    if last_defended != optimal_path[-1]:
                        optimal_path.append(last_defended)

        return optimal, optimal_path
