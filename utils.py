from ctypes import *
import numpy as np

class Tree:
    def __init__(self, nodes, edges, nodes_positions=None, edges_as_matrix=False):
        self.nodes = nodes
        if(nodes_positions is not None):
            self.nodes_positions = np.array([nodes_positions[i] for i in self.nodes])
        self.height = 0
        self.__is_directed__ = False
        self.root = None
        self.height = 0
        
        if edges_as_matrix:
            self.edges = edges
        else:
            self.edges = np.zeros((len(nodes), len(nodes)))

            for o, d in edges:
                self.edges[o, d] = 1
                self.edges[d, o] = 1
    
    @property
    def is_directed(self):
        return self.__is_directed__

    def add_firefighter_position(self, pos):
        if (self.nodes_positions.shape[0] == self.nodes.shape[0]):
            self.nodes_positions = np.concatenate((self.nodes_positions, [pos]), axis=0)
        else:
            self.nodes_positions[self.nodes.shape[0]] = pos

    def __subtree_to_directed__(self, tree, node, visited):
        height = 0
        max_height = 0
        visited[node] = True

        for o in tree.nodes:
            if not visited[o] and tree.edges[node][o] != 0:
                tree.edges[o][node] = 0
                height = self.__subtree_to_directed__(tree, o, visited)
                if(height > max_height):
                    max_height = height

        return 1 + max_height
                
    def to_directed(self, root):
        
        d_tree = Tree(
                    np.copy(self.nodes), 
                    np.copy(self.edges), 
                    None, 
                    edges_as_matrix=True
                )
        
        d_tree.nodes_positions = np.copy(self.nodes_positions)

        visited = [False] * d_tree.nodes.shape[0]
        height = self.__subtree_to_directed__(d_tree, root, visited)
        d_tree.__is_directed__ = True
        d_tree.root = root
        d_tree.height = height
    
        return d_tree, height
    
    def get_path_to_root(self, node):
        assert self.__is_directed__, "The tree must be converted to directed"

        ancestors = [node]
        c_node = node

        for _ in range(self.nodes.shape[0]):
            if c_node != self.root:
                c_node = np.argwhere(self.edges.T[c_node] == 1)[0, 0]
            else:
                break
            
            ancestors.append(c_node)
            
            
        return np.array(ancestors)

    def get_subtree_nodes(self, node):
        assert self.__is_directed__, "The tree must be converted to directed"

        nodes = [node]
        c_node = node

        nodes_idx = 0

        for _ in range(self.nodes.shape[0]):
            if(nodes_idx >= len(nodes)):
                break

            c_node = nodes[nodes_idx]
            next_nodes = np.argwhere(self.edges[c_node] == 1).flatten()
            nodes += next_nodes.tolist()
            nodes_idx += 1

        return np.array(nodes)
        
class TREE(Structure):
    _fields_ = [
        ("n_nodes", c_uint8),
        ("height", c_uint8),
        ("n_leaves", c_uint8),
        ("nodes", POINTER(c_uint8)),
        ("nodes_x", POINTER(c_float)),
        ("nodes_y", POINTER(c_float)),
        ("nodes_z", POINTER(c_float)),
        ("egdes", POINTER(POINTER(c_float)))
    ]


def tree_to_structure(tree):
    assert tree.is_directed, "The tree must be converted to directed"

    n_nodes = tree.nodes.shape[0]
    n_positions = tree.nodes_positions.shape[0]
    n_leaves = (int) (np.argwhere(tree.edges.sum(axis=-1) == 0).flatten().shape[0])

    return TREE(
        n_nodes,
        tree.height,
        n_leaves,
        (c_uint8 * n_nodes)(*tree.nodes.tolist()),
        (c_float * n_positions)(*tree.nodes_positions[:, 0].tolist()),
        (c_float * n_positions)(*tree.nodes_positions[:, 1].tolist()),
        (c_float * n_positions)(*tree.nodes_positions[:, 2].tolist()),

        (POINTER(c_float) * n_nodes)(*[ (c_float * n_nodes)(*r) for r in tree.edges.tolist()])
        )
        

