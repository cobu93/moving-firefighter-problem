
typedef struct {    
    int n_nodes;
    int height;
    int n_leaves;
    int* nodes;
    float* nodes_x;
    float* nodes_y;
    float* nodes_z;
    float** edges;
} TREE;

uint subtree_to_directed(TREE* tree, uint node, uint* visited){
    uint i, height;
    uint max_height = 0;
    
    visited[node] = 1;

    for(i = 0; i < tree->n_nodes; i++){
        if(visited[i] == 0 && tree->edges[node][i] != 0){
            tree->edges[i][node] = 0;
            height = subtree_to_directed(tree, i, visited);
            if(height > max_height){
                max_height = height;
            }
        }
    }

    return 1 + max_height;
}

uint tree_to_directed(TREE* tree, uint root){
    uint height;
    uint* visited = (uint*) calloc(tree->n_nodes, sizeof(uint));

    height = subtree_to_directed(tree, root, visited);
    free(visited);
    return height;
}

uint compute_n_leaves(TREE* tree){
    uint i, j, leaves = 0;
    float n_edges;

    for(i=0; i < tree->n_nodes; i++){
        n_edges = 0;
        for(j=0; j < tree->n_nodes; j++){
            n_edges += tree->edges[i][j];
        }

        if(n_edges == 0){
            leaves++;
        }
    }

    return leaves;
}

uint compute_subtree(TREE* tree, uint node, uint* subtree){
    uint cardinality = 1;
    uint i;
    
    subtree[node] = 1;
    
    for(i = 0; i < tree->n_nodes; i++){
        // If there is an edge node->i
        if(tree->edges[node][i] != 0){
            cardinality += compute_subtree(tree, i, subtree);
        }
    }

    return cardinality;
}

uint compute_parents(TREE* tree, uint node, uint* parents){
    uint cardinality = 1;
    uint i;
    
    parents[node] = 1;
    
    for(i = 0; i < tree->n_nodes; i++){
        // If there is an edge node->i
        if(tree->edges[i][node] != 0){
            cardinality += compute_parents(tree, i, parents);
        }
    }

    return cardinality;
}