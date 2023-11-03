
typedef struct {    
    int n_nodes;
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