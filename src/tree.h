#include "tdefinitions.h"

typedef struct {    
    nsize_t n_nodes;
    nsize_t height;
    nsize_t n_leaves;
    nsize_t* nodes;
    float* nodes_x;
    float* nodes_y;
    float* nodes_z;
    float** edges;
} TREE;


nsize_t compute_subtree(TREE* tree, nsize_t node, bool* subtree){
    nsize_t cardinality = 1;
    nsize_t i;
    
    subtree[node] = 1;
    
    for(i = 0; i < tree->n_nodes; i++){
        // If there is an edge node->i
        if(tree->edges[node][i] != 0){
            cardinality += compute_subtree(tree, i, subtree);
        }
    }

    return cardinality;
}

nsize_t compute_parents(TREE* tree, nsize_t node, bool* parents){
    nsize_t cardinality = 1;
    nsize_t i;
    
    parents[node] = 1;
    
    for(i = 0; i < tree->n_nodes; i++){
        // If there is an edge node->i
        if(tree->edges[i][node] != 0){
            cardinality += compute_parents(tree, i, parents);
        }
    }

    return cardinality;
}