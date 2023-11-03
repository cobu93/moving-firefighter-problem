#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "tree.h"
#include "memory.h"

typedef unsigned int uint;


float tau(
    TREE* tree,
    uint origin,
    uint destiny
){
    float distance;
    distance = pow(tree->nodes_x[origin] - tree->nodes_x[destiny], 2);
    distance += pow(tree->nodes_y[origin] - tree->nodes_y[destiny], 2);
    distance += pow(tree->nodes_z[origin] - tree->nodes_z[destiny], 2);

    return sqrt(distance);
}

uint hops(
    TREE* tree, 
    uint root, 
    uint node,
    uint* hops_memory
    ){

    if(node == root){
        return 0;
    }

    if(hops_memory[node] != 0){
        return hops_memory[node];
    }

    uint i, j;
    uint c_node = node;
    uint hops = 0;


    for(i = 0; i < tree->n_nodes; i++){

        if(c_node == root){
            break;
        }

        for(j = 0; j < tree->n_nodes; j++){
            if(tree->edges[j][c_node] != 0){
                c_node = j;
                hops++;
                break;
            }
        }

    }

    hops_memory[node] = hops;

    return hops;
}

void get_feasible(
    TREE* tree,
    uint* forest, 
    uint root,
    uint firefighter_position, 
    float current_time,
    float t_propagation,
    uint* feasible,
    float* distances,
    uint* hops_memory
    ){
    
    uint i, n_hops;
    float distance;

    for(i = 0; i < tree->n_nodes; i++){
        distance = 0;
        n_hops = 0;
        if(forest[i] == 1){
            distance = tau(tree, firefighter_position, i);
            n_hops = hops_memory[i];
        }
        
        distances[i] = distance;
        feasible[i] = 0;
        
        if(current_time + distance < n_hops * t_propagation){
            feasible[i] = 1;
        }
    }
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

void compute_subforest(uint* forest, uint* subtree, uint n_nodes, uint* subforest, uint* n_subforest, uint* cardinality){

    uint i;

    (*n_subforest) = 0;
    (*cardinality) = 0;

    for(i = 0; i < n_nodes; i++){
        if(forest[i] == 1 && subtree[i] == 0){
            subforest[i] = 1;
            (*n_subforest)++;
        }

        if(subtree[i] == 1 && forest[i] == 1){
            (*cardinality)++;
        }
    }
}


int mfp_dp_opt(
    TREE* tree, 
    uint height, 
    uint* forest, 
    uint n_forest, 
    uint root, 
    uint firefighter_position, 
    float current_time, 
    float t_propagation,
    MEMORY* memory, 
    uint* hops_memory,
    uint** subtree_memory
    ){
    
    if(n_forest == 0){
        return 0;
    }

    if(current_time > height * t_propagation){
        return 0;
    }    

    uint n_nodes = (uint) tree->n_nodes;

    int memory_value;
    memory_value = search_in_memory(memory, forest, n_nodes, firefighter_position, current_time);
    
    if(memory_value >= 0){
        return memory_value;
    }

    uint* feasible = (uint*) malloc(n_nodes * sizeof(uint));    
    float* distances = (float*) malloc(n_nodes * sizeof(float));    
    get_feasible(tree, forest, root, firefighter_position, current_time, t_propagation, feasible, distances, hops_memory);    
    
    uint i;
    int n_subtree, n_subforest, cardinality;
    int val, max_val = 0;
    uint* subtree;
    uint* subforest;

    subforest = (uint*) malloc(n_nodes * sizeof(uint));
    
    for(i = 0; i < n_nodes; i++){
        if(feasible[i] == 1){
            memset(subforest, 0, n_nodes * sizeof(uint));
            n_subforest = 0;
            cardinality = 0;

            compute_subforest(forest, subtree_memory[i], n_nodes, subforest, &n_subforest, &cardinality);   

            val = cardinality + mfp_dp_opt(
                                    tree, 
                                    height, 
                                    subforest, 
                                    n_subforest, 
                                    root, 
                                    i, 
                                    current_time + distances[i],
                                    t_propagation,
                                    memory,
                                    hops_memory,
                                    subtree_memory
                                );

            if(val > max_val){
                max_val = val;
            }
        }
    }

    free(subforest);
    free(feasible);
    free(distances);

    int result = add_to_memory(memory, max_val, forest, n_forest, firefighter_position, current_time);

    if(result < 0){
        delete_memory(&memory);
        free(forest);
        free(hops_memory);
    }

    return max_val;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void mfp_dp_solver(TREE tree, int root, int firefighter_position, float t_propagation, int* optimal){

    uint height, n_forest;
    uint* forest;
    uint* hops_memory;
    uint** subtree_memory;
    uint i, j = 0;

    uint opt = 0;

    height = tree_to_directed(&tree, root);
    forest = (uint*)malloc(tree.n_nodes * sizeof(uint));
    hops_memory = (uint*) calloc(tree.n_nodes, sizeof(uint));
    subtree_memory = (uint**) calloc(tree.n_nodes, sizeof(uint*));

    for(i = 0; i < tree.n_nodes; i++){
        subtree_memory[i] = (int*) calloc(tree.n_nodes, sizeof(uint));
        compute_subtree(&tree, i, subtree_memory[i]);
        hops_memory[i] = hops(&tree, root, i, hops_memory);        
    }

    n_forest = tree.n_nodes - 1;
    
    for(i = 0; i < tree.n_nodes; i++){
        if(i == root){
            forest[i] = 0;
        }
        else{
            forest[i] = 1;
        }
    }
   
    
    MEMORY* memory;
    init_memory(&memory);

    opt = mfp_dp_opt(&tree, height, forest, n_forest, (uint) root, (uint) firefighter_position, 0, t_propagation, memory, hops_memory, subtree_memory);

    if(opt >= tree.n_nodes){
        opt = -1;
    }

    delete_memory(&memory);
    free(forest);
    free(hops_memory);
    free(subtree_memory);

    (*optimal) = opt;
}