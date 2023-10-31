#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {    
    int n_nodes;
    int* nodes;
    float* nodes_x;
    float* nodes_y;
    float* nodes_z;

    int n_edges;
    int* edges_o;
    int* edges_d;
} TREE;

float tau(
    TREE tree,
    int origin,
    int destiny
){
    float distance;
    distance = pow(tree.nodes_x[origin] - tree.nodes_x[destiny], 2);
    distance += pow(tree.nodes_y[origin] - tree.nodes_y[destiny], 2);
    distance += pow(tree.nodes_z[origin] - tree.nodes_z[destiny], 2);

    return sqrt(distance);
}

int hops(
    TREE tree, 
    int root, 
    int node){

    

    return 100;

}

int* get_feasible(
    TREE tree,
    int* forest, 
    int n_forest, 
    int root,
    int firefighter_position, 
    int timestep
    ){
    
    int i, n_hops;
    float distance;

    for(i = 0; i < n_forest; i++){
        distance = tau(tree, firefighter_position, forest[i]);
        n_hops = hops(tree, root, forest[i]);
        if(timestep + distance < n_hops){
            printf("From %d to %d are %f [%d]\n", firefighter_position, forest[i], distance, n_hops);
        }
    }
}


float mfp_dp_opt(TREE tree, int* forest, int n_forest, int root, int firefighter_position, int timestep){
    get_feasible(tree, forest, n_forest, root, firefighter_position, timestep);
    return 2;
}


float mfp_dp_solver(TREE tree, int root, int firefighter_position){

    int timestep = 0;
    int n_forest = tree.n_nodes - 1;
    int* forest = (int*)malloc(n_forest * sizeof(int));
    int i, j = 0;

    float opt = 0;

    for(i = 0; i < tree.n_nodes; i++){
        if(tree.nodes[i] != root){
            forest[j] = tree.nodes[i];
            j++;
        }
    }

    opt = mfp_dp_opt(tree, forest, n_forest, root, firefighter_position, 0);
    free(forest);

    return opt;
}