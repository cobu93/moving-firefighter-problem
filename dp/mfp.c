#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h> 
#include <unistd.h>

#include "tree.h"
#include "memory.h"

typedef unsigned int uint;
typedef long long int lli;

_Atomic int max_threads_allowed = 10;
_Atomic lli threads_cnt = 0;
pthread_mutex_t mutex_n_threads;

typedef struct {
    TREE* tree; 
    uint height; 
    uint* forest; 
    uint n_forest; 
    uint root; 
    uint firefighter_position; 
    float current_time; 
    float t_propagation;
    MEMORY* memory; 
    uint* hops_memory;
    uint** subtree_memory;
    uint recursion_limit;
    int* opt_path;
    uint* opt_value;
} OPT_PARAMS;


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


void* mfp_dp_opt(void* i_params){

    OPT_PARAMS* p = (OPT_PARAMS*) i_params;

    if(p->recursion_limit <= 0){
        *(p->opt_value) = 0;
        return NULL;
    }
    
    if(p->n_forest == 0){
        *(p->opt_value) = 0;
        return NULL;
    }

    if(p->current_time > p->height * p->t_propagation){
        *(p->opt_value) = 0;
        return NULL;
    }    

    uint n_nodes = (uint) p->tree->n_nodes;

    int memory_value;
    memory_value = search_in_memory(
                        p->memory, 
                        p->forest, 
                        n_nodes, 
                        p->firefighter_position, 
                        p->current_time
                    );
    
    if(memory_value >= 0){
        *(p->opt_value) = memory_value;
        return NULL;
    }

    uint* feasible = (uint*) malloc(n_nodes * sizeof(uint));    
    float* distances = (float*) malloc(n_nodes * sizeof(float));    
    
    get_feasible(
            p->tree, 
            p->forest, 
            p->root, 
            p->firefighter_position, 
            p->current_time, 
            p->t_propagation, 
            feasible, 
            distances, 
            p->hops_memory
        );    
    
    uint i, j;
    int n_subtree; 
    int val, max_val = 0;
    uint* subtree;

    pthread_t* threads = (pthread_t*) malloc(n_nodes * sizeof(pthread_t));
    int* cardinalities = (int*) malloc(n_nodes * sizeof(int));
    int* opt_values = (int*) malloc(n_nodes * sizeof(int));


    OPT_PARAMS** sub_params = (OPT_PARAMS**) malloc(n_nodes * sizeof(OPT_PARAMS*));

    for(i = 0; i < n_nodes; i++){
        if(feasible[i] == 1){

            sub_params[i] = (OPT_PARAMS*) malloc(sizeof(OPT_PARAMS));         
            sub_params[i]->tree = p->tree; 
            sub_params[i]->height = p->height; 
            sub_params[i]->forest = (uint*) calloc(n_nodes, sizeof(uint));
            sub_params[i]->n_forest = 0;
            sub_params[i]->root = p->root; 
            sub_params[i]->firefighter_position = i; 
            sub_params[i]->current_time = p->current_time + distances[i]; 
            sub_params[i]->t_propagation = p->t_propagation;
            sub_params[i]->memory = p->memory; 
            sub_params[i]->hops_memory = p->hops_memory;
            sub_params[i]->subtree_memory = p->subtree_memory;
            sub_params[i]->recursion_limit = p->recursion_limit - 1;
            sub_params[i]->opt_path = (int*) malloc((p->recursion_limit - 1) * sizeof(int));
            sub_params[i]->opt_value = &(opt_values[i]);
            

            cardinalities[i] = 0;
    
            compute_subforest(
                        p->forest, 
                        p->subtree_memory[i], 
                        n_nodes, 
                        sub_params[i]->forest, 
                        &(sub_params[i]->n_forest), 
                        &(cardinalities[i])
                    );   
           
            for(j = 0; j < p->recursion_limit - 1; j++){
                sub_params[i]->opt_path[j] = -1;
            }

            int continue_exec = 0;

            pthread_mutex_lock(&mutex_n_threads);
            continue_exec = max_threads_allowed <= 0? 0: 1;
            pthread_mutex_unlock(&mutex_n_threads);

            while(!continue_exec){
                pthread_mutex_lock(&mutex_n_threads);
                continue_exec = max_threads_allowed <= 0? 0: 1;
                pthread_mutex_unlock(&mutex_n_threads);
                usleep(1000);
            }

            pthread_mutex_lock(&mutex_n_threads);
            max_threads_allowed--;
            threads_cnt++;
            printf("Reserved\n");
            pthread_mutex_unlock(&mutex_n_threads);


            if (pthread_create(
                    &threads[i], 
                    NULL, 
                    mfp_dp_opt, 
                    sub_params[i]) != 0){
                printf("Error: Cannot create thread %d\n", i);
                break;
            }
        }
    }

    for(i = 0; i < n_nodes; i++){
        if(feasible[i] == 1){

            if(pthread_join(threads[i], NULL) != 0){
                printf("Error: Cannot join thread %d\n", i);
            }

            pthread_mutex_lock(&mutex_n_threads);
            max_threads_allowed++;            
            printf("Freed\n");
            pthread_mutex_unlock(&mutex_n_threads);

            free(sub_params[i]->forest);

            val = cardinalities[i] + *(sub_params[i]->opt_value);

            if(val > max_val){
                max_val = val;
                p->opt_path[0] = i;

                for(j = 0; j < p->recursion_limit - 1; j++){
                    p->opt_path[j + 1] = sub_params[i]->opt_path[j];
                }

            }

            free(sub_params[i]->opt_path);
            free(sub_params[i]);
        }
    }
    
    free(opt_values);
    free(sub_params);
    free(cardinalities);
    free(threads);
    free(feasible);
    free(distances);

    int result = add_to_memory(
                    p->memory, 
                    max_val, 
                    p->forest, 
                    p->n_forest, 
                    p->firefighter_position, 
                    p->current_time
                );

    if(result < 0){
        delete_memory(&(p->memory));
        free(p->forest);
        free(p->hops_memory);
    }

    *(p->opt_value) = max_val;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void mfp_dp_solver(TREE tree, int root, int firefighter_position, float t_propagation, int* optimal, int* opt_path){

    uint n_forest, height, n_leaves;
    uint* forest;
    uint* hops_memory;
    uint** subtree_memory;
    uint i, j = 0;

    uint opt = 0;

    //height = tree_to_directed(&tree, root);
    //n_leaves = compute_n_leaves(&tree);

    height = (uint) tree.height;
    n_leaves = (uint) tree.n_leaves;
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
   
    for(i = 0; i < n_leaves; i++){
        opt_path[i] = -1;
    }
    
    MEMORY* memory;
    init_memory(&memory);

    if (pthread_mutex_init(&mutex_n_threads, NULL) != 0) { 
        printf("Threads mutex init has failed\n"); 
    } 

    OPT_PARAMS params;

    params.tree = &tree; 
    params.height = height; 
    params.forest = forest; 
    params.n_forest = n_forest; 
    params.root = (uint) root; 
    params.firefighter_position = (uint) firefighter_position; 
    params.current_time = 0; 
    params.t_propagation = t_propagation;
    params.memory = memory; 
    params.hops_memory = hops_memory;
    params.subtree_memory = subtree_memory;
    params.recursion_limit = n_leaves;
    params.opt_path = opt_path;
    params.opt_value = &opt;

    mfp_dp_opt(&params);

    if(opt >= tree.n_nodes){
        opt = -1;
    }

    printf("Max threads allowed %d\n", max_threads_allowed);
    printf("Run threads %d\n", threads_cnt);

    pthread_mutex_destroy(&mutex_n_threads);

    
    delete_memory(&memory);
    free(forest);
    free(hops_memory);
    free(subtree_memory);

    (*optimal) = opt;
    
}
