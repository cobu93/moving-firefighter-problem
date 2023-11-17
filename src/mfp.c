#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h> 
#include <unistd.h>

#include "tdefinitions.h"
#include "tree.h"
#include "memory.h"


// Maximum hreads allowed.
#define MAX_THREADS_ALLOWED 2000000000

/*
*
* Definitions
*
*/
_Atomic long long int __threads_available__ = MAX_THREADS_ALLOWED;
pthread_mutex_t __mutex_n_threads__;
void* mfp_dp_opt(void* i_params);

/*
Structure to feed the recursive dynamic programming function
*/

typedef struct {
    TREE* tree; // Tree
    nsize_t height; // Tree's height
    bool* forest; // Forest for current step
    nsize_t n_forest; // Forest's size
    nsize_t root; // Tree's root
    nsize_t firefighter_position; // Current FF position
    float current_time; // Current time
    float t_propagation; // Propagation time
    MEMORY* memory; // Memory for memoization (If null no memoization considered)
    nsize_t* hops_memory; // Precomputed hops for each node in tree
    bool** subtree_memory; // Precomputed subtree memory for each node
    bool** parents_memory; // Precomputed parents memory for each node
    nsize_t recursion_limit; // Recursion limit (Number of leaves)
    msize_t* opt_path; // Optimal path for recursion step
    nsize_t* opt_value; // Saved nodes for optimal path
} OPT_PARAMS;

typedef struct {
    OPT_PARAMS* old_params;
    nsize_t node;
    float distance;
    msize_t* opt_path; 
    nsize_t* opt_value;
    nsize_t* cardinality;
} SF_OPT_PARAMS;



/*
*
* Utils functions
*
*/
bool reserve_thread(){
    /*
    Description: Attempt to reserve a thread
    Outputs: 1 if success else 0
    */

    bool reserved = 0;

    pthread_mutex_lock(&__mutex_n_threads__);
    if(__threads_available__ > 0){
        __threads_available__--;
        reserved = 1;
        }    
    pthread_mutex_unlock(&__mutex_n_threads__);

    
    return reserved;
}

void free_thread(){
    /*
    Description: Free space for new thread reservation
    Outputs: None
    */

    pthread_mutex_lock(&__mutex_n_threads__);
    __threads_available__++;            
    pthread_mutex_unlock(&__mutex_n_threads__);
}



float tau(
    TREE* tree,
    uint8_t origin,
    uint8_t destiny
){
    /*
    Description: Compute distance between two points
    Outputs: Distance between origin and destiny
    */

    float distance;
    distance = pow(tree->nodes_x[origin] - tree->nodes_x[destiny], 2);
    distance += pow(tree->nodes_y[origin] - tree->nodes_y[destiny], 2);
    distance += pow(tree->nodes_z[origin] - tree->nodes_z[destiny], 2);

    return sqrt(distance);
}

nsize_t hops(
    TREE* tree, 
    nsize_t root, 
    nsize_t node,
    nsize_t* hops_memory
    ){
    
    /*
    Description: Compute the number of hops between root and node
    Outputs: Number of hops, and fill hops memory with that value
    */

    if(node == root){
        return 0;
    }

    if(hops_memory[node] != 0){
        return hops_memory[node];
    }

    nsize_t i, j;
    nsize_t c_node = node;
    nsize_t hops = 0;


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
    bool* forest, 
    nsize_t root,
    nsize_t firefighter_position, 
    float current_time,
    float t_propagation,
    bool* feasible,
    float* distances,
    nsize_t* hops_memory
    ){

    /*
    Description: Compute the feasible nodes given the FF position
    Outputs: Fills feasible array with feasible nodes
    */
    
    nsize_t i, n_hops;
    float distance;

    for(i = 0; i < tree->n_nodes; i++){
        distance = -1;
        n_hops = 0;
        if(forest[i] == 1){
            distance = tau(tree, firefighter_position, i);
            n_hops = hops_memory[i];

            if(current_time + distance < n_hops * t_propagation){
                feasible[i] = 1;
            }
            else{
                feasible[i] = 0;
            }
            
        }
        else{
            feasible[i] = 0;
        }

        distances[i] = distance;
        
    }
}

void compute_subforest(
    bool* forest, 
    bool* subtree, 
    bool* parents, 
    nsize_t n_nodes, 
    bool* subforest, 
    nsize_t* n_subforest, 
    nsize_t* cardinality){

    /*
    Description: Compute the subforest given the subtree and parents
    Outputs: Fills subforest, subforest size (n_subforest) and cardinality of subtree.
    */

    nsize_t i;

    (*n_subforest) = 0;
    (*cardinality) = 0;

    for(i = 0; i < n_nodes; i++){
        subforest[i] = 0;

        if(forest[i] == 1 && subtree[i] == 0 && parents[i] == 0){
            subforest[i] = 1;
            (*n_subforest)++;
        }

        if(subtree[i] == 1 && forest[i] == 1){
            (*cardinality)++;
        }
    }
}

void* sf_mfp_dp_opt(void* params){

    /*
    Description: Recursive function for MFP (Dynamic programming)
    Outputs: Fills optimal value and optimal path.
    */

    SF_OPT_PARAMS* sfp  = (SF_OPT_PARAMS*) params;
    OPT_PARAMS* p = sfp->old_params;

    nsize_t j;
    nsize_t n_nodes = p->tree->n_nodes;

    // Create params
    OPT_PARAMS* sub_params = (OPT_PARAMS*) malloc(sizeof(OPT_PARAMS));
    sub_params->tree = p->tree; 
    sub_params->height = p->height; 
    sub_params->forest = (bool*) malloc(n_nodes * sizeof(bool));
    sub_params->n_forest = 0;
    sub_params->root = p->root; 
    sub_params->firefighter_position = sfp->node; 
    sub_params->current_time = p->current_time + sfp->distance; 
    sub_params->t_propagation = p->t_propagation;
    sub_params->memory = p->memory; 
    sub_params->hops_memory = p->hops_memory;
    sub_params->subtree_memory = p->subtree_memory;
    sub_params->parents_memory = p->parents_memory;
    sub_params->recursion_limit = p->recursion_limit - 1;
    sub_params->opt_path = sfp->opt_path;
    sub_params->opt_value = sfp->opt_value;

    // Initialize optimal path with -1           
    for(j = 0; j < p->recursion_limit - 1; j++){
        sub_params->opt_path[j] = -1;
    }

    *(sfp->cardinality) = 0;
    
    // Compute subforest
    compute_subforest(
                p->forest, 
                p->subtree_memory[sfp->node],
                p->parents_memory[sfp->node], 
                n_nodes, 
                sub_params->forest, 
                &(sub_params->n_forest), 
                sfp->cardinality
            );  

    mfp_dp_opt(sub_params);

    free(sub_params->forest);
    free(sub_params);

}

void* mfp_dp_opt(void* i_params){

    /*
    Description: Recursive function for MFP (Dynamic programming)
    Outputs: Fills optimal value and optimal path.
    */

    OPT_PARAMS* p = (OPT_PARAMS*) i_params;
    nsize_t n_nodes = p->tree->n_nodes;

    // Base case when the path is longer than he recursion limit (n leaves)
    // A path can't be longer than the number of leaves
    if(p->recursion_limit <= 0){
        *(p->opt_value) = 0;
        return NULL;
    }

    // Base case when any other node can be saved
    if(p->n_forest == 0){
        *(p->opt_value) = 0;
        return NULL;
    }

    // Base case when time doesn't allow more saved nodes
    if(p->current_time > p->height * p->t_propagation){
        *(p->opt_value) = 0;
        return NULL;
    }    


    if(p->memory != NULL){
        msize_t memory_value;
        memory_value = search_in_memory(
                            p->memory, 
                            p->forest, 
                            n_nodes, 
                            p->firefighter_position, 
                            p->current_time
                        );
        
        // Base case if value exists in memory
        if(memory_value >= 0){
            *(p->opt_value) = (nsize_t) memory_value;
            return NULL;
        }
    }

    bool* feasible = (bool*) malloc(n_nodes * sizeof(bool));    
    float* distances = (float*) malloc(n_nodes * sizeof(float));    
    
    // Compute feasible nodes
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
    
    nsize_t i, j;
    nsize_t val, max_val = 0;
    bool* subtree;

    pthread_t* threads = (pthread_t*) malloc(n_nodes * sizeof(pthread_t));
    bool* threaded = (bool*) malloc(n_nodes * sizeof(bool));
    nsize_t* cardinalities = (nsize_t*) malloc(n_nodes * sizeof(nsize_t));
    nsize_t* opt_values = (nsize_t*) malloc(n_nodes * sizeof(nsize_t));

    SF_OPT_PARAMS** sub_params = (SF_OPT_PARAMS**) malloc(n_nodes * sizeof(SF_OPT_PARAMS*));

    for(i = 0; i < n_nodes; i++){
        threaded[i] = 0;
        // Only execute for feasible nodes
        if(feasible[i] == 1){

            // Create params
            sub_params[i] = (SF_OPT_PARAMS*) malloc(sizeof(OPT_PARAMS));
            sub_params[i]->old_params = p;
            sub_params[i]->node = i;
            sub_params[i]->distance = distances[i];
            sub_params[i]->opt_path = (msize_t*) malloc((p->recursion_limit - 1) * sizeof(msize_t));
            sub_params[i]->opt_value = &(opt_values[i]);
            sub_params[i]->cardinality = &(cardinalities[i]);
    
            // Try to execute in a thread            
            // If threaded execution fails, executes sequentially
            if(reserve_thread() == 1){
                if (pthread_create(
                        &threads[i], 
                        NULL, 
                        sf_mfp_dp_opt, 
                        sub_params[i]) != 0){                    
                    free_thread();                    
                    sf_mfp_dp_opt(sub_params[i]);
                }
                else{
                    threaded[i] = 1;
                }
            }
            else{
                sf_mfp_dp_opt(sub_params[i]);
            }
        }
    }

    for(i = 0; i < n_nodes; i++){
        if(feasible[i] == 1){
            if(threaded[i] == 1){
                if(pthread_join(threads[i], NULL) != 0){
                    printf("ERROR: Cannot join thread %d\n", i);
                }

                free_thread();
            }

            val = cardinalities[i] + *(sub_params[i]->opt_value);

            if(val > max_val){
                max_val = val;
                p->opt_path[0] = (msize_t) i;
                
                for(j = 0; j < p->recursion_limit - 1; j++){
                    p->opt_path[j + 1] = sub_params[i]->opt_path[j];
                }

            }

            free(sub_params[i]->opt_path);
            free(sub_params[i]);
        }
    }

    free(sub_params);
    free(threaded);    
    free(opt_values);
    free(cardinalities);
    free(threads);
    free(feasible);
    free(distances);


    if(p->memory != NULL){
        add_to_memory(
                p->memory, 
                max_val, 
                p->forest, 
                p->n_forest, 
                p->firefighter_position, 
                p->current_time
            );        
    }

    *(p->opt_value) = max_val;
}


void mfp_dp_solver(
    TREE tree, 
    int root, 
    int firefighter_position, 
    float t_propagation, 
    int* optimal, 
    int* opt_path,
    int use_memoization
    ){
        
    /*
    Description: Wrapper function for MFP solver
    Outputs: Fills optimal value and optimal path.
    */

    nsize_t n_forest, height, n_leaves;
    bool* forest;
    nsize_t* hops_memory;
    bool** subtree_memory;
    bool** parents_memory;
    nsize_t i, j = 0;
    nsize_t opt = 0;

    height = (nsize_t) tree.height;
    n_leaves = (nsize_t) tree.n_leaves;
    forest = (bool*) malloc(tree.n_nodes * sizeof(bool));
    hops_memory = (nsize_t*) calloc(tree.n_nodes, sizeof(nsize_t));
    subtree_memory = (bool**) calloc(tree.n_nodes, sizeof(bool*));
    parents_memory = (bool**) calloc(tree.n_nodes, sizeof(bool*));


    // Initialize subtrees and parents memory
    for(i = 0; i < tree.n_nodes; i++){
        subtree_memory[i] = (bool*) calloc(tree.n_nodes, sizeof(bool));
        compute_subtree(&tree, i, subtree_memory[i]);

        parents_memory[i] = (bool*) calloc(tree.n_nodes, sizeof(bool));
        compute_parents(&tree, i, parents_memory[i]);

        hops_memory[i] = hops(&tree, root, i, hops_memory);        
    }

    // Initialize subforest
    n_forest = tree.n_nodes - 1;
    
    for(i = 0; i < tree.n_nodes; i++){
        if(i == root){
            forest[i] = 0;
        }
        else{
            forest[i] = 1;
        }
    }

    // Initialize optimal path   
    for(i = 0; i < n_leaves + 1; i++){
        opt_path[i] = -1;
    }
    
    // Initialize memory
    MEMORY* memory = NULL;
    if(use_memoization == 1){
        init_memory(&memory);
    }

    // Initialize threading mutex
    if (pthread_mutex_init(&__mutex_n_threads__, NULL) != 0) { 
        printf("ERROR: Threads mutex init has failed\n"); 
    } 


    // Build starting parameters
    OPT_PARAMS params;
    params.tree = &tree; 
    params.height = height; 
    params.forest = forest; 
    params.n_forest = n_forest; 
    params.root = (nsize_t) root; 
    params.firefighter_position = (nsize_t) firefighter_position; 
    params.current_time = 0; 
    params.t_propagation = t_propagation;
    params.memory = memory; 
    params.hops_memory = hops_memory;
    params.subtree_memory = subtree_memory;
    params.parents_memory = parents_memory;
    params.recursion_limit = n_leaves + 1;
    params.opt_path = (msize_t*) opt_path;
    params.opt_value = &opt;

    // Optimize
    mfp_dp_opt(&params);

    // Validate number of threads consistency
    if(__threads_available__ != MAX_THREADS_ALLOWED){
        printf("ERROR: Something went wrong with threads :S\n");
    }

    pthread_mutex_destroy(&__mutex_n_threads__);


    for(i = 0; i < tree.n_nodes; i++){
        free(subtree_memory[i]);
        free(parents_memory[i]);
    }
    
    if(memory != NULL){
        delete_memory(&memory);
    }

    free(forest);
    free(hops_memory);
    free(subtree_memory);
    free(parents_memory);

    (*optimal) = opt;    
}




/*
*
* Greedy solver
*
*/

typedef struct {
    nsize_t* forest;
    bool* subtree; 
    bool* parents; 
    nsize_t n_nodes; 
    bool* subforest; 
    nsize_t* n_subforest; 
    nsize_t* cardinality;
} SUBFOREST_PARAMS;

void* compute_subforest_thread(void* i_params){

    /*
    Description: Wrapper function for subforest threaded
    Outputs: Fills optimal value and optimal path.
    */

    SUBFOREST_PARAMS* p = (SUBFOREST_PARAMS*) i_params;

    compute_subforest(
                p->forest, 
                p->subtree,
                p->parents, 
                p->n_nodes, 
                p->subforest, 
                p->n_subforest, 
                p->cardinality
            );  
}

void mfp_greedy_solver(TREE tree, int root, int firefighter_position, float t_propagation, int* optimal, int* opt_path){

    /*
    Description: MFP greedy solver
    Outputs: Fills optimal value and optimal path.
    */

    nsize_t n_forest, height, n_leaves;
    nsize_t* hops_memory;
    bool* forest;
    bool* feasible;
    float* distances;
    
    bool** subtree_memory;
    bool** parents_memory;
    nsize_t i, j = 0;

    float current_time;
    float next_time;
    nsize_t current_position;

    int opt = 0;

    height = (nsize_t) tree.height;
    n_leaves = (nsize_t) tree.n_leaves;
    nsize_t n_nodes = (nsize_t) tree.n_nodes;

    forest = (bool*) malloc(n_nodes * sizeof(bool));
    hops_memory = (nsize_t*) calloc(n_nodes, sizeof(nsize_t));
    subtree_memory = (bool**) calloc(n_nodes, sizeof(bool*));
    parents_memory = (bool**) calloc(n_nodes, sizeof(bool*));
    feasible = (bool*) malloc(n_nodes * sizeof(bool));    
    distances = (float*) malloc(n_nodes * sizeof(float)); 

    // Prefilling memory
    for(i = 0; i < n_nodes; i++){
        subtree_memory[i] = (bool*) calloc(n_nodes, sizeof(bool));
        compute_subtree(&tree, i, subtree_memory[i]);

        parents_memory[i] = (bool*) calloc(n_nodes, sizeof(bool));
        compute_parents(&tree, i, parents_memory[i]);

        hops_memory[i] = hops(&tree, root, i, hops_memory);        
    }
    
    n_forest = n_nodes - 1;
    
    for(i = 0; i < n_nodes; i++){
        if(i == root){
            forest[i] = 0;
        }
        else{
            forest[i] = 1;
        }
    }
   
    for(i = 0; i < n_leaves + 1; i++){
        opt_path[i] = -1;
    }
    
    
    if (pthread_mutex_init(&__mutex_n_threads__, NULL) != 0) { 
        printf("ERROR: Threads mutex init has failed\n"); 
    } 

    current_time = 0;
    next_time = 0;
    current_position = (nsize_t) firefighter_position;
    

    nsize_t* cardinalities = (nsize_t*) malloc(n_nodes * sizeof(nsize_t));
    nsize_t* n_subforests = (nsize_t*) malloc(n_nodes * sizeof(nsize_t));

    nsize_t max_cardinality;

    pthread_t* threads = (pthread_t*) malloc(n_nodes * sizeof(pthread_t));
    bool* threaded = (bool*) malloc(n_nodes * sizeof(bool));

    SUBFOREST_PARAMS** sub_params = (SUBFOREST_PARAMS**) malloc(n_nodes * sizeof(SUBFOREST_PARAMS*));

    for(i = 0; i < n_leaves + 1; i++){  

        current_time = next_time;
        
        get_feasible(
                &tree,
                forest, 
                (nsize_t) root, 
                current_position, 
                current_time, 
                t_propagation, 
                feasible, 
                distances, 
                hops_memory
            );

        for(j = 0; j < n_nodes; j++){  
            cardinalities[j] = 0;
            threaded[j] = 0;

            if(feasible[j] == 1){

                sub_params[j] = (SUBFOREST_PARAMS*) malloc(sizeof(SUBFOREST_PARAMS));  

                sub_params[j]->forest = forest;
                sub_params[j]->subtree = subtree_memory[j];
                sub_params[j]->parents = parents_memory[j]; 
                sub_params[j]->n_nodes = n_nodes; 
                sub_params[j]->subforest = (bool*) malloc(n_nodes * sizeof(bool)); 
                sub_params[j]->n_subforest = &(n_subforests[j]); 
                sub_params[j]->cardinality = &(cardinalities[j]);

                
                if(reserve_thread() == 1){
                    if (pthread_create(
                            &threads[j], 
                            NULL, 
                            compute_subforest_thread, 
                            sub_params[j]) != 0){
                        free_thread();                    
                        compute_subforest_thread(sub_params[j]);
                    }
                    else{
                        threaded[j] = 1;
                    }
                }
                else{
                    compute_subforest_thread(sub_params[j]);
                }
            }
        }

        max_cardinality = 0;
        for(j = 0; j < n_nodes; j++){  
            if(feasible[j] == 1){

                if(threaded[j] == 1){
                    if(pthread_join(threads[j], NULL) != 0){
                        printf("ERROR: Cannot join thread %d\n", i);
                    }

                    free_thread();
                }

                if(cardinalities[j] > max_cardinality){
                    max_cardinality = cardinalities[j];
                    current_position = j;
                    next_time = current_time + distances[j];
                    memcpy(forest, sub_params[j]->subforest, n_nodes * sizeof(bool));
                    opt_path[i] = j;
                }

                free(sub_params[j]->subforest);
                free(sub_params[j]);
            }
        }

        opt += max_cardinality;
    }


    if(__threads_available__ != MAX_THREADS_ALLOWED){
        printf("ERROR: Something went wrong with threads :S\n");
    }

    pthread_mutex_destroy(&__mutex_n_threads__);

    for(i = 0; i < tree.n_nodes; i++){
        free(subtree_memory[i]);
        free(parents_memory[i]);
    }

    free(sub_params);
    free(threads);
    free(threaded);
    free(cardinalities);
    free(n_subforests);
    free(feasible);
    free(distances);
    free(forest);
    free(hops_memory);
    free(subtree_memory);
    free(parents_memory);

    (*optimal) = opt;
    
}
