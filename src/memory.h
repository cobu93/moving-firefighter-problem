#include <pthread.h> 
#include <math.h>
#include "tdefinitions.h"

pthread_mutex_t __mutex_memory__; 

typedef struct
{
    char* key;
    uint key_size;
    nsize_t value;
} HASH_ITEM;

typedef struct NODE_{
    HASH_ITEM* item;
    struct NODE_* next_node;
    struct NODE_* previous_node;
} NODE;

typedef struct {  
    NODE* first;  
    NODE* last;
} MEMORY;

float round_decimals(float val, uint decimals){
    float factor = pow(10, decimals);
    return roundf(val * factor) / factor;
}


uint get_key(nsize_t location, bool* forest, float time, nsize_t n_nodes, char** key){

    uint time_int = (int) ((time + 1) * pow(10, 4)); 
    uint time_size = (int)log10(time_int) + 1;
    uint location_size = (int)log10(location + 1) + 1;
    nsize_t i;

    uint key_size = location_size + n_nodes + time_size + 1;
    (*key) = (char*) malloc(key_size * sizeof(char));
    sprintf(*key, "%d%s%d", location, (char*) forest, time_int);

    return key_size;
}

msize_t search_in_memory(MEMORY* memory, bool* forest, nsize_t n_nodes, nsize_t location, float time){
    char* key;
    uint key_size = get_key(location, forest, time, n_nodes, &key);
    
    pthread_mutex_lock(&__mutex_memory__);
    NODE* left = memory->first;
    NODE* right = memory->last;
    NODE* last_left = NULL;
    NODE* last_right = NULL;
    msize_t result = -1;
    
    while(result < 0 && left != right && left != last_right && right != last_left){

        if(key_size == left->item->key_size){
            if(memcmp(left->item->key, key, key_size * sizeof(char)) == 0){
                result = left->item->value;
                break;
            }
        }

        if(key_size == right->item->key_size){
            if(memcmp(right->item->key, key, key_size * sizeof(char)) == 0){
                result = right->item->value;
                break;
            }
        }

        last_left = left;
        last_right = right;

        left = left->next_node;
        right = left->previous_node;
    }
        
    free(key);
    pthread_mutex_unlock(&__mutex_memory__);

    if(result >= 0){
        printf("Smooth baby... Found in memory\n");
    }

    return result;    
}

int add_to_memory(MEMORY* memory, nsize_t value, bool* forest, uint n_nodes, int location, float time){

    pthread_mutex_lock(&__mutex_memory__);
    NODE* node = (NODE*) malloc(sizeof(NODE));
    node->item = (HASH_ITEM*) malloc(sizeof(HASH_ITEM));

    uint key_size = get_key(location, forest, time, n_nodes, &(node->item->key));
    node->item->key_size = key_size;
    node->item->value = value;
    node->previous_node = NULL;
    node->next_node = memory->first;
    memory->first->previous_node = node;
    memory->first = node;
    pthread_mutex_unlock(&__mutex_memory__);

    return 0;
}


void init_memory(MEMORY** memory){
    
    if (pthread_mutex_init(&__mutex_memory__, NULL) != 0) { 
        printf("ERROR: Mutex init has failed\n"); 
    } 

    *memory = (MEMORY*) malloc(sizeof(MEMORY));

    (*memory)->first = (NODE*) malloc(sizeof(NODE));
    (*memory)->first->item = (HASH_ITEM*) malloc(sizeof(HASH_ITEM));
    (*memory)->first->item->key = NULL;
    (*memory)->first->item->key_size = 0;
    
    (*memory)->last = (NODE*) malloc(sizeof(NODE));
    (*memory)->last->item = (HASH_ITEM*) malloc(sizeof(HASH_ITEM));
    (*memory)->last->item->key = NULL;
    (*memory)->last->item->key_size = 0;

    (*memory)->first->next_node = (*memory)->last;
    (*memory)->first->previous_node = NULL;

    (*memory)->last->next_node = NULL;
    (*memory)->last->previous_node = (*memory)->first;
}

void delete_memory(MEMORY** memory){

    pthread_mutex_destroy(&__mutex_memory__);

    NODE* left = (*memory)->first;
    NODE* right = (*memory)->last;
    NODE* next_left = NULL;
    NODE* next_right = NULL;
    
    while(1){

        next_left = left->next_node;
        next_right = right->previous_node;

        free(left->item->key);
        free(left->item);
        free(left);

        if(left == right){
            break;
        }

        free(right->item->key);
        free(right->item);
        free(right);

        if(next_left == right && next_right == left){
            break;
        }
        
        left = next_left;
        right = next_right;       

    }
}
