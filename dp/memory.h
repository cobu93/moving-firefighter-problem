#include <pthread.h> 

typedef unsigned int uint;
typedef long long int lli;

pthread_mutex_t mutex; 

typedef struct
{
    lli key;
    lli value;
} HASH_ITEM;

typedef struct
{
    HASH_ITEM* items;
    lli size;
} MAP;

typedef struct {  
    MAP* forest_map;
    MAP* location_map;
    MAP* time_map;

    int*** storage;
} MEMORY;


lli get_map_value(MAP* map, lli key){
    int i;
    
    for(i = 0; i < map->size; i++){
        if(map->items[i].key == key){
            return map->items[i].value;
        }
    }

    return -1;
}

lli add_map_value(MAP* map, lli key, lli value){
    int i;

    HASH_ITEM* items_backup;
    items_backup = (HASH_ITEM*) malloc(map->size * sizeof(HASH_ITEM));
    memcpy(items_backup, map->items, map->size * sizeof(HASH_ITEM));
    free(map->items);

    map->size++;
    map->items = (HASH_ITEM*) malloc(map->size * sizeof(HASH_ITEM));
    memcpy(map->items, items_backup, (map->size - 1) * sizeof(HASH_ITEM));
    free(items_backup);

    map->items[map->size - 1].key = key;
    map->items[map->size - 1].value = value;

    return value;
}

void free_storage(int**** storage, long long int  dim1, long long int  dim2, long long int dim3){
    long long int i, j, k;
    
    for(i = 0; i < dim1; i++){        
        for(j = 0; j < dim2; j++){
            free((*storage)[i][j]);
        }
        free((*storage)[i]);
    }

    if(dim1 > 0){
        free(*storage);
    }
    
}

lli location_to_key(int location){
    return (lli) location;
}

lli forest_to_key(uint* forest, uint n_nodes){
    lli key = 0;
    uint i;

    for(i = 0; i < n_nodes; i++){
        if(forest[i] > 0){
            key = key | (0x01 << i);
        }
    }

    return key;
}

lli time_to_key(float time){
    lli key = 0;
    memcpy(&key, &time, sizeof(float));
    return key;
}

int search_in_memory(MEMORY* memory, uint* forest, uint n_nodes, int location, float time){
    lli location_idx, forest_idx, time_idx;
    lli location_key, forest_key, time_key;

    location_key = location_to_key(location);
    location_idx = get_map_value(memory->location_map, location_key);
    if(location_idx < 0){ return -1; }

    forest_key = forest_to_key(forest, n_nodes);
    forest_idx = get_map_value(memory->forest_map, forest_key);
    if(forest_idx < 0){ return -1; }

    time_key = time_to_key(time);
    time_idx = get_map_value(memory->time_map, time_key);
    if(time_idx < 0){ return -1; }

    return memory->storage[location_idx][forest_idx][time_idx];
}

int add_to_memory(MEMORY* memory, uint value, uint* forest, uint n_nodes, int location, float time){
    
    lli location_idx, forest_idx, time_idx;
    lli location_key, forest_key, time_key;
    lli location_size, forest_size, time_size;
    uint i, j, k;

    pthread_mutex_lock(&mutex);

    location_size = memory->location_map->size;
    location_key = location_to_key(location);
    location_idx = get_map_value(memory->location_map, location_key);
    if(location_idx < 0){
        location_idx = add_map_value(memory->location_map, location_key, memory->location_map->size);
    }

    forest_size = memory->forest_map->size;
    forest_key = forest_to_key(forest, n_nodes);
    forest_idx = get_map_value(memory->forest_map, forest_key);
    if(forest_idx < 0){
        forest_idx = add_map_value(memory->forest_map, forest_key, memory->forest_map->size);
    }

    time_size = memory->time_map->size;
    time_key = time_to_key(time);
    time_idx = get_map_value(memory->time_map, time_key);
    if(time_idx < 0){
        time_idx = add_map_value(memory->time_map, time_key, memory->time_map->size);
    }

    int*** backup;

    backup = realloc(memory->storage, memory->location_map->size * sizeof(int**));
    if(backup == NULL){ 
        printf("Memory allocation failed!\n"); 
        return -1; 
    }        
    memory->storage = backup;

    if(location_size != memory->location_map->size){
        memory->storage[memory->location_map->size - 1] = NULL;
    }

    for(i = 0; i < memory->location_map->size; i++){
        int** tmp = realloc(memory->storage[i], memory->forest_map->size * sizeof(int*));
        if(tmp == NULL){ 
            printf("Memory allocation failed!\n"); 
            return -1; 
        }               
        memory->storage[i] = tmp;

        if(forest_size != memory->forest_map->size){
            memory->storage[i][memory->forest_map->size - 1] = NULL;
        }        
    }

    if(location_size != memory->location_map->size){
        for(j = 0; j < memory->forest_map->size; j++){
            memory->storage[memory->location_map->size - 1][j] = NULL;
        }
    }
    
    for(i = 0; i < memory->location_map->size; i++){
        for(j = 0; j < memory->forest_map->size; j++){

            int initialize = (memory->storage[i][j] == NULL);

            int* tmp = realloc(memory->storage[i][j], memory->time_map->size * sizeof(int));
            if(tmp == NULL){ 
                printf("Memory allocation failed!\n"); 
                return -1; 
            }
            memory->storage[i][j] = tmp;

            if(initialize){
                for(k = 0; k < memory->time_map->size; k++){
                    memory->storage[i][j][k] = -1;
                }
            }
            else{
                if(time_size != memory->time_map->size){
                    memory->storage[i][j][memory->time_map->size - 1] = -1;
                }
            }
        }
    }

    memory->storage[location_idx][forest_idx][time_idx] = value;
    pthread_mutex_unlock(&mutex);
    return 0;
}


void init_memory(MEMORY** memory){
    
    if (pthread_mutex_init(&mutex, NULL) != 0) { 
        printf("Mutex init has failed\n"); 
    } 

    *memory = (MEMORY*) malloc(sizeof(MEMORY));

    (*memory)->forest_map = (MAP*) malloc(sizeof(MAP));
    (*memory)->location_map = (MAP*) malloc(sizeof(MAP));
    (*memory)->time_map = (MAP*) malloc(sizeof(MAP));

    (*memory)->forest_map->items = (HASH_ITEM*) malloc(sizeof(HASH_ITEM));
    (*memory)->location_map->items = (HASH_ITEM*) malloc(sizeof(HASH_ITEM));
    (*memory)->time_map->items = (HASH_ITEM*) malloc(sizeof(HASH_ITEM));

    (*memory)->forest_map->size = 0;
    (*memory)->location_map->size = 0;
    (*memory)->time_map->size = 0;

    (*memory)->storage = NULL;
}

void delete_memory(MEMORY** memory){

    pthread_mutex_destroy(&mutex);

    free_storage(&((*memory)->storage), (*memory)->location_map->size, (*memory)->forest_map->size, (*memory)->time_map->size);

    free((*memory)->forest_map->items);
    free((*memory)->location_map->items);
    free((*memory)->time_map->items);

    free((*memory)->forest_map);
    free((*memory)->location_map);
    free((*memory)->time_map);

    free(*memory);
}
