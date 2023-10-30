#ifndef LISTA

#define LISTA
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct {  
   int completado;
   int order;     
   MPI_Request *request_send;
   double *solutionSend;    
} typeofelementL;

struct nodeL {
        struct nodeL *next;
        struct nodeL *prev;        
        typeofelementL *element;
};

typedef struct nodeL *positionL;

struct list_buffer {
        int length;
        int number_solutions;
        int size_solutions;
        positionL *init;
        positionL *end;
};

typedef struct list_buffer *list;


void create_list(list *, int, int);
void destroy_list(list *);
unsigned isEmpty(list);
unsigned length(list);
void add_node(list *, positionL*, int );
void delete_node(list *, positionL*);
void init_list(list *, int);
positionL return_pos_n_migration(list, int);
double* return_reception_node(list *, positionL*);
int n_number_of_nodes_to_replace(list *);
unsigned length(list);
positionL int_to_positionL(list, int);
MPI_Request* get_request_node(list *, positionL*);
void update_list(list *);

//#endif
#endif


