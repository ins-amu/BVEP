/*
 * File:   dynamic_list.c
 * Author: david
 *
 * Created on 14 de enero de 2013, 16:00
 */
//#ifdef MPI2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


#ifdef MPI2
#include <mpi.h>
#include <dynamic_list.h>

void create_list(list *l, int number_solutions, int size_solutions) {
    (*l) = NULL;
    (*l) = (list ) malloc(1 * sizeof (struct list_buffer));
    (*l)->init = (positionL *) malloc(1 * sizeof (struct nodeL ));
    (*l)->end = (*l)->init;
    ((struct nodeL *)(*l)->end)->next = NULL;
    (*l)->length = 0;
    (*l)->number_solutions = number_solutions;
    (*l)->size_solutions = size_solutions;
}



int get_length(list *l) {
    return (*l)->length;
}

unsigned isEmpty(list l){
    if (l->length == 0)
        return 1;
    return 0;
}




void destroy_list(list *l) {
    positionL *p1;
    MPI_Request *reqM;
    p1 = (*l)->init;
    int flag;

    while ((positionL *) ((struct nodeL *) p1)->next != NULL) {
        reqM = get_request_node(l, p1);
        flag = 0;
	while (flag == 0) {
 	    MPI_Cancel( reqM);
            MPI_Wait(reqM, MPI_STATUS_IGNORE);
            MPI_Test(reqM, &flag, MPI_STATUS_IGNORE);
	}
        delete_node(l, p1);
    }

    if ((*l)->init != NULL) {
        free((*l)->init);
        p1 = NULL;
    }

    
    free(*l);
    *l = NULL;
}


unsigned length(list l) {
    return (l->length);
}

void add_node(list *l, positionL *p, int mig) {
    positionL *q;
    if (p != NULL) {
        q = (positionL *) ((struct nodeL *) p)->next;
        ((struct nodeL *) p)->next = (struct nodeL *) malloc(1 * sizeof (struct nodeL));
        (((struct nodeL *) p)->next)->element = (typeofelementL *) malloc(1 * sizeof (typeofelementL));

        (((struct nodeL *) p)->next)->element->completado = 0;
        (((struct nodeL *) p)->next)->element->order = mig;
        
        (((struct nodeL *) p)->next)->element->solutionSend = (double *) malloc( (double) ((*l)->number_solutions * (*l)->size_solutions) *  sizeof (double));
        (((struct nodeL *) p)->next)->element->request_send = (MPI_Request *) malloc(1 *  sizeof (MPI_Request));

        (((struct nodeL *) ((struct nodeL *) p))->next)->next = (struct nodeL *) q;
        if ((((struct nodeL *) ((struct nodeL *) p))->next)->next == NULL) 
            (*l)->end = (positionL *) ((struct nodeL *) p)->next;

        (*l)->length = (*l)->length + 1;
    }


}

void delete_node(list *l, positionL *p) {

    positionL *q;
    double *solSEND;
    MPI_Request *req;
    typeofelementL *el;
            

    q = (positionL *) ((struct nodeL *) p)->next;

    solSEND =  (((struct nodeL *) q)->element->solutionSend);
    req = (((struct nodeL *) q)->element->request_send);
    el = ((struct nodeL *) q)->element;
    
    ((struct nodeL *) p)->next = ((struct nodeL *) q)->next;

    if (((struct nodeL *) p)->next == NULL)
        (*l)->end = p;

    (*l)->length = (*l)->length - 1;
    
    
    if (req != NULL) {
        free(req);
        req = NULL;
    }
    
    if ( solSEND != NULL) {
        free(solSEND);
        solSEND = NULL;
    }
    
    if (el != NULL) {
        free(el);
        el = NULL;
    }
    
    if ( q != NULL ) {
        free(q);
        q = NULL;
    }

}




MPI_Request *  get_request_node(list *l, positionL *pos){
    
    if (pos == NULL) {
       return ((struct nodeL *)(*l)->init)->next->element->request_send;
    }
    else {
       return ((struct nodeL *)pos)->next->element->request_send;	
    }
    
}




double * return_reception_node(list *l, positionL *p1){ 

    if (p1 == NULL) 
        return ((struct nodeL *)(*l)->init)->next->element->solutionSend;
    else 
        return ((struct nodeL *) p1)->next->element->solutionSend;
}





void update_list(list *l) {
    positionL *p1;
    int p1_int, p2_int;
    
    p1_int = 0;
    p2_int = (*l)->length;
    
    p1 = (*l)->init;
   
    while (p1_int < p2_int ) {
        if (((struct nodeL *) p1)->next->element->completado == 1) {
            delete_node(l, p1);
        }
        else {
            p1 = (positionL *) ((struct nodeL *) p1)->next;
        }
        p1_int++;
    }
}
 
#endif
