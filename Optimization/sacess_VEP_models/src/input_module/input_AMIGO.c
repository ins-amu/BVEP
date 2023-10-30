#include <structure_paralleltestbed.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "AMIGO_model.h"
#include "AMIGO_problem.h"
#include "AMIGO_pe.h"
#include <hdf5.h>

void insert_group_dataset_int_value(hid_t group_id, const char *dataname,  int *value ) {
        hid_t dataset_id, driver;
        double value_d;
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif  
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver, &value_d);
        *value = (int) value_d;
}

void insert_group_dataset_double_value(hid_t group_id, const char *dataname,  double *value, herr_t  *status) {
        hid_t dataset_id, driver;
        double value_d;
         
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif  
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        *status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver, &value_d);
        *value = value_d;
}


void insert_group_dataset_vector_double_value(hid_t group_id, const char *dataname,  double *value ) {
        hid_t dataset_id, space,driver;
        double *value_aux;
        hsize_t dim[2];
        int ndims, i;
        herr_t      status;
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif  
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);

        value_aux = (double *) malloc (dim[0]*dim[1]* sizeof (double));
        status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver,value_aux);   
        for (i=0;i<dim[0]; i++) {
            value[i] = value_aux[i];
        }
        
        free(value_aux);
        value_aux = NULL;
}


void insert_group_dataset_vector_double_value_reference(hid_t group_id, const char *dataname,  double *value, int exp_num ) {
        hid_t dataset_id, space,driver;
        double *value_aux;
        hobj_ref_t  *rdata;  
        hsize_t dim[2];
        int ndims, i;
        herr_t      status;
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif  
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        rdata = (hobj_ref_t *) malloc (dim[0]*dim[1]* sizeof (hobj_ref_t));
        status = H5Dread (dataset_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, driver,rdata); 
        dataset_id = H5Rdereference (dataset_id, H5R_OBJECT, &rdata[exp_num]);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        value_aux = (double *) malloc(dim[0]*dim[1]*sizeof(double));
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver,value_aux );
        free(rdata);
        rdata = NULL;
        for (i=0;i<dim[0]; i++) {
            value[i] = value_aux[i];
        }
        free(value_aux);
        value_aux = NULL;
        H5Dclose(dataset_id);
        
}



void insert_group_dataset_matrix_double_value_reference(hid_t group_id, const char *dataname,  double **value, int exp_num,  int mi, int mj) {
        hid_t dataset_id, space,driver;
        double *value_aux;
        hobj_ref_t  *rdata;  
        hsize_t dim[2];
        int ndims, i, j;
        herr_t      status;
        int counter;
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif  
    
        counter = 0;
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        rdata = (hobj_ref_t *) malloc (dim[0]*dim[1]* sizeof (hobj_ref_t));
        status = H5Dread (dataset_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, driver,rdata); 
        dataset_id = H5Rdereference (dataset_id, H5R_OBJECT, &rdata[exp_num]);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        value_aux = (double *) malloc(dim[0]*dim[1]*sizeof(double));
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver,value_aux );
        free(rdata);
        rdata = NULL;
        for (i=0;i<mi; i++) {
            for (j=0;j< mj; j++) {
                value[i][j] = value_aux[counter];
                counter++;
            }
        }
        free(value_aux);
        value_aux = NULL;
        H5Dclose(dataset_id);
}



int return_size___group_dataset_matrix_double_value_reference(hid_t group_id, const char *dataname,  double **value, int exp_num) {
        hid_t dataset_id, space, driver;
        hobj_ref_t  *rdata;  
        hsize_t dim[2];
        int ndims;
        herr_t      status;
        
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif  
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        rdata = (hobj_ref_t *) malloc (dim[0]*dim[1]* sizeof (hobj_ref_t));
        status = H5Dread (dataset_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, driver,rdata); 
        dataset_id = H5Rdereference (dataset_id, H5R_OBJECT, &rdata[exp_num]);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        
        
        free(rdata);
        rdata = NULL;
        H5Dclose(dataset_id);
        return dim[0]*dim[1];
}




void insert_group_dataset_vector_integer_value2(hid_t group_id, const char *dataname,  int *value ) {
        hid_t dataset_id, space, driver;
        double *value_aux; 
        hsize_t dim[2];
        int ndims, i;
        herr_t      status;
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif   
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        value_aux = (double *) malloc (dim[0]*dim[1]* sizeof (double));
        status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver,value_aux);   
        for (i=0;i<dim[0]; i++) {
            value[i] =  (int) value_aux[i] - 1;
        }
        
        free(value_aux);
        value_aux = NULL;
}

void insert_group_dataset_vector_integer_value2_reference(hid_t group_id, const char *dataname,  int *value, int exp_num ) {
        hid_t dataset_id, space, driver;
        double *value_aux;
        hobj_ref_t  *rdata;  
        hsize_t dim[2];
        int ndims, i;
        herr_t      status;
        
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif    
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        rdata = (hobj_ref_t *) malloc (dim[0]*dim[1]* sizeof (hobj_ref_t));
        status = H5Dread (dataset_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, driver,rdata);   
        dataset_id = H5Rdereference (dataset_id, H5R_OBJECT, &rdata[exp_num]);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        value_aux = (double *) malloc(dim[0]*dim[1]*sizeof(double) );
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver, value_aux); 
        for (i=0;i<dim[0]; i++) {
            value[i] =  (int) ( (int) value_aux[i]) - 1;
        }
        free(rdata);
        rdata = NULL;        
        free(value_aux);
        value_aux = NULL;
}


void insert_group_dataset_vector_double_value2(hid_t group_id, const char *dataname,  double *value ) {
        hid_t dataset_id, space, driver;
        double *value_aux;
        hsize_t dim[2];
        int ndims, i;
        herr_t      status;
        
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif       
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        //printf ("%d - %d\n", dim[0], dim[1]);
        value_aux = (double *) malloc (dim[0]*dim[1]* sizeof (double));
        status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver,value_aux);   
        for (i=0;i<dim[0]; i++) {
            value[i] =  value_aux[i] - 1;
        }
        
        free(value_aux);
        value_aux = NULL;
}


void insert_group_dataset_int_value_reference(hid_t group_id, const char *dataname,  int *value, int exp_num ) {
        hid_t dataset_id, space, driver;
        double value_d;
        hobj_ref_t  *rdata;  
        hsize_t dim[2];
        int ndims;
        herr_t      status;
        
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif       
    
        dataset_id = H5Dopen2(group_id,dataname, driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        rdata = (hobj_ref_t *) malloc (dim[0]*dim[1]* sizeof (hobj_ref_t));
        status = H5Dread (dataset_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, driver,rdata);   
        dataset_id = H5Rdereference (dataset_id, H5R_OBJECT, &rdata[exp_num]);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver, &value_d);
        free(rdata);
        rdata = NULL;
        *value = (int) value_d;
        
}


AMIGO_model* hdf5AllocateAMIGOmodel(hid_t privstruct_id,hid_t inputs_id, int exp_num){

	AMIGO_model* amigo_model;
	int n_states ,n_observables,n_pars,n_opt_pars,n_times,n_opt_ics,n_controls,n_controls_t, index_observables_dim[2],i;
        double *index_observables_d, max_num_steps_d, max_error_test_fails_d;
	hid_t  pesol_id;
        hid_t model_id;
        hid_t  exp_id;
        hid_t  ivpsol_id;
        hid_t  nlpsol_id;
        hid_t  dataset_id;
        hid_t  space, driver;
        hsize_t     dim[2];
        herr_t      status, status2;
        hobj_ref_t  *rdata;   
        int         ndims,size_index;
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif
    
        model_id = H5Gopen2(inputs_id, "model/", driver);
        pesol_id = H5Gopen2(inputs_id, "PEsol/", driver);
        exp_id = H5Gopen2(inputs_id, "exps/", driver);
        ivpsol_id = H5Gopen2(inputs_id, "ivpsol/", driver);
        nlpsol_id = H5Gopen2(inputs_id, "nlpsol/", driver);
        
        insert_group_dataset_int_value(model_id, "n_st", &n_states);
        insert_group_dataset_int_value_reference(privstruct_id,"n_obs",&n_observables,exp_num);
        insert_group_dataset_int_value(model_id, "n_par", &n_pars);
        insert_group_dataset_int_value(pesol_id, "n_theta", &n_opt_pars);
        insert_group_dataset_int_value_reference(privstruct_id,"n_s",&n_times,exp_num);
        insert_group_dataset_int_value_reference(pesol_id,"n_local_theta_y0",&n_opt_ics,exp_num);
        insert_group_dataset_int_value(model_id, "n_stimulus", &n_controls);
        
        
        dataset_id = H5Dopen2(privstruct_id, "t_con", driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        rdata = (hobj_ref_t *) malloc (dim[0]*dim[1]* sizeof (hobj_ref_t));
        status = H5Dread (dataset_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, driver,rdata);   
        dataset_id = H5Rdereference (dataset_id, H5R_OBJECT, &rdata[exp_num]);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        n_controls_t = dim[0];
        free(rdata);
        rdata = NULL;
        
        
        amigo_model=allocate_AMIGO_model(n_states,n_observables,n_pars,
		n_opt_pars,n_times,n_opt_ics,n_controls, n_controls_t, exp_num);

        dim[0]=-1;
        dim[1]=-1;
        dataset_id = H5Dopen2(privstruct_id, "index_observables", driver);
        space = H5Dget_space (dataset_id);
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        
        rdata = (hobj_ref_t *) malloc (dim[0]*dim[1]* sizeof (hobj_ref_t));
        status = H5Dread (dataset_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, driver,rdata);   
        dataset_id = H5Rdereference (dataset_id, H5R_OBJECT, &rdata[exp_num]);
        space = H5Dget_space (dataset_id);
        
        ndims = H5Sget_simple_extent_dims (space, dim, NULL);
        
        
        
        index_observables_d = (double *) malloc( dim[0] * dim[1] *sizeof(double));
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver, index_observables_d);        
        
        index_observables_dim[0] = dim[1];        
        index_observables_dim[1] = dim[0];
        free(rdata);
        rdata = NULL;
	        
        size_index = 0;
        for (i = 0; i < n_observables; i++){
            size_index = size_index + index_observables_d[i];
        }
        
        //printf("2index_observables_dim %d - %d\n",dim[0], dim[1]);
        if( (index_observables_dim[0]<n_observables  && index_observables_dim[1]<n_observables) || 
            size_index == 0 ){
		amigo_model->use_obs_func=1;
		amigo_model->use_sens_obs_func=1;

	}else{

		for (i = 0; i < n_observables; i++){
			amigo_model->index_observables[i]=index_observables_d[i]-1;
                        //printf("index_observables_d %d - %d\n",i, index_observables_d[i] );
		}
	}        
        
        
        free(index_observables_d);
        index_observables_d = NULL;
            
        if ( n_pars > 0)
                insert_group_dataset_vector_double_value(model_id,"par",&amigo_model->pars[0] );
        if (n_opt_pars>0) {
                insert_group_dataset_vector_integer_value2(pesol_id,"index_global_theta", &(amigo_model->index_opt_pars[0]) );
                insert_group_dataset_vector_double_value(pesol_id,"global_theta_guess", &(amigo_model->pars_guess[0]) );
                insert_group_dataset_vector_double_value(pesol_id,"global_theta_min", &(amigo_model->pars_LB[0]) );
                insert_group_dataset_vector_double_value(pesol_id,"global_theta_max", &(amigo_model->pars_UB[0]) );
        }
  
          

        insert_group_dataset_vector_double_value_reference(privstruct_id, "t_in", &(amigo_model->t0 ), exp_num);
        
        insert_group_dataset_vector_double_value_reference(privstruct_id, "t_f", &(amigo_model->tf), exp_num);

        
        if (n_times>0)
                insert_group_dataset_vector_double_value_reference(privstruct_id, "t_int", &(amigo_model->t[0]), exp_num);
        
 
        
        if (n_states>0)
                insert_group_dataset_vector_double_value_reference(exp_id, "exp_y0", &(amigo_model->y0[0]), exp_num);
        
        
        if (n_opt_ics>0) {
                insert_group_dataset_vector_integer_value2_reference(pesol_id,"index_local_theta_y0", &(amigo_model->index_opt_ics[0]), exp_num );
                insert_group_dataset_vector_double_value_reference(pesol_id,"local_theta_y0_guess", &(amigo_model->y0_guess[0]), exp_num );
                insert_group_dataset_vector_double_value_reference(pesol_id,"local_theta_y0_min", &(amigo_model->y0_LB[0]), exp_num );
                insert_group_dataset_vector_double_value_reference(pesol_id,"local_theta_y0_max", &(amigo_model->y0_UB[0]), exp_num );
        }

        if (n_controls_t>0) 
                insert_group_dataset_vector_double_value_reference(privstruct_id,"t_con", &(amigo_model->controls_t[0]), exp_num );
        if ( (n_controls_t>=0) && (n_controls>=0) ) {
                insert_group_dataset_matrix_double_value_reference(privstruct_id, "u", amigo_model->controls_v, exp_num, n_controls, n_controls_t-1 );
                insert_group_dataset_matrix_double_value_reference(privstruct_id, "pend", amigo_model->slope, exp_num, n_controls, n_controls_t-1 );
        }
        
        if (n_observables>0) {
            if (return_size___group_dataset_matrix_double_value_reference > 0) {
                insert_group_dataset_matrix_double_value_reference(privstruct_id, "exp_data", amigo_model->exp_data, exp_num, n_observables, n_times );           
            }
        }
        
        
        insert_group_dataset_matrix_double_value_reference(exp_id, "Q", amigo_model->Q, exp_num, n_observables, n_times);
        insert_group_dataset_vector_double_value_reference(privstruct_id, "w_obs", amigo_model->w_obs, exp_num);
        
        //Simulation Related Parameter
        
        insert_group_dataset_double_value(ivpsol_id, "rtol",  &(amigo_model->reltol), &status2 );
        if (status2 < 0 ) {
            amigo_model->reltol=1e-6;
        }
        insert_group_dataset_double_value(ivpsol_id, "atol",  &(amigo_model->atol), &status2 );
        if (status2 < 0 ) {
            amigo_model->atol=1e-6;
        }   
        insert_group_dataset_double_value(ivpsol_id, "max_step_size",  &(amigo_model->max_step_size), &status2 );
        if (status2 < 0 ) {
            amigo_model->max_step_size=DBL_MAX;
        }        
        insert_group_dataset_double_value(ivpsol_id, "ivp_maxnumsteps",  &max_num_steps_d, &status2 );
        if (status2 < 0 )  {
            amigo_model->max_num_steps = 100000;
        }
        else {
            amigo_model->max_num_steps=(int) max_num_steps_d ;
        }
        //printf("amigo_model->max_num_steps %d\n", amigo_model->max_num_steps);
	//printf("amigo_model->max_step_size %lf\n", amigo_model->max_step_size);
        //////
        insert_group_dataset_double_value(ivpsol_id, "max_error_test_fails",  &max_error_test_fails_d, &status2 );
        if (status2 < 0 )  {
            amigo_model->max_error_test_fails = 50;
        }
        else {
            amigo_model->max_error_test_fails=(int) max_error_test_fails_d ;
        }
        //printf("amigo_model->max_error_test_fails %d\n",amigo_model->max_error_test_fails);        
            
        insert_group_dataset_double_value(nlpsol_id, "mkl_tol", &(amigo_model->mkl_tol) , &status2 );   
       // printf("%lf\n",amigo_model->mkl_tol);
        
        
        
        H5Gclose(ivpsol_id);
        H5Gclose(exp_id);
        H5Gclose(model_id);
        H5Gclose(pesol_id);
        H5Gclose(nlpsol_id);
        H5Dclose(dataset_id);
        
	return(amigo_model);
}


AMIGO_problem* openMatFileAMIGO(const char* file){
    hid_t       idfile;  
    hid_t       driver;
    hid_t       group_inputs;
    hid_t       group_privstruct;
    hid_t       dataset_id;
    hid_t       nlpsol_id;
    herr_t      status, status2;
    int i;
    double n_exp_d;
    int n_exp;
    AMIGO_problem* amigo_problem;
    AMIGO_model** amigo_models;
    double cvodes_gradient_d, mkl_gradient_d, iterprint_d;
    const char *file2;
    
    
#ifdef  MPI2
    driver = H5P_DEFAULT;
#else
    driver = H5P_DEFAULT;    
#endif
            
    H5Eset_auto( H5E_DEFAULT, NULL, NULL);
    file2 = (char *) malloc(100*sizeof(char));
    file2 = file;
    
    idfile = H5Fopen((const char*) file2, (unsigned) H5F_ACC_RDONLY, driver);
    
    group_inputs = H5Gopen2(idfile, "/inputs/", driver);
    group_privstruct = H5Gopen2(idfile, "/privstruct/", driver);
    nlpsol_id = H5Gopen2(idfile, "nlpsol/", driver);
    
    dataset_id = H5Dopen2(group_privstruct, "n_exp", driver);
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, driver, &n_exp_d);

    amigo_models=(AMIGO_model**)malloc(sizeof(AMIGO_model*)*n_exp_d);
    n_exp = (int) n_exp_d;

    for (i = 0; i < n_exp; i++){
	amigo_models[i]=hdf5AllocateAMIGOmodel(group_privstruct,group_inputs,i);
    }
    
    amigo_problem=allocate_AMIGO_problem(n_exp,amigo_models);
    amigo_problem->n_exp = n_exp;
    
    insert_group_dataset_double_value(nlpsol_id, "cvodes_gradient", &cvodes_gradient_d, &status2);
    if (status2 < 0) {
        amigo_problem->cvodes_gradient = 1;
    } else {
        amigo_problem->cvodes_gradient = (int) cvodes_gradient_d;
    }
    
    insert_group_dataset_double_value(nlpsol_id, "mkl_gradient", &mkl_gradient_d, &status2);
    if (status2 < 0) {
        amigo_problem->mkl_gradient = 0;
    } else {
        amigo_problem->mkl_gradient = (int) mkl_gradient_d;
    }
    
    insert_group_dataset_double_value(nlpsol_id, "iterprint", &iterprint_d, &status2);
    if (status2 < 0) {
        amigo_problem->verbose = 0;
    } else {
        amigo_problem->verbose = (int) iterprint_d;
    }
    
    
    
    status = H5Fclose(idfile);
    status = H5Gclose(group_inputs);
    status = H5Gclose(group_privstruct);
    status = H5Dclose(dataset_id);
    status = H5Gclose(nlpsol_id);
    
    return(amigo_problem);
}
