MODULE modsacess
    USE iso_c_binding
    USE scattersearchtypes
    USE scattersearchfunctions
    USE localsolver
    USE parallelscattersearchfunctions
   
#ifdef MPI2 
        
CONTAINS

    FUNCTION sacess(exp1, fitnessfunction, results1, maxfunevals, ftarget) RESULT (outresult)
    
        ! Declaración de variables
        USE common_functions
        USE qsort_module
        
        IMPLICIT NONE
        TYPE(local_solver_help_vars) :: local_solver_var
        TYPE(algorithm_common_vars) :: common_vars
        TYPE(ess_common_vars) :: ess_comm
        TYPE(time_ess) :: time
        TYPE(populations) :: pop
        TYPE(master_migration) :: migration_master
        TYPE(C_PTR), INTENT(INOUT) :: exp1, results1
        REAL (C_DOUBLE), INTENT(INOUT) :: ftarget
        INTEGER (C_LONG), INTENT(INOUT) :: maxfunevals
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        INTEGER(C_LONG) :: nfunevaltotal
        INTEGER :: i, j
        TYPE(opts) :: opts1, default1
        TYPE(problem) :: problem1
        TYPE(resultsf) :: results
        INTEGER :: xbest_in_refset, bks_in_population_flag
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xbest,fbest, fbest_lastiter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xl_log, xu_log, randomV
        !INTEGER, INTENT(IN) :: openmp

        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: fbest_new
        INTEGER :: counter
        INTEGER, DIMENSION(:), ALLOCATABLE :: indvect
        INTEGER :: lb_p, auxsize
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: st_p
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: factor, v1, v2, v3, new_comb
        INTEGER(C_INT) :: outresult
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: NAN, INF, threshold_slave
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: aux_val_cooperative
        INTEGER :: tam
        REAL(C_DOUBLE) :: cputime1
        
        !COOPERATIVE VARIABLES
        INTEGER :: flag , cooperativempitestess, OMP_NUM_THREADS
        REAL(C_DOUBLE) :: initmpi, calctimempi
        INTEGER(C_INT) :: getopenmpoption, openmp_pos

!-------------------------------------------------------------------------------------------------------!
! INITIALIZE STRUCTS AND PARAMETERS
!-------------------------------------------------------------------------------------------------------!
        CALL setnan(NAN)
        CALL setdblmax(INF)

        ! [SACESS] INITIALIZE PROBLEM1, OPTS1, COMMON_VARS%NVAR STRUCTS
        CALL problem_specifications(exp1, problem1,opts1,common_vars%nvar,ftarget,maxfunevals) 
        ! INITIALIZE COMMON_VARS STRUCTS
        CALL init_common_vars(exp1, common_vars, opts1)
        ! INITIALIZE TIME STRUCTS
        CALL init_time_vars(time)
        ! CHECKBOUNDS : Check if bounds have the same dimension
        CALL checkbounds(problem1,common_vars%status)
        ! INIT BEST VARIABLES
        CALL initbestvars(problem1,xbest,fbest,common_vars%nvar)
        ! INIT OLD BEST VARIABLES        
        CALL initoldbestvars(common_vars%BKS_x,common_vars%BKS_f,common_vars%nvar)
        ! INIT INT BIN VARIABLES
        CALL initintbinvars(problem1)
        ! INIT RANDOM GENERATOR        
        CALL initrngrandomparallel(exp1,common_vars%idp)
        ! CREATE MIGRATION MASTER: CREATE VARS TO USE IN SELF-ADAPTATION MODE
        CALL create_migration_master(migration_master,common_vars%NPROC) 
        ! CALC THE SIZE OF THE POPULATION
        CALL calc_dim_refset(opts1%globaloptions%dim_refset, common_vars%nvar,&
                opts1%useroptions%iterprint,common_vars%idp,common_vars%NPROC)
        ! [SACESS] IF THE CONFIGURATION IS HETEROGENEOUS, A DIFFERENT CONFIGURATION IS ASIGNED FOR EACH PROCESSROR
        CALL hete_param_eSS2(exp1, problem1, opts1, common_vars%NPROC,opts1%globaloptions%dim_refset, &
                                                                        common_vars%idp,common_vars%nvar)
       ! CALC THE SIZE OF THE NDIVERSE SET
        CALL calc_ndiverse_hyper(opts1%globaloptions%ndiverse, common_vars%nvar, &
                        opts1%useroptions%iterprint, common_vars%idp, opts1%globaloptions%dim_refset)
        ! INITIALIZATION OF OUTPUT VARS
        CALL initoutputvars(exp1)
        ! ENABLE DISTRIBUTED CRITERIA
        CALL setdistcriteria(exp1,1)
        ! CHARGE THE PARAMETERS FOR MPI COOPERATION 
        CALL chargecooperativeparametersfortran(exp1, opts1%globaloptions%dim_refset, tam, common_vars%idp, &
                        ftarget, maxfunevals)  
        ! CHARGE THE PARAMETERS FOR TOPOLOGY COMMUNICATIONS
        CALL createcooperativetopologyess(exp1)
        ! PRINT SLAVE INITIAL INFORMATION 
        CALL slave_information(common_vars%idp, opts1%globaloptions%dim_refset, opts1%globaloptions%ndiverse, 0)
        ! INITIALIZE TIME STRUCTS
        time%starttime = initmpi()
        CALL saveinittime(exp1,time%starttime)
        ! INIT LOCAL SOLVER FORTRAN OPTIONS
        CALL initlocaloptions(opts1)
        ! INIT INEQUELITY CONSTRAINTS FORTRAN OPTIONS 
        CALL init_inequality_constraints(problem1%neq, problem1%CL, problem1%CU)
        problem1%ineq = problem1%neq + problem1%ineq
        ! CALC NCONST VARIABLE
        CALL calcnconst(problem1,common_vars%nconst) 
        ! CHECK OUTPUT OBJECTIVE FUNCTION
        CALL check_output_obj_funct(problem1, opts1, common_vars%status, common_vars%nconst)
        ! PASAMOS A LOGARITMO #SUB INIT_LOG
        ALLOCATE(xl_log(size(problem1%XL)))
        ALLOCATE(xu_log(size(problem1%XU)))
        xl_log = problem1%XL
        xu_log = problem1%XU


        bks_in_population_flag=0
        ! LOGARITHMIC OPTIONS
        if (ALLOCATED(opts1%useroptions%log_var)) then   
            CALL converttolog2(xl_log, common_vars%nvar,opts1%useroptions%log_var)
            CALL converttolog2(xu_log, common_vars%nvar,opts1%useroptions%log_var)
        end if
        ! CREATE LOCAL SEARCH STRUCTS
        CALL build_aux_local_search(problem1,local_solver_var, opts1, exp1)
        ! CALC THE SIZE OF SERIALIZE MESSAGE
        CALL sizeseriallize(common_vars%nvar, common_vars%nconst, common_vars%sizeser)

!-------------------------------------------------------------------------------------------------------!        
! OPEN MASTER ASYNCHRONOUS BUFFER && CREATE WINDOWS    
!-------------------------------------------------------------------------------------------------------!  
        CALL asynchinitmasterandwindows( exp1, common_vars%sizeser )     
!-------------------------------------------------------------------------------------------------------!  
!-------------------------------------------------------------------------------------------------------!        
! CALC INIT SOLUTION SET AND REFSET GLOBAL 
!-------------------------------------------------------------------------------------------------------!     
#ifdef OPENMP
        CALL calc_init_solt_and_refset_global(opts1,exp1,pop,common_vars,xl_log,xu_log,&
                                                        problem1,fitnessfunction,1)
#else
        CALL calc_init_solt_and_refset_global(opts1,exp1,pop,common_vars,xl_log,xu_log,&
                                                        problem1,fitnessfunction,0)
#endif
            
        common_vars%lastevals = common_vars%nfuneval
!-------------------------------------------------------------------------------------------------------!         
! CHECK AND INITIALIZE THE BEST VALUE IN REFSET
!-------------------------------------------------------------------------------------------------------!
        if (common_vars%idp .NE. 0) then
                  xbest = pop%refset%x(:,1)
                  fbest = pop%refset%fpen(1)
                  CALL check_best_value(pop%refset, opts1, xbest, fbest )
        end if
!-------------------------------------------------------------------------------------------------------!
! SUBTOURINE PRINT_INITIAL_POPULATION
!-------------------------------------------------------------------------------------------------------!
        CALL optimization_begins(common_vars%idp)
        if ((opts1%useroptions%iterprint .eq. 1) .AND. (common_vars%idp .NE. 0)) then
                cputime1 = calctimeMPI(exp1,time%starttime)
                WRITE(*, '(A14, I3, A26, I10, A8, D25.10, A10, F12.2, A11, E10.1)') "[PROCESSOR ID=", &
                    common_vars%idp, "] Initial Pop: NFunEvals: ", common_vars%nfuneval, &
                    " Bestf: ", fbest(1), " CPUTime: ", cputime1, " variance: ", variance_vector(pop%refset%fpen)
        endif      
!-------------------------------------------------------------------------------------------------------!
! SUBTOURINE CREATE RESULTS STRUCT AND REORDER POPULATION VARIABLES
!-------------------------------------------------------------------------------------------------------! 
        CALL create_results(results, opts1, xbest, fbest, common_vars%nfuneval, cputime1, common_vars%nvar)
        if (common_vars%idp .NE. 0 ) CALL sort_refset(pop%refset,opts1, indvect)
        if (common_vars%idp .EQ. 0 ) then
            pop%refset%x = NAN
            pop%refset%fpen = NAN
            pop%refset%f = NAN
            pop%refset%penalty = NAN
            if (common_vars%nconst .GT. 0) then
                pop%refset%nlc = NAN         
            end if            
        end if
        ALLOCATE(pop%refset_change(opts1%globaloptions%dim_refset))
        pop%refset_change = 0
!-------------------------------------------------------------------------------------------------------!
! SUBTOURINE INITIALIZE COMBINATION METHOD
!-------------------------------------------------------------------------------------------------------!
        CALL create_combination_scheme_ess(opts1, common_vars, ess_comm, problem1)
!-------------------------------------------------------------------------------------------------------!        
! INIT ASYNCHRONOUS STOPPING CRITERIA         
!-------------------------------------------------------------------------------------------------------!        
        CALL initcooperativestoppingcriteriaess( exp1 )
!-------------------------------------------------------------------------------------------------------!
! INIT OUTPUT VARIABLES
!-------------------------------------------------------------------------------------------------------!
        cputime1 = calctimeMPI(exp1,time%starttime)
        CALL initprintfile(exp1, fbest, common_vars%parallel, common_vars%idp, cputime1, common_vars%nfuneval, 1)    
        CALL printgant(exp1,cputime1,1)
!-------------------------------------------------------------------------------------------------------!         
! INIT : Algorithm main loop
!-------------------------------------------------------------------------------------------------------! 
    do while (common_vars%out_solver .EQ. 0)

!----------------------------------------------------------------------------------------------------------------------------------
! MIGRATION ASYNCHRONOUS REGION            
!----------------------------------------------------------------------------------------------------------------------------------   
!----------------------------------------------------------------------------------------------------------------------------------
!       MASTER CODE            
!----------------------------------------------------------------------------------------------------------------------------------   
          if ( common_vars%idp .EQ. 0 ) then

             CALL asynchronous_master_acess (exp1,opts1,problem1,pop%refset,common_vars%fin,common_vars%nfuneval,fbest,&
                                    common_vars%stopOptimization,common_vars%sizeser, &
                                    common_vars%nconst, common_vars%nvar, time%starttime, migration_master%vector_proc, &
                                    migration_master%ocurrences_send, common_vars  )   

             cputime1 = calctimeMPI(exp1,time%starttime)
             CALL printgant(exp1,cputime1,1)  
!----------------------------------------------------------------------------------------------------------------------------------
!       SLAVE CODE            
!----------------------------------------------------------------------------------------------------------------------------------        
          else
                cputime1 = calctimeMPI(exp1,time%starttime)
                CALL printinititeration(exp1, common_vars%iter, cputime1)

                if ((common_vars%iter .NE. 0) ) then
                    cputime1 = calctimeMPI(exp1,time%starttime)
                   
                    CALL asynchronous_slave_acess (exp1,opts1,problem1,pop%refset,common_vars, common_vars%nconst, &
                           common_vars%BKS_x, common_vars%BKS_f, &
                           time, pop%refset_change, fitnessfunction, xl_log, xu_log, &
                           common_vars%nfuneval, fbest, xbest, & 
                           common_vars%ncounter_slave_recv_sol, common_vars%ncounter_slave_send_sol,&
                           common_vars%lastevals, bks_in_population_flag)
                           
                    CALL adapt_slave(exp1,opts1,problem1,fitnessfunction,common_vars,common_vars%nfuneval,xl_log,xu_log, &
                          ess_comm%hyper_x_L,ess_comm%hyper_x_U,ess_comm%MaxSubSet,ess_comm%MaxSubSet2,ess_comm%ppp,&
                          ess_comm%index1,ess_comm%index2, ess_comm%index,common_vars%nconst, pop%refset,&
                          fbest,xbest,pop%refset_change,local_solver_var,time, &
                          common_vars%ncounter_slave_recv_sol, common_vars%ncounter_slave_send_sol, &
                          common_vars%pending_adaptation,common_vars%lastevals, bks_in_population_flag,1)        

                    CALL select_better_solution(results,opts1,pop%refset,xbest,fbest,local_solver_var%use_bestx, &
                                common_vars%nfuneval,cputime1,common_vars%fin, exp1)

                else 
                    common_vars%BKS_f(1)   = fbest(1)
                    common_vars%BKS_x(:,1) = xbest
                end if

          end if
!----------------------------------------------------------------------------------------------------------------------------------
!       METAHEURISTIC CODE. (1) DELETE DUPLICATED (2) CREATE CHILDSET, MEMBERS_UPDATE, CANDIDATESET AND CANDIDATE_UPDATE
!----------------------------------------------------------------------------------------------------------------------------------           
          if ( common_vars%fin .GT. 0 ) common_vars%out_solver = 1 
          if ( (common_vars%fin .LT. 1) .AND. (common_vars%idp .NE. 0 )) then
            
          CALL create_childset ( pop%childset, common_vars%nvar, CEILING(ess_comm%MaxSubSet2), common_vars%nconst )   
          counter = 1
          if (.not.ALLOCATED(pop%members_update)) then
                ALLOCATE(pop%members_update(opts1%globaloptions%dim_refset))
                pop%members_update =  1
          end if  

          CALL check_duplicated_replace (exp1, fitnessfunction, problem1,pop%refset,opts1, &
                        pop%parents_index1, pop%parents_index2, &
                        xl_log, xu_log, pop%refset_change, pop%members_update, &
                        ess_comm%index1, ess_comm%index2, &
                        common_vars%nfuneval, common_vars%nvar, common_vars%nconst )
    
          CALL create_candidate(pop%candidateset, pop%refset)
          if (.not.ALLOCATED(pop%candidate_update)) then
                ALLOCATE(pop%candidate_update(opts1%globaloptions%dim_refset))
                pop%candidate_update = 0
          end if

!----------------------------------------------------------------------------------------------------------------------------------
! COMBINATION SCATTER SEARCH
!----------------------------------------------------------------------------------------------------------------------------------
          CALL combine_scatter_search(exp1, opts1, common_vars, ess_comm, pop, new_comb)

          !print *, opts1%globaloptions%n_stuck
!----------------------------------------------------------------------------------------------------------------------------------
! EVALUATION SCATTER SEARCH
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef OPENMP
          CALL eval_solution(exp1,opts1,common_vars,CEILING(ess_comm%MaxSubSet2),pop,&
                                        fitnessfunction,problem1,1,counter,new_comb,ess_comm%index)
#else
          CALL eval_solution(exp1,opts1,common_vars,CEILING(ess_comm%MaxSubSet2),pop,&
                                        fitnessfunction,problem1,0,counter,new_comb,ess_comm%index)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! CHECK PARALLEL STOPPING CRITERIA  SUBROUTINE           
!----------------------------------------------------------------------------------------------------------------------------------              
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results, opts1, pop%refset, xbest, fbest,local_solver_var%use_bestx, &
                        common_vars%nfuneval, cputime1, common_vars%fin, exp1)
            CALL printiterationcess(exp1, fbest(1), common_vars%nfuneval,cputime1, 0, 0)
            CALL stoppingcriteriamodule(exp1,common_vars,opts1,results,fbest,xbest,time)
!----------------------------------------------------------------------------------------------------------------------------------
!       METAHEURISTIC CODE. (3) BEYOND METHOD (4) REMOVE STUCK MEMEBRS 
!----------------------------------------------------------------------------------------------------------------------------------
            ALLOCATE(pop%members_to_update(size(pop%candidate_update)))
            pop%members_to_update = pop%members_to_update * 0
            CALL index_of_ones(pop%candidate_update, pop%members_to_update)
            CALL reajust_index(pop%members_to_update)

            CALL apply_beyond_to_members_to_update( exp1, opts1, fitnessfunction, problem1, pop%members_to_update,&
                         common_vars%nfuneval, pop%refset, pop%candidateset,&
                         common_vars%nconst, common_vars%nrand, common_vars%nvar)                 
              
            if ( .not. ALLOCATED(fbest_lastiter) ) then
                ALLOCATE(fbest_lastiter(size(fbest)))
            end if
            fbest_lastiter=fbest
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results, opts1, pop%refset, xbest, fbest, local_solver_var%use_bestx, &
                                     common_vars%nfuneval, cputime1, common_vars%fin, exp1)    

            common_vars%iter=common_vars%iter+1
            local_solver_var%n_minimo=local_solver_var%n_minimo+1 
         
            if (common_vars%idp .NE. 0) then
                CALL printiterationcess(exp1, fbest(1), common_vars%nfuneval, cputime1, 0, 0) 
            end if
            
            CALL printiterationcesslog( exp1, fbest(1), common_vars%nfuneval, cputime1, common_vars%iter)
            CALL update_refset_change(pop%members_update, pop%refset_change)  
            CALL remove_possible_stuck_members_in_refset(problem1, opts1, exp1, fitnessfunction, common_vars%nconst, &
                        pop%refset, pop%refset_change, common_vars%nfuneval, common_vars%nvar, xl_log, xu_log )               
            CALL check_the_improvement_of_fbest(results, opts1, fbest, xbest, fbest_lastiter, &
                        common_vars%nreset, cputime1, common_vars%fin, common_vars%nfuneval)
!----------------------------------------------------------------------------------------------------------------------------------
!       METAHEURISTIC CODE. (5) LOCAL SOLVER
!---------------------------------------------------------------------------------------------------------------------------------- 
            CALL local_solver_method1(exp1,opts1,common_vars,local_solver_var,time,pop,&
                                   problem1,results,fitnessfunction,xbest,fbest,fbest_lastiter) 
            fbest_lastiter=fbest
!----------------------------------------------------------------------------------------------------------------------------------
! CHECK PARALLEL STOPPING CRITERIA SUBROUTINE           
!----------------------------------------------------------------------------------------------------------------------------------              
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results, opts1, pop%refset, xbest, fbest,local_solver_var%use_bestx,&
                        common_vars%nfuneval,cputime1,common_vars%fin, exp1)
            CALL printiterationcess(exp1, fbest(1), common_vars%nfuneval, cputime1, 0, 0)
            CALL stoppingcriteriamodule(exp1,common_vars,opts1,results,fbest,xbest,time)
!----------------------------------------------------------------------------------------------------------------------------------                         
!       METAHEURISTIC CODE. (6) FINALIZE ALGORITHM (7) DEALLOCATE VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
            CALL finalize_algorithm(exp1,problem1,results,opts1,fitnessfunction,pop%refset,xbest,fbest,&
                        common_vars%fin,common_vars%nfuneval,common_vars%NPROC, &
                        cputime1, common_vars%nconst, common_vars%nvar,local_solver_var,common_vars%idp)            

            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL printiterationcess(exp1, fbest(1), common_vars%nfuneval, cputime1, 0, 0) 

            ! SUBROUTINE END ITERATION OUTPUT
            if (common_vars%fin .GT. 0) then
                CALL printiterationcess(exp1, fbest(1), common_vars%nfuneval, cputime1, 0, 0)
                cputime1 = calctimeMPI(exp1, time%starttime)
                CALL printgant(exp1, cputime1, 1)                       
            end if

            CALL destroy_refsettype(pop%candidateset)
            CALL destroy_refsettype(pop%childset)
            if (ALLOCATED(pop%parents_index1)) DEALLOCATE(pop%parents_index1)
            if (ALLOCATED(pop%parents_index2)) DEALLOCATE(pop%parents_index2)
            if (ALLOCATED(pop%candidate_update)) DEALLOCATE(pop%candidate_update)
            if (ALLOCATED(pop%members_update)) DEALLOCATE(pop%members_update)
            if (ALLOCATED(new_comb))  DEALLOCATE(new_comb)
            if (ALLOCATED(pop%members_to_update))  DEALLOCATE(pop%members_to_update)
            if (opts1%useroptions%iterprint .eq. 1 ) then
                cputime1 = calctimeMPI(exp1,time%starttime)
            endif
            IF (ALLOCATED(indvect) ) DEALLOCATE(indvect)
            CALL printrefset(exp1, pop%refset%fpen, opts1%globaloptions%dim_refset, cputime1, &
                                                                        common_vars%iter, pop%refset_change ) 
            CALL printcooperative( exp1, pop%refset%cooperative, opts1%globaloptions%dim_refset, cputime1,  common_vars%iter ) 
            CALL printenditeration(exp1)

          end if
          
        end do
!---------------------------------------------------------------------------------------------------------
! END METAHEURISTIC
!---------------------------------------------------------------------------------------------------------

        CALL mpibarrieress()
        CALL ending_solver_message(common_vars%idp, common_vars%fin)
        CALL mpibarrieress()    
        
!---------------------------------------------------------------------------------------------------------
! GATHER THE RESULTS IN MPI EXECUTION CASE
!---------------------------------------------------------------------------------------------------------
        cputime1 = calctimeMPI(exp1,time%starttime)
        CALL select_better_solution(results, opts1, pop%refset, xbest, fbest, local_solver_var%use_bestx, &
                        common_vars%nfuneval, cputime1, common_vars%fin, exp1)  
        CALL returnminlocelement(exp1, fbest(1), fbest_new, i, common_vars%idp)
        CALL cooperativebcastelement(exp1,xbest,common_vars%nvar,i)
        CALL returnsumelementlong(exp1, common_vars%nfuneval, nfunevaltotal)
        CALL returnavgelementintdouble(exp1, common_vars%iter, results%totaliter, common_vars%NPROC-1)
        CALL returnminelement(exp1, results%timevtr, results%totalvtr)
        cputime1 = calctimeMPI(exp1,time%starttime)
        outresult = 1
        results%timetotal = REAL(cputime1,KIND=C_DOUBLE)
        fbest=fbest_new
        CALL detranslationinterface( xbest, exp1 )
        CALL printmasterocurrencesend( exp1, migration_master%ocurrences_send, common_vars%NPROC, 1, &
                        opts1%globaloptions%dim_refset, &
                        opts1%localoptions%n2,opts1%localoptions%balance)

!---------------------------------------------------------------------------------------------------------
        CALL mpibarrieress() 
        CALL seed_recount(exp1, common_vars%idp)
        CALL mpibarrieress()   
                
        CALL updateresultsess(exp1, results1, results%timetotal, nfunevaltotal, fbest(1), xbest, results%totaliter,results%totalvtr)
        
        CALL settotaltime(results1, results%timetotal)
        CALL setparalleltime(results1, time%timeparallel)
        CALL setlocalsolvertime(results1, time%localsolvertime)
        

        if (common_vars%idp .EQ. 0) then
        !        CALL printsolution(exp1, xbest, fbest)
                CALL  printresults(exp1, fbest, nfunevaltotal,results%timetotal, results%totalvtr,results%totaliter)
        end if
        
        
        ! Free memory
        CALL destroyasynchinitmasterandwindows(exp1)
        CALL destroy_opts(opts1)
        CALL destroy_opts(default1)
        CALL destroy_problem(problem1)

        
        if (ALLOCATED(xl_log)) DEALLOCATE(xl_log)
        if (ALLOCATED(xu_log)) DEALLOCATE(xu_log)   
        if (ALLOCATED(common_vars%BKS_x)) DEALLOCATE(common_vars%BKS_x)
        if (ALLOCATED(common_vars%BKS_f)) DEALLOCATE(common_vars%BKS_f)
        if (ALLOCATED(randomV)) DEALLOCATE(randomV)        
        if (ALLOCATED(pop%solutions)) DEALLOCATE(pop%solutions)

        CALL deallocate_scheme_ess(ess_comm)       
  
        if (ALLOCATED(indvect)) DEALLOCATE(indvect)
        if (ALLOCATED(factor)) DEALLOCATE(factor)
        if (ALLOCATED(pop%refset_change)) DEALLOCATE(pop%refset_change)
        if (ALLOCATED(opts1%localoptions%finish)) DEALLOCATE(opts1%localoptions%finish)
        if (ALLOCATED(xbest)) DEALLOCATE(xbest)
        if (ALLOCATED(fbest)) DEALLOCATE(fbest)
        if (ALLOCATED(fbest_lastiter)) DEALLOCATE(fbest_lastiter)
        CALL destroy_local_solver_vars(local_solver_var)     
        CALL destroy_migration_master(migration_master) 
        CALL mpibarrieress()
        CALL destroy_refsettype(pop%solutionset)
        CALL mpibarrieress()
        CALL destroy_refsettype(pop%refset)
        CALL mpibarrieress()
        CALL destroy_resultsf(results)

        CALL mpibarrieress()

        
    END FUNCTION sacess
#endif

    
END MODULE modsacess
