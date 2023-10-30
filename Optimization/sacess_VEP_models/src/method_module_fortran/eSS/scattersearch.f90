MODULE scattersearch
    USE iso_c_binding
    USE scattersearchtypes
    USE scattersearchfunctions
    USE localsolver
     
CONTAINS

    FUNCTION sscattersearch(exp1, fitnessfunction, results1, maxfunevals, ftarget) RESULT (outresult)
    
        ! Declaración de variables
        USE common_functions
        USE qsort_module   
        IMPLICIT NONE
        TYPE(time_ess) :: time
        TYPE(algorithm_common_vars) :: common_vars
        TYPE(local_solver_help_vars) :: local_solver_var

        TYPE(C_PTR), INTENT(INOUT) :: exp1, results1
        INTEGER (C_LONG), INTENT(INOUT) :: maxfunevals
        REAL (C_DOUBLE), INTENT(INOUT) :: ftarget
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj

        INTEGER :: stopoptimization, fin, nreset
        INTEGER(C_LONG) :: nfuneval
        INTEGER :: i, j, nconst
        REAL(C_DOUBLE) :: cputime1
        INTEGER :: clock_rate, clock_start, clock_stop
        TYPE(opts) :: opts1, default1
        TYPE(problem) :: problem1
        TYPE(resultsf) :: results
        
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::  F02
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xbest, fbest, fbest_lastiter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xl_log, xu_log, randomV

        
        INTEGER, DIMENSION(:), ALLOCATABLE :: ncomb1, refset_change
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: ncomb2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: MaxSubSet, MaxSubSet2
        
        TYPE(Refsettype) :: solutionset, refset, candidateset, childset
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: solutions
        INTEGER :: counter
      
        INTEGER, DIMENSION(:), ALLOCATABLE :: indvect
        INTEGER :: nrand
        INTEGER, DIMENSION(:), ALLOCATABLE :: index1, index2, index, diff_index
        INTEGER :: lb_p, auxsize
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: seed
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: st_p
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: ppp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: hyper_x_L, hyper_x_U, &
                                                                                                factor, v1, v2, v3, new_comb
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: parents_index1, parents_index2
        INTEGER, DIMENSION(:), ALLOCATABLE :: members_update, candidate_update, members_to_update
        INTEGER(C_INT) :: outresult
        REAL(C_DOUBLE) :: resultado
        
        CALL problem_specifications(exp1, problem1,opts1,common_vars%nvar,ftarget,maxfunevals) 

        ! INIT VARS
        time%localsolvertime = 0.0
        time%timeparallel = 0.0
        common_vars%idp = 0
        common_vars%parallel = 0
        stopoptimization = 0        
        nfuneval = 0
        fin = 0
        nreset = 0
        common_vars%NPROC = 1
        
        ! INIT THE TIME VARS
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        CALL SYSTEM_CLOCK(count=clock_start)
        ! CHECKBOUNDS: Check if bounds have the same dimension   
        CALL checkbounds(problem1,common_vars%status)
        !INIT BEST VARIABLES
        CALL initbestvars(problem1,xbest,fbest,common_vars%nvar)
        !INIT INT BIN VARIABLES
        CALL initintbinvars(problem1)
        ! INIT RANDOM GENERATOR
        CALL initrngrandomserial(exp1)
        ! CALC THE REFSET SIZE
        CALL calc_dim_refset(opts1%globaloptions%dim_refset, common_vars%nvar, opts1%useroptions%iterprint,&
        common_vars%idp,common_vars%NPROC) 
        ! CALC THE NDIVERSE SIZE        
        CALL calc_ndiverse(opts1%globaloptions%ndiverse, common_vars%nvar, opts1%useroptions%iterprint,common_vars%idp)
        ! INITIALIZATION OF LOCAL SOLVER
        CALL initlocalsolvervarsess(exp1)
        ! OUTPUT VARS INITIALIZATION
        CALL initoutputvars(exp1)
        ! ENABLE DISTRIBUTED CRITERIA
        CALL setdistcriteria(exp1,0)
        
        if (opts1%globaloptions%ndiverse < opts1%globaloptions%dim_refset) then
            opts1%globaloptions%ndiverse = opts1%globaloptions%dim_refset
        end if

        ! INIT LOCAL OPTIONS
        CALL initlocaloptions(opts1) 
        ! INIT INEQUALITY CONSTRAINTS
        CALL init_inequality_constraints(problem1%neq, problem1%CL, problem1%CU)
        problem1%ineq = problem1%neq + problem1%ineq
        ! CALC CONST NUMBER
        CALL calcnconst(problem1,nconst) 
        ! CHECK OUTPUT OBJECTIVE FUNCTION
        CALL check_output_obj_funct(problem1, opts1,common_vars%status, nconst)

        ! PASAMOS A LOGARITMO
        ALLOCATE(xl_log(size(problem1%XL)))
        ALLOCATE(xu_log(size(problem1%XU)))
        xl_log = problem1%XL
        xu_log = problem1%XU
             
        if (ALLOCATED(opts1%useroptions%log_var)) then    
            CALL converttolog2(xl_log, common_vars%nvar,opts1%useroptions%log_var)
            CALL converttolog2(xu_log, common_vars%nvar,opts1%useroptions%log_var)
        end if

        
        CALL build_aux_local_search(problem1,local_solver_var, opts1, exp1)

        !Initialize variables
        local_solver_var%stage_1=1
        local_solver_var%n_minimo = 0;
        if (opts1%localoptions%empty .eq. 0) local_solver_var%n_critico = opts1%localoptions%n1
        ! Possible combinations among the opts1%globaloptions%dim_refset elements in Refset, taken in pairs
        ALLOCATE(ncomb1(opts1%globaloptions%dim_refset))
        i=1
        ncomb1 = (/(i, i = 1, opts1%globaloptions%dim_refset)/)
        call nchoosek_fpar(ncomb1, 2, ncomb2)
        

        MaxSubSet = (opts1%globaloptions%dim_refset**2 - opts1%globaloptions%dim_refset)/2
        MaxSubSet2 = 2 * MaxSubSet

        ! Creas inicialmente o conxunto soluciónW
        CALL create_init_solutions(exp1,opts1,solutions,common_vars%nvar,xl_log,xu_log);

        
        if ((problem1%int_var .GT. 0) .OR. (problem1%bin_var .GT. 1)) then
            CALL ssm_round_int(solutions, problem1%int_var + problem1%bin_var, problem1%XL, problem1%XU)
        end if

        ! avaliamos o conxunto solución e creamos a estructura solutionset
        CALL evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1, solutions, common_vars%nvar,&
        nfuneval, nconst, xl_log, xu_log)
        
        
        
        
        ! Creamos o conxunto Refset
        CALL create_refset ( exp1, solutionset, refset, common_vars%nvar, opts1, nconst, opts1%globaloptions%dim_refset )


        !Check best value in refset
        xbest=refset%x(:,1)
        fbest=refset%fpen(1)
        CALL check_best_value(refset, opts1, xbest, fbest )
        
        common_vars%iter = 0
        if (opts1%useroptions%iterprint .eq. 1) then
            CALL SYSTEM_CLOCK(count=clock_stop)
            cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)
            WRITE(*, '(A25, I10, A6, D40.30, A10, E15.4, A8, E15.4)') "Initial Pop: NFunEvals: ", nfuneval, &
            " Bestf: ", fbest(1), &
            " CPUTime: ", cputime1, " fbest: ", variance_vector(fbest)
            
        end if
        

        ! Creatos la estructura results y la inicializamos con los primeros valores
        CALL create_results(results, opts1, xbest, fbest, nfuneval, cputime1, common_vars%nvar)

        ! Sort the Refset
        CALL sort_refset(refset,opts1, indvect)

        if (opts1%globaloptions%combination .eq. 1) then         
            nrand = common_vars%nvar
        else if (opts1%globaloptions%combination .eq. 2) then
            nrand = 1
        end if

        local_solver_var%use_bestx = 0

        auxsize = choose(ncomb1)
        ALLOCATE(index1(auxsize))
        ALLOCATE(index2(auxsize))
        ALLOCATE(index(auxsize))
        ALLOCATE(diff_index(auxsize))
        index1 = ncomb2(:, 1)
        index2 = ncomb2(:, 2)

        index = ncomb2(:, 2)
        
        
        CALL FUSION_VECTOR_INT(index1, index)

        diff_index = (index2 - index1)
        lb_p = 1 !Minimum ppp
        st_p = REAL(0.75d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        !Defines maximum ppp. max(ppp)= lb_p+ st_p
        !Example: lb_p=1.25 and st_p=0.5 varies ppp from
        !%1.25 to 1.75

        ALLOCATE (ppp(auxsize))
        ppp = st_p * REAL((diff_index - 1),KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/ &
        REAL((opts1%globaloptions%dim_refset - 2),KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) + &
        REAL(lb_p,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        

        !hyper_x_L=repmat(problem1%XL,MaxSubSet,1);
        !hyper_problem1%XU=repmat(problem1%XU,MaxSubSet,1);
        ALLOCATE(hyper_x_L(common_vars%nvar,CEILING(MaxSubSet)))
        ALLOCATE(hyper_x_U(common_vars%nvar,CEILING(MaxSubSet)))
        i=1
        j=1
        hyper_x_L = reshape((/ ((problem1%XL(i), i = 1,common_vars%nvar), j = 1,  CEILING(MaxSubSet)) /), &
        (/ common_vars%nvar, CEILING(MaxSubSet) /))
        hyper_x_U = reshape((/ ((problem1%XU(i), i = 1,common_vars%nvar), j = 1,  CEILING(MaxSubSet)) /), &
        (/ common_vars%nvar, CEILING(MaxSubSet) /))


        ALLOCATE(refset_change(opts1%globaloptions%dim_refset))
        refset_change = 0

        
        F02 = fbest(1)
        if (opts1%useroptions%iterprint .eq. 1 ) then
            CALL SYSTEM_CLOCK(count=clock_stop)
            cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)
            print *, "CALL initprintfile(exp1, F02, common_vars%parallel, common_vars%idp, cputime1, nfuneval,0)"
            CALL initprintfile(exp1, F02, common_vars%parallel, common_vars%idp, cputime1, nfuneval,0)                
        endif           
        
        ! Algorithm main loop
        do while (fin .eq. 0)
            
            CALL create_childset ( childset, common_vars%nvar, CEILING(MaxSubSet2), nconst )   
            
            counter = 1
            if (.not.ALLOCATED(members_update)) then
                ALLOCATE(members_update(opts1%globaloptions%dim_refset))
                members_update =  1
            end if
            
    
            CALL check_duplicated_replace (exp1, fitnessfunction, problem1,refset,opts1, parents_index1, parents_index2, &
                        xl_log, xu_log, refset_change, members_update, index1, index2, nfuneval,&
                        common_vars%nvar, nconst )
        
            CALL create_candidate(candidateset, refset)

            if (.not.ALLOCATED(candidate_update)) then
                ALLOCATE(candidate_update(opts1%globaloptions%dim_refset))
                candidate_update = 0
            end if

            ALLOCATE( factor(common_vars%nvar, size(ppp)))
            factor = reshape((/ ((ppp(j), i = 1, common_vars%nvar), j = 1, size(ppp)) /), (/ common_vars%nvar,size(ppp) /))
            factor = factor * (parents_index2 - parents_index1)/ REAL(1.5d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) 

            
            ALLOCATE( v1(common_vars%nvar, size(ppp)))
            ALLOCATE( v2(common_vars%nvar, size(ppp)))
            ALLOCATE( v3(common_vars%nvar, size(ppp)))


                       
            v1 = parents_index1 - factor
            v2 = parents_index2 - factor
            v3 = REAL(2.0d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  * parents_index2  - parents_index1 - factor  
            
 
            CALL check_vector_in_hyper_bounds(exp1,opts1, v1, hyper_x_L, hyper_x_U )
            CALL check_vector_in_hyper_bounds(exp1,opts1, v3, hyper_x_L, hyper_x_U )
            CALL generate_new_comb_matrix(exp1, new_comb, ppp, MaxSubSet, nrand, v1,v2,v3)

                        
 
            CALL update_candidateset_with_new_comb( exp1,opts1, fitnessfunction,problem1,new_comb,candidateset,childset,&
                candidate_update, members_update,MaxSubSet2,nrand,nconst,nfuneval,counter, index )
             
        
!----------------------------------------------------------------------------------------------------------------------------------
! CHECK PARALLEL STOPPING CRITERIA           
!----------------------------------------------------------------------------------------------------------------------------------              
                
            CALL SYSTEM_CLOCK(count=clock_stop)
            cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)     
            CALL select_better_solution(results, opts1, refset, xbest, fbest, local_solver_var%use_bestx,nfuneval, &
                                                        cputime1, fin, exp1)              
            CALL printiteration(exp1, common_vars%iter, fbest(1), nfuneval, cputime1) 
            fin = check_stopping_criteria (problem1, opts1, cputime1, stopOptimization, fbest, nfuneval)                
!----------------------------------------------------------------------------------------------------------------------------------
   
            
            CALL reajustchild(childset, counter, common_vars%nvar,  nconst)
            
            
            ALLOCATE(members_to_update(size(candidate_update)))
            members_to_update = members_to_update * 0
            CALL index_of_ones(candidate_update, members_to_update)
            CALL reajust_index(members_to_update)
       
            CALL apply_beyond_to_members_to_update( exp1, opts1, fitnessfunction, problem1, members_to_update,  nfuneval, &
                    refset, candidateset, nconst, nrand, common_vars%nvar)
            
            if ( .not. ALLOCATED(fbest_lastiter) ) then
                ALLOCATE(fbest_lastiter(size(fbest)))
            end if
            fbest_lastiter=fbest
            cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)                 
            CALL select_better_solution(results, opts1, refset, xbest, fbest, local_solver_var%use_bestx, nfuneval,cputime1, fin,&
            exp1) 

            common_vars%iter=common_vars%iter+1
            local_solver_var%n_minimo=local_solver_var%n_minimo+1
            call cpu_time(cputime1) 

            if (opts1%useroptions%iterprint .eq. 1) then
                WRITE(*, '(A13, I10, A6, I10, A9, D40.30, A12, f9.2, A12, E15.4,A5, E15.4)') "Iterations: ", &
                common_vars%iter, " NFunEvals: ",&
                nfuneval, & 
                " Bestf: ", fbest(1), &
                " CPUTime: ", cputime1, &
                " variance: ", variance_vector(refset%fpen), &
                " VTR", problem1%vtr
            end if
            
            CALL update_refset_change(members_update, refset_change)  
            
   
            CALL remove_possible_stuck_members_in_refset(problem1, opts1, exp1, fitnessfunction, nconst, &
                        refset, refset_change,nfuneval, common_vars%nvar, xl_log, xu_log )

            
            CALL check_the_improvement_of_fbest(results, opts1, fbest, xbest, fbest_lastiter, nreset, cputime1, fin, nfuneval)

            if (( opts1%localoptions%empty .eq. 0 )  .and. (fin .eq. 0)) then
                if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. (opts1%localoptions%bestx .gt. 0) .and. &
                                    (local_solver_var%n_minimo .ge. local_solver_var%n_critico) ) then    
                 if ( local_solver_var%use_bestx .gt. 0 ) then
                    CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                                childset,xbest,fbest,nfuneval,fbest_lastiter,refset,refset_change,0)                    
                    local_solver_var%n_minimo=0
                 end if
                else 
                 if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. (local_solver_var%n_minimo .ge.&
                 local_solver_var%n_critico) ) then
                    CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                                childset,xbest,fbest,nfuneval,fbest_lastiter,refset,refset_change,0)   
                    local_solver_var%n_minimo=0
                 end if 
                end if
            end if
            
            fbest_lastiter=fbest
            
            ! check the stopping criteria    
            CALL SYSTEM_CLOCK(count=clock_stop)
            cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)  
            fin = check_stopping_criteria (problem1, opts1, cputime1, stopoptimization, fbest, nfuneval)
            
            CALL finalize_algorithm(exp1,problem1,results,opts1,fitnessfunction,refset,xbest,fbest,fin,nfuneval,&
            common_vars%NPROC, &
                    cputime1, nconst, common_vars%nvar, local_solver_var,common_vars%idp)
            
            
            if (opts1%useroptions%iterprint .eq. 1 ) then
                CALL SYSTEM_CLOCK(count=clock_stop)
                cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)
                CALL printiterationcesslog(exp1, fbest(1), nfuneval, cputime1, common_vars%iter)
                CALL printiteration(exp1, common_vars%iter, fbest(1), nfuneval, cputime1) 
            endif            
           
            CALL addlocalscounteress(exp1)
            if (fin > 0) then
                if (opts1%useroptions%iterprint .eq. 1 ) then
                    CALL SYSTEM_CLOCK(count=clock_stop)
                    cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)
                endif
                CALL printiteration(exp1, common_vars%iter, fbest(1), nfuneval,cputime1)
            end if
            
            CALL destroy_refsettype(candidateset)
            CALL destroy_refsettype(childset)
            
            
            if (ALLOCATED(parents_index1)) DEALLOCATE(parents_index1)
            if (ALLOCATED(parents_index2)) DEALLOCATE(parents_index2)

            if (ALLOCATED(factor))  DEALLOCATE(factor)
            if (ALLOCATED(v1))  DEALLOCATE(v1)
            if (ALLOCATED(v2))  DEALLOCATE(v2)
            if (ALLOCATED(v3))  DEALLOCATE(v3)
            if (ALLOCATED(candidate_update)) DEALLOCATE(candidate_update)
            if (ALLOCATED(members_update)) DEALLOCATE(members_update)
            if (ALLOCATED(new_comb))  DEALLOCATE(new_comb)
            if (ALLOCATED(members_to_update))  DEALLOCATE(members_to_update)
            
        end do
        
        ! SUBROUTINE DESTROY_PROBLEMS
        outresult = 1
        CALL SYSTEM_CLOCK(count=clock_stop)
        cputime1 = REAL(clock_stop-clock_start)/REAL(clock_rate)  
        results%timetotal = REAL(cputime1,KIND=C_DOUBLE)
        resultado = fbest(1)
        
        
        CALL seed_recount(exp1, common_vars%idp)

        
        CALL updateresultsess(exp1, results1, results%timetotal, nfuneval, fbest(1), xbest, common_vars%iter, results%timetotal)
        CALL settotaltime(results1, results%timetotal)
        CALL setlocalsolvertime(results1, time%localsolvertime)

        ! Free memory
        CALL destroy_opts(opts1)
        CALL destroy_opts(default1)
        CALL destroy_problem(problem1)

        
        if (ALLOCATED(xl_log)) DEALLOCATE(xl_log)
        if (ALLOCATED(xu_log)) DEALLOCATE(xu_log)   
        if (ALLOCATED(randomV)) DEALLOCATE(randomV)
        if (ALLOCATED(solutions)) DEALLOCATE(solutions)
        if (ALLOCATED(ncomb1)) DEALLOCATE(ncomb1)
        if (ALLOCATED(ncomb2)) DEALLOCATE(ncomb2)
        if (ALLOCATED(indvect)) DEALLOCATE(indvect)
        if (ALLOCATED(index1)) DEALLOCATE(index1)
        if (ALLOCATED(index2)) DEALLOCATE(index2)
        if (ALLOCATED(index)) DEALLOCATE(index)
        if (ALLOCATED(diff_index)) DEALLOCATE(diff_index)
        if (ALLOCATED(diff_index)) DEALLOCATE(diff_index)
        if (ALLOCATED(factor)) DEALLOCATE(factor)
        if (ALLOCATED(ppp)) DEALLOCATE(ppp)
        if (ALLOCATED(hyper_x_L)) DEALLOCATE(hyper_x_L)
        if (ALLOCATED(hyper_x_U)) DEALLOCATE(hyper_x_U)
        if (ALLOCATED(refset_change)) DEALLOCATE(refset_change)
        if (ALLOCATED(opts1%localoptions%finish)) DEALLOCATE(opts1%localoptions%finish)
        if (ALLOCATED(xbest)) DEALLOCATE(xbest)
        if (ALLOCATED(fbest)) DEALLOCATE(fbest)
        if (ALLOCATED(fbest_lastiter)) DEALLOCATE(fbest_lastiter)     
        
        CALL destroy_local_solver_vars(local_solver_var)     
        CALL destroy_refsettype(solutionset)
        CALL destroy_refsettype(refset)
        CALL destroy_resultsf(results)
        
    END FUNCTION sscattersearch

END MODULE scattersearch
