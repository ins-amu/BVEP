MODULE modcess
    USE iso_c_binding
    USE scattersearchtypes
    USE scattersearchfunctions
    USE localsolver
    USE parallelscattersearchfunctions
#ifdef MPI2
   
        
CONTAINS
    
    FUNCTION cess(exp1, fitnessfunction, results1, maxfunevals, ftarget) RESULT (outresult)
    
        ! Declaración de variables
        USE common_functions
        USE qsort_module   
        IMPLICIT NONE
        TYPE(local_solver_help_vars) :: local_solver_var
        TYPE(time_ess) :: time
        TYPE(algorithm_common_vars) :: common_vars
        TYPE(C_PTR), INTENT(INOUT) :: exp1, results1
        REAL (C_DOUBLE), INTENT(INOUT) :: ftarget
        INTEGER (C_LONG), INTENT(INOUT) :: maxfunevals
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj

        INTEGER(C_INT) :: returnseedcounter, counterseed
        INTEGER :: stopoptimization, fin, stopfin, nreset,  l, rest
        INTEGER(C_LONG) :: nfuneval, checkcooperativemigrationcriteriacess, nfunevaltotal
        INTEGER :: i, j, nconst
        REAL(C_DOUBLE) :: cputime1
        INTEGER :: clock_start
        TYPE(opts) :: opts1, default1
        TYPE(problem) :: problem1
        TYPE(resultsf) :: results
        
        INTEGER, DIMENSION(:), ALLOCATABLE :: dim_refset_array
        INTEGER :: concat_dim_refset, jj
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: F02
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: seed
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xbest, fbest, fbest_lastiter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: xl_log, xu_log, randomV, fbestg
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: xbestg

        INTEGER, DIMENSION(:), ALLOCATABLE :: ncomb1, refset_change
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: ncomb2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: MaxSubSet, MaxSubSet2
        
        TYPE(Refsettype) :: solutionset, refset, candidateset, childset, solutionset_global
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: solutions
        INTEGER :: counter
      
        INTEGER, DIMENSION(:), ALLOCATABLE :: indvect
        INTEGER :: nrand
        INTEGER, DIMENSION(:), ALLOCATABLE :: index1, index2, index, diff_index
        INTEGER :: lb_p, auxsize
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: st_p
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: ppp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: hyper_x_L, hyper_x_U, &
                                                                                                factor, v1, v2, v3, new_comb
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: parents_index1, parents_index2
        INTEGER, DIMENSION(:), ALLOCATABLE :: members_update, candidate_update, members_to_update
        INTEGER(C_INT) :: outresult        
        
        
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::  NAN, INF
        INTEGER(C_INT) :: error1, setnan
        INTEGER :: tam
        INTEGER(C_INT) :: getopenmpoption, openmp_pos

        !COOPERATIVE VARIABLES
        REAL(C_DOUBLE) :: initmpi, calctimempi
        INTEGER(C_LONG) :: mig
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: servector, servector_global

!-------------------------------------------------------------------------------------------------------!
! SETTINGS AND INICIALIZATION 
!-------------------------------------------------------------------------------------------------------!
        error1 = setnan(NAN)
        time%timeparallel = 0.0
        time%localsolvertime =  0.0
        CALL setdblmax(INF)
        CALL problem_specifications(exp1, problem1,opts1,common_vars%nvar,ftarget,maxfunevals) 
        

       
! Inicializanse e seteanse algunhas variables 
        stopoptimization = 0        
        nfuneval = 0
        fin = 0
        stopfin = 0
        nreset = 0
        CALL setNPROC(exp1, common_vars%NPROC)            
        ! Check if bounds have the same dimension   CHECKBOUNDS
        CALL checkbounds(problem1,common_vars%status)
        !INIT BEST VARIABLES
        CALL initbestvars(problem1,xbest,fbest,common_vars%nvar)
        !INIT INT BIN VARIABLES
        CALL initintbinvars(problem1)
        CALL chargeid(exp1,common_vars%idp)
        CALL initrngrandomparallel(exp1,common_vars%idp)
        ! DISABLE DISTRIBUTED CRITERIA
        CALL setdistcriteria(exp1,0)
        
        
! Calculase o tamanho do refset
        CALL calc_dim_refset(opts1%globaloptions%dim_refset, common_vars%nvar, opts1%useroptions%iterprint,&
        common_vars%idp,common_vars%NPROC) 
! Asignase a cada procesador unha configuracion diferente
        CALL hete_param_eSS2(exp1, problem1, opts1, common_vars%NPROC, opts1%globaloptions%dim_refset, &
        common_vars%idp, common_vars%nvar)
! Calculase o tamanho do conxunto de solucions iniciais a partir do cal se creara o refset 
        CALL calc_ndiverse(opts1%globaloptions%ndiverse, common_vars%nvar, opts1%useroptions%iterprint,common_vars%idp)
        CALL initlocalsolvervarsess(exp1)
        CALL initoutputvars(exp1)
        
        if (opts1%globaloptions%ndiverse < opts1%globaloptions%dim_refset) then
            opts1%globaloptions%ndiverse = opts1%globaloptions%dim_refset
        end if
        
        
        common_vars%parallel = 1
        CALL slave_information(common_vars%idp, opts1%globaloptions%dim_refset, opts1%globaloptions%ndiverse, 0)        
        time%starttime = initmpi()
        CALL saveinittime(exp1,time%starttime)
        CALL chargecooperativeparametersfortran(exp1, opts1%globaloptions%dim_refset, tam, common_vars%idp, ftarget, maxfunevals)  
        CALL createcooperativetopologyess(exp1)

        ! INIT LOCAL OPTIONS
        CALL initlocaloptions(opts1) 
        CALL init_inequality_constraints(problem1%neq, problem1%CL, problem1%CU)
   
        problem1%ineq = problem1%neq + problem1%ineq
       
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
        

!-------------------------------------------------------------------------------------------------------!        
! CALC INIT SOLUTION SET AND REFSET GLOBAL       
!-------------------------------------------------------------------------------------------------------!  
        
        
        CALL create_init_solutions(exp1,opts1,solutions,common_vars%nvar,xl_log,xu_log);     
        

        if ((problem1%int_var .GT. 0) .OR. (problem1%bin_var .GT. 1)) then
            CALL ssm_round_int(solutions, problem1%int_var + problem1%bin_var, problem1%XL, problem1%XU)
        end if
#ifdef OPENMP            
        openmp_pos = getopenmpoption(exp1)
        ! avaliamos o conxunto solución e creamos a estructura solutionset
        if ( openmp_pos .EQ. 1) then
        CALL evaluate_solutions_set_parallel(exp1,fitnessfunction,solutionset,problem1,opts1, solutions,&
        common_vars%nvar, nfuneval, nconst,xl_log, xu_log )
        else 
        CALL evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1, solutions, common_vars%nvar, &
                nfuneval, nconst,xl_log, xu_log)
        end if
#else        
        CALL evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1, solutions, common_vars%nvar, &
                nfuneval, nconst,xl_log, xu_log)
#endif
          
        ! Creamos o conxunto Refset   
        CALL create_refset (exp1, solutionset, refset, common_vars%nvar, opts1, nconst, opts1%globaloptions%dim_refset )

       
       
!-------------------------------------------------------------------------------------------------------!         
   
        !Initialize variables
        local_solver_var%stage_1=1
        local_solver_var%n_minimo = 0;
        if (opts1%localoptions%empty .eq. 0) local_solver_var%n_critico = opts1%localoptions%n1

        ! Possible combinations among the opts1%globaloptions%dim_refset elements in Refset, taken in pairs
        ALLOCATE(ncomb1(opts1%globaloptions%dim_refset))
        ncomb1 = (/(i, i = 1, opts1%globaloptions%dim_refset)/)
        call nchoosek_fpar(ncomb1, 2, ncomb2)
        

        MaxSubSet = (opts1%globaloptions%dim_refset**2 - opts1%globaloptions%dim_refset)/2
        MaxSubSet2 = 2 * MaxSubSet        

        !Check best value in refset
        xbest=refset%x(:,1)
        fbest=refset%fpen(1)
        CALL check_best_value(refset, opts1, xbest, fbest )
        common_vars%iter = 0

!-------------------------------------------------------------------------------------------------------!        
! CALC INIT SOLUTION SET AND REFSET GLOBAL       
!-------------------------------------------------------------------------------------------------------!         
        CALL sizeseriallize(common_vars%nvar, nconst, common_vars%sizeser)
        CALL optimization_begins(common_vars%idp)    
        
        if ((opts1%useroptions%iterprint .eq. 1) .AND. (common_vars%idp .NE. 0)) then
                cputime1 = calctimeMPI(exp1,time%starttime)
                WRITE(*, '(A14, I3, A26, I10, A8, D25.10, A10, F12.2, A11, E10.1)') "[PROCESSOR ID=", &
                    common_vars%idp, "] Initial Pop: NFunEvals: ", nfuneval, &
                    " Bestf: ", fbest(1), " CPUTime: ", cputime1, " variance: ", variance_vector(refset%fpen)
        endif    
        

        ! Creatos la estructura results y la inicializamos con los primeros valores
        CALL create_results(results, opts1, xbest, fbest, nfuneval, cputime1, common_vars%nvar)
        

!-------------------------------------------------------------------------------------------------------!
        
        ! Sort the Refset
        CALL sort_refset(refset,opts1, indvect)

        if (common_vars%idp .EQ. 0 ) then
            refset%x = NAN
            refset%fpen = NAN
            refset%f = NAN
            refset%penalty = NAN
            if (nconst .GT. 0) then
                refset%nlc = NAN         
            end if            
        end if

                
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
        hyper_x_L = reshape((/ ((problem1%XL(i), i = 1,common_vars%nvar), j = 1,  CEILING(MaxSubSet)) /),&
        (/ common_vars%nvar, CEILING(MaxSubSet) /))
        hyper_x_U = reshape((/ ((problem1%XU(i), i = 1,common_vars%nvar), j = 1,  CEILING(MaxSubSet)) /),&
        (/ common_vars%nvar, CEILING(MaxSubSet) /))


        ALLOCATE(refset_change(opts1%globaloptions%dim_refset))
        refset_change = 0
        
!-------------------------------------------------------------------------------------------------------!        
        cputime1 = calctimeMPI(exp1,time%starttime)
        F02 = fbest(1)
        CALL initprintfile(exp1, F02, common_vars%parallel, common_vars%idp, cputime1, nfuneval,0)           
        !CALL printgant(exp1,cputime1,1) 

        
!-------------------------------------------------------------------------------------------------------!         
! INIT : Algorithm main loop
!-------------------------------------------------------------------------------------------------------! 
        do while (stopfin .eq. 0)
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL printinititeration(exp1, common_vars%iter, cputime1)

            !CALL printgant(exp1,cputime1) 
                
            if (common_vars%idp .NE. 0) then 
                CALL printrefset(exp1, refset%fpen, opts1%globaloptions%dim_refset,cputime1, common_vars%iter, refset_change )
                F02 = fbest(1)
                CALL printcheckvtr(exp1, F02, nfuneval, cputime1, results%timevtr ) 
            end if
!----------------------------------------------------------------------------------------------------------------------------------
! MIGRATION REGION            
!----------------------------------------------------------------------------------------------------------------------------------   
        CALL printinititeration(exp1, common_vars%iter, cputime1)
        mig =  checkcooperativemigrationcriteriacess(exp1, cputime1);

        if (mig .EQ. 1) THEN
                 
            
            cputime1 = calctimeMPI(exp1,time%starttime)
            if (common_vars%idp .NE. 0) then 
                    fbest_lastiter=fbest     
                    F02=fbest(1)
                    CALL printiterationcess(exp1, F02, nfuneval,cputime1, 0, 1) 
              end if
            !Primeiro comprobamos o criterio de parada
            CALL setcountstopsvaressmaster(exp1, fin, stopfin) ! Aqui sincroniza

            cputime1 = calctimeMPI(exp1,time%starttime)
            if ((common_vars%idp .EQ. 0).AND. (opts1%useroptions%iterprint .eq. 1)) then
                print * , ""
                WRITE (*, '(A45,F10.2)') "[MASTER] START MIGRATION OPERATIONS. CPUTime:", cputime1 
            end if
            if (stopfin .NE. 1 ) then
                ALLOCATE(xbestg (common_vars%nvar, common_vars%NPROC ))
                ALLOCATE(fbestg (common_vars%NPROC ) )
                ALLOCATE(dim_refset_array (common_vars%NPROC) )

                !print *, "cooperativegathertelementint"
                CALL cooperativegathertelementint(exp1,opts1%globaloptions%dim_refset, 1, dim_refset_array) 
                CALL cooperativebcastelementint(exp1,dim_refset_array, common_vars%NPROC)
                concat_dim_refset = 0
                do l=2,common_vars%NPROC
                    concat_dim_refset =  concat_dim_refset + dim_refset_array(l)
                end do
                
                !print *, "create_refset_empty"
                CALL create_refset_empty ( solutionset_global, common_vars%nvar, opts1, nconst, concat_dim_refset )
                ALLOCATE(servector( common_vars%sizeser , opts1%globaloptions%dim_refset))  
                ALLOCATE(servector_global( common_vars%sizeser , concat_dim_refset) )
                CALL serializerefset(refset, common_vars%nvar, opts1%globaloptions%dim_refset, servector, rest ) 
                
                !print *, "migrationsynchcooperativesentcess"
                CALL migrationsynchcooperativesentcess(exp1, common_vars%sizeser, servector, servector_global, dim_refset_array, &
                                concat_dim_refset)
 
                CALL deserializerefset(solutionset_global, common_vars%nvar, concat_dim_refset, servector_global )    
                
                CALL cooperativegathertelement(exp1,xbest, common_vars%nvar, xbestg ) 
                CALL cooperativegathertelement(exp1,fbest, 1, fbestg ) 
                 
                CALL cooperativebcastelement(exp1, xbestg, common_vars%nvar*(common_vars%NPROC), 0 )
                CALL cooperativebcastelement(exp1, fbestg, common_vars%NPROC, 0 )
                CALL incrementmaxtimemig(exp1)
                
                if (common_vars%idp .EQ. 0 ) then
                    if ((common_vars%idp .EQ. 0).AND. (opts1%useroptions%iterprint .eq. 1)) then
                        cputime1 = calctimeMPI(exp1,time%starttime)
                        WRITE (*, '(A49)') "[MASTER] BCAST THE BEST SOLUTIONS FOR EACH SLAVE:" 
                        do jj=2,common_vars%NPROC
                            WRITE (*, '(A15, I3, A6, F15.4)') "      SLAVE ID=", jj-1 ," - FX=", fbestg(jj)
                            
                        end do
                    end if
                                
                    DEALLOCATE(servector)
                    DEALLOCATE(servector_global)
                    DEALLOCATE(xbestg)
                    DEALLOCATE(fbestg)
                    DEALLOCATE(dim_refset_array) 
                else
                  if (common_vars%idp .NE. 0) then
                         DEALLOCATE(solutions)
                         CALL create_init_solutions(exp1,opts1,solutions,common_vars%nvar,xl_log,xu_log)
                  end if
                  if (ALLOCATED(problem1%X0)) DEALLOCATE(problem1%X0)
                  if (ALLOCATED(problem1%F0)) DEALLOCATE(problem1%F0)

                  ALLOCATE(problem1%X0(common_vars%nvar, concat_dim_refset + common_vars%NPROC - 1))
                  ALLOCATE(problem1%F0(concat_dim_refset + common_vars%NPROC -1 ))
             
                  problem1%X0(:,1:common_vars%NPROC-1) = xbestg(:,2:common_vars%NPROC)
                  problem1%X0(:,common_vars%NPROC:concat_dim_refset+common_vars%NPROC-1) =  solutionset_global%x
                  problem1%F0(1:common_vars%NPROC-1) = fbestg(2:common_vars%NPROC)
                  problem1%F0(common_vars%NPROC:concat_dim_refset+common_vars%NPROC-1) =  solutionset_global%fpen     

                  if ((problem1%int_var .GT. 0) .OR. (problem1%bin_var .GT. 1)) then
                    CALL ssm_round_int(solutions, problem1%int_var + problem1%bin_var, problem1%XL, problem1%XU)
                  end if
            
                  CALL destroy_refsettype(solutionset)
                  CALL destroy_refsettype(refset)
                
                  CALL unique_solution (problem1, problem1%X0, problem1%F0, common_vars%nvar )
 
                  if (common_vars%idp .NE. 0) then
                  CALL evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1, solutions, common_vars%nvar, &
                              nfuneval, nconst,xl_log, xu_log)
                  end if     

                  CALL create_refset ( exp1, solutionset, refset, common_vars%nvar, opts1, nconst,&
                  opts1%globaloptions%dim_refset )                       
               
                  DEALLOCATE(servector)
                  DEALLOCATE(servector_global)
                  DEALLOCATE(xbestg)
                  DEALLOCATE(fbestg)
                  DEALLOCATE(dim_refset_array)
                
                end if

                CALL destroy_refsettype(solutionset_global)
                
            
                cputime1 = calctimeMPI(exp1,time%starttime)
           
                local_solver_var%stage_1 = 1
                local_solver_var%stage_2 = 0
                refset_change = 0
                local_solver_var%use_bestx = 0

                if (common_vars%idp .NE. 0) then
                      if (ALLOCATED(indvect)) DEALLOCATE(indvect)
                      CALL sort_refset(refset,opts1, indvect)               
                      fbest=refset%fpen(1)
                      xbest=refset%x(:,1)                
                      F02 = fbest(1)
                      CALL printiterationcess(exp1, F02, nfuneval,cputime1, 0, 2) 
                end if 
                
                IF (ALLOCATED(local_solver_var%initial_points)) DEALLOCATE(local_solver_var%initial_points)
                IF (ALLOCATED(local_solver_var%initial_points_values)) DEALLOCATE(local_solver_var%initial_points_values)
                IF (ALLOCATED(local_solver_var%local_solutions)) DEALLOCATE(local_solver_var%local_solutions)
                IF (ALLOCATED(local_solver_var%local_solutions_values)) DEALLOCATE(local_solver_var%local_solutions_values)        
        
            end if
            cputime1 = calctimeMPI(exp1,time%starttime)
            !CALL printgant(exp1,cputime1,3)    
            if ((common_vars%idp .EQ. 0).AND. (opts1%useroptions%iterprint .eq. 1)) then
                WRITE (*, '(A43,F10.2)') "[MASTER] END MIGRATION OPERATIONS. CPUtime:", cputime1             
            end if
        endif
        mig = 0
            
            
            
        if ( (stopfin .NE. 1) .AND. (common_vars%idp .NE. 0 )) then
!----------------------------------------------------------------------------------------------------------------------------------
! SEQUENTIAL ALGORITHM
!----------------------------------------------------------------------------------------------------------------------------------
     
            CALL create_childset ( childset, common_vars%nvar, CEILING(MaxSubSet2), nconst )   
            counter = 1
            if (.not.ALLOCATED(members_update)) then
                ALLOCATE(members_update(opts1%globaloptions%dim_refset))
                members_update =  1
            end if

! ELIMINACION DE DUPLICADOS E DE SOLUCION MOI PARECIDAS QUE ESTAN NO REFSET                        
            CALL check_duplicated_replace (exp1, fitnessfunction, problem1,refset,opts1, parents_index1, parents_index2, &
                        xl_log, xu_log, refset_change, members_update, index1, index2, nfuneval, common_vars%nvar, nconst )     

! CREATION DAS SOLUCION CANDIDATAS 
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

          
! EVALUACION DAS SOLUCION CANDIDATAS 
#ifdef OPENMP
            openmp_pos = getopenmpoption(exp1) 
            CALL update_candidateset_with_new_comb_parallel( exp1,opts1, fitnessfunction,problem1,new_comb,candidateset,&
                childset,&
                candidate_update, members_update,CEILING(MaxSubSet2),nrand,nconst,nfuneval,counter, index)
#else        
            CALL update_candidateset_with_new_comb( exp1,opts1,fitnessfunction,problem1,new_comb,candidateset,childset,&
                candidate_update,members_update,MaxSubSet2,nrand,nconst,nfuneval,counter, index )
#endif                

            CALL reajustchild(childset, counter, common_vars%nvar, nconst)
            
            ALLOCATE(members_to_update(size(candidate_update)))
            members_to_update = members_to_update * 0
            CALL index_of_ones(candidate_update, members_to_update)
            CALL reajust_index(members_to_update)
            
!----------------------------------------------------------------------------------------------------------------------------------
! CHECK PARALLEL STOPPING CRITERIA           
!----------------------------------------------------------------------------------------------------------------------------------                          
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results,opts1,refset,xbest,fbest,local_solver_var%use_bestx,nfuneval,cputime1,fin, exp1) 
            CALL printiterationcess(exp1, fbest(1), nfuneval,cputime1, 0, 0) 
            fin = check_stopping_criteria (problem1, opts1,cputime1, stopOptimization, fbest, nfuneval)

!----------------------------------------------------------------------------------------------------------------------------------
            
! APLICACION DO METODO DE MELLORA 
           CALL apply_beyond_to_members_to_update( exp1, opts1, fitnessfunction, problem1, members_to_update,  nfuneval, &
                  refset, candidateset, nconst, nrand, common_vars%nvar)
                             
            if ( .not. ALLOCATED(fbest_lastiter) ) then
                ALLOCATE(fbest_lastiter(size(fbest)))
            end if
            fbest_lastiter=fbest
            cputime1 = calctimeMPI(exp1,time%starttime)
            CALL select_better_solution(results,opts1,refset,xbest,fbest,local_solver_var%use_bestx,nfuneval,cputime1,fin, exp1) 

            common_vars%iter=common_vars%iter+1
            local_solver_var%n_minimo=local_solver_var%n_minimo+1
                   

            CALL printiterationcesslog( exp1, fbest(1),nfuneval, cputime1, common_vars%iter)

            CALL update_refset_change(members_update, refset_change)  
            
                             
            CALL remove_possible_stuck_members_in_refset(problem1, opts1, exp1, fitnessfunction, nconst, &
                        refset, refset_change,nfuneval, common_vars%nvar, xl_log, xu_log )
                           
            
            CALL check_the_improvement_of_fbest(results, opts1, fbest, xbest, fbest_lastiter, nreset, cputime1, fin, nfuneval)
  
            if ((opts1%localoptions%empty .eq. 0 )  .and. (fin .eq. 0))then

              F02=fbest(1)
              if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. (opts1%localoptions%bestx .gt. 0) .and. &
                       (local_solver_var%n_minimo .ge. local_solver_var%n_critico) ) then    

                  if ( local_solver_var%use_bestx .gt. 0 ) then
                    cputime1 = calctimeMPI(exp1,time%starttime)
                    !CALL printgant(exp1,cputime1,1)
                    CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                                childset,xbest,fbest,nfuneval,fbest_lastiter,refset,refset_change,0)
                    cputime1 = calctimeMPI(exp1,time%starttime)                                
                    !CALL printgant(exp1,cputime1,2)                    
                    local_solver_var%n_minimo=0
                  end if
              else 

                  if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. &
                            (local_solver_var%n_minimo .ge. local_solver_var%n_critico) ) then

                    cputime1 = calctimeMPI(exp1,time%starttime)
                    !CALL printgant(exp1,cputime1,1)                      
                    CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                                childset,xbest,fbest,nfuneval,fbest_lastiter,refset,refset_change,0)
                    cputime1 = calctimeMPI(exp1,time%starttime)                                
                    !CALL printgant(exp1,cputime1,2)           
                    local_solver_var%n_minimo=0
                  end if 
              end if
              F02=fbest(1)
                
            end if


            fbest_lastiter=fbest 
            cputime1 = calctimeMPI(exp1,time%starttime)
                
            fin = check_stopping_criteria (problem1, opts1, cputime1, stopoptimization, fbest, nfuneval)
            CALL finalize_algorithm(exp1,problem1,results,opts1,fitnessfunction,refset,xbest,fbest,fin,nfuneval,common_vars%NPROC, &
                    cputime1, nconst, common_vars%nvar, local_solver_var, common_vars%idp)
            
            F02 = fbest(1)
            if (common_vars%idp .NE. 0) CALL printiterationcess(exp1, F02, nfuneval,cputime1, 0, 0)
            CALL printrefset(exp1, refset%fpen, opts1%globaloptions%dim_refset,cputime1, common_vars%iter, refset_change )
 
            CALL addlocalscounteress(exp1)
            if (fin > 0) then
                cputime1 = calctimeMPI(exp1,time%starttime)             
                CALL printiterationcess(exp1, fbest(1), nfuneval,cputime1, 0, 0)
                cputime1 = calctimeMPI(exp1,time%starttime)
                !CALL printgant(exp1,cputime1,1)                       
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

            if (opts1%useroptions%iterprint .eq. 1 ) then
                cputime1 = calctimeMPI(exp1,time%starttime)
            end if

            CALL printrefset(exp1, refset%fpen,  opts1%globaloptions%dim_refset, cputime1, common_vars%iter, refset_change )

          end if
        end do  
        
!---------------------------------------------------------------------------------------------------------        
! GATHER THE RESULTS IN MPI EXECUTION CASE
!---------------------------------------------------------------------------------------------------------
#ifdef MPI2
        CALL returnminlocelement(exp1, fbest, F02, i, common_vars%idp)
        CALL cooperativebcastelement(exp1,xbest,common_vars%nvar,i)
        CALL returnsumelementlong(exp1, nfuneval, nfunevaltotal)
        CALL  returnavgelementintdouble( exp1, common_vars%iter, results%totaliter, common_vars%NPROC-1)
        CALL  returnminelement(exp1, results%timevtr, results%totalvtr)
        
        CALL mpibarrieress()
        CALL ending_solver_message(common_vars%idp, fin)
        CALL mpibarrieress()  
        
        cputime1 = calctimeMPI(exp1,time%starttime)
        outresult = 1
        results%timetotal = REAL(cputime1,KIND=C_DOUBLE)
        CALL select_better_solution(results, opts1, refset, xbest, fbest,local_solver_var%use_bestx, nfuneval, cputime1, fin, exp1)
        fbest(1) = F02
#endif       

        CALL detranslationinterface( xbest, exp1 )
        
        CALL mpibarrieress() 
        CALL seed_recount(exp1, common_vars%idp)
        CALL mpibarrieress()          
        
        CALL updateresultsess(exp1, results1, results%timetotal, nfunevaltotal, fbest, xbest, results%totaliter, results%totalvtr)
        CALL settotaltime(results1, results%timetotal)
        CALL setparalleltime(results1, time%timeparallel)
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
        
    END FUNCTION cess
#endif
    
END MODULE modcess
