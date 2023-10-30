!#define TIMEPAR 1
!
MODULE scattersearchfunctions
    USE iso_c_binding
    USE scattersearchtypes
    USE common_functions
    USE qsort_module
    USE funcevalinterface
    USE localsolver
#ifdef OPENMP
    USE omp_lib
#endif
CONTAINS
    
! ----------------------------------------------------------------------
! SUBROUTINE PROBLEM_SPECIFICATIONS
! ----------------------------------------------------------------------
    SUBROUTINE problem_specifications(exp1,problem1,opts1,nvar,ftarget,maxfunevals) 
        IMPLICIT NONE
        TYPE(problem), INTENT(INOUT) :: problem1
        TYPE(opts), INTENT(INOUT) :: opts1
        INTEGER, INTENT(INOUT) :: nvar
        INTEGER(C_LONG), INTENT(INOUT) :: maxfunevals
        REAL(C_DOUBLE), INTENT(INOUT) ::ftarget
        INTEGER(C_INT) :: empty1, empty2, empty3, chargeboundsnconst, chargeproblemargs
        INTEGER(C_INT) :: error, chargeuseroptions, chargelocaloptions, chargeglobaloptions, chargebounds, chargedimension
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        INTEGER :: i,status, logvar
        TYPE(opts) ::  default1
        CHARACTER(kind=C_CHAR, len=20), POINTER  :: del, fin, sol
        CHARACTER(kind=C_CHAR, len=20), TARGET  :: del2, fin2, sol2
        REAL (C_DOUBLE) :: tolf,tolx,prob_bound, tolc
        REAL (C_DOUBLE)  :: thfactor, maxdistfactor
        REAL (C_DOUBLE) :: balance
        INTEGER :: returninitsol,ii,numbersol,counter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: XX0
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: FF0
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: MAX_DBL
        CALL setdblmax(MAX_DBL)
 
        CALL create_new_problem(problem1)
        CALL create_new_opts(opts1)
        
        
        del2 = opts1%globaloptions%delete
        if (ALLOCATED(opts1%localoptions%finish)) THEN
            fin2 = opts1%localoptions%finish(1)
        END IF
        fin2=""
        sol2 = opts1%globaloptions%delete
        del => del2
        fin => fin2
        sol => sol2
        error = chargedimension(exp1,nvar) 
        problem1%empty = 0
        opts1%empty=0
        opts1%useroptions%maxeval = maxfunevals
        if (.not. ALLOCATED(problem1%vtr) ) ALLOCATE(problem1%vtr(1))
        problem1%vtr(1) = ftarget

        if (.not. ALLOCATED(opts1%useroptions%log_var )) ALLOCATE(opts1%useroptions%log_var(nvar))
        empty1 = chargeuseroptions(exp1,opts1%useroptions%maxtime, opts1%useroptions%weight, tolc, &
            prob_bound, opts1%useroptions%nstuck_solution, opts1%useroptions%strategy, &
            opts1%useroptions%inter_save, opts1%useroptions%iterprint, opts1%useroptions%plot,  opts1%useroptions%log_var , &
            opts1%useroptions%init_point, logvar )

        numbersol = returninitsol(exp1)
        counter=1
! INIT SOLUTIONS
        IF (numbersol .GT. 0) THEN
                ALLOCATE(XX0(nvar , numbersol))
                ALLOCATE(FF0(numbersol))
                do ii=1,numbersol
                        CALL loadinitpoints(exp1, XX0(:,ii), ii-1,  FF0 )
                end do

                ALLOCATE(problem1%X0(nvar , numbersol))
                do ii=1,numbersol
                        if (FF0(ii) .NE. MAX_DBL ) then
                        problem1%X0(:,counter)=XX0(:,ii)
                        counter=counter+1
                        end if
                end do

                ALLOCATE(problem1%F0(counter-1))
                do ii=1,numbersol
                        if (FF0(ii) .EQ. MAX_DBL ) then
                        problem1%X0(:,counter)=XX0(:,ii)
                        counter=counter+1
                        end if
                end do

                counter=1
                do ii=1,numbersol
                if (FF0(ii) .NE. MAX_DBL ) then
                        problem1%F0(counter)=FF0(ii)
                        counter=counter+1
                end if
                end do

                DEALLOCATE(XX0)
                DEALLOCATE(FF0)


                if ( logvar .eq. 1) then
                        do ii=1,numbersol
                              !CALL converttolog2(problem1%X0(:,ii), nvar,
                              !opts1%useroptions%log_var)
                              !CALL converttolog3(problem1%X0(:,ii), nvar,
                              !opts1%useroptions%log_var, exp1)
                              CALL convtrans(problem1%X0(:,ii), nvar, exp1)
                        end do
                end if


        END IF
! END SOLUTIONS


        if ( logvar .eq. 0 ) DEALLOCATE(opts1%useroptions%log_var) 
            
            opts1%useroptions%tolc = tolc 
            opts1%useroptions%prob_bound = prob_bound
        if ( empty1 .eq. 1 ) then
            opts1%useroptions%empty = 0
        else 
            opts1%useroptions%empty = 1
        end if
       
        empty2 = chargeglobaloptions(exp1,opts1%globaloptions%dim_refset,opts1%globaloptions%ndiverse,opts1%globaloptions%initiate,&
            opts1%globaloptions%combination, opts1%globaloptions%regenerate, del, &
            opts1%globaloptions%intens, tolf, opts1%globaloptions%diverse_criteria, &
            tolx,opts1%globaloptions%n_stuck)
            
            opts1%globaloptions%tolf = tolf
            opts1%globaloptions%tolx = tolx
            
        if ( empty2 .eq. 1 ) then 
            opts1%globaloptions%empty = 0
        else 
            opts1%globaloptions%empty = 1  
        end if

        empty3 = chargelocaloptions(exp1,opts1%localoptions%tol,opts1%localoptions%iterprint,opts1%localoptions%n1, &
            opts1%localoptions%n2,balance, fin, opts1%localoptions%bestx, &
            opts1%localoptions%merit_filter, opts1%localoptions%distance_filter, thfactor, &
            maxdistfactor, opts1%localoptions%wait_maxdist_limit, opts1%localoptions%wait_th_limit, &
            sol, opts1%localoptions%threshold_local)
            
            opts1%localoptions%thfactor = thfactor
            opts1%localoptions%maxdistfactor = maxdistfactor
            opts1%localoptions%balance = balance
            
       
        if ( empty3 .eq. 1 ) then 
            opts1%localoptions%empty = 0
        else
            opts1%localoptions%empty = 1   
        end if
        
        if ( (empty3 .eq. 1) .or. (empty2 .eq. 1) .or. (empty1 .eq. 1) ) then
            opts1%empty = 0
        else
            opts1%empty = 1
        end if
            
        if ( opts1%localoptions%empty .eq. 0 ) then
            opts1%localoptions%solver = sol
            if ( fin .NE.  "" ) THEN
            if (.not.ALLOCATED(opts1%localoptions%finish))  ALLOCATE(opts1%localoptions%finish(1))
                opts1%localoptions%finish(1) = fin
            END IF            
            opts1%globaloptions%delete = del
            opts1%localoptions%extrap%empty= 1
        end if
       
        if (.not. ALLOCATED(problem1%XU))  ALLOCATE(problem1%XU(nvar))
        if (.not. ALLOCATED(problem1%XL))  ALLOCATE(problem1%XL(nvar))
        error = chargeproblemargs(exp1,problem1%ineq,problem1%int_var,problem1%bin_var,problem1%neq)

        if ( problem1%ineq .GT. 0 )  then 
            ALLOCATE(problem1%CU(problem1%ineq))
            ALLOCATE(problem1%CL(problem1%ineq))
            error = chargeboundsnconst(exp1, problem1%XU, problem1%XL, nvar, problem1%CU, problem1%CL, problem1%ineq )

        else
            error = chargebounds(exp1, problem1%XU, problem1%XL, nvar )
        end if
        
            
        IF ((opts1%useroptions%maxeval < 0) .AND. (opts1%useroptions%maxtime < 0)) THEN
            PRINT *, 'WARNING:Either opts.maxeval or opts.maxtime must be defined as a stop criterion '
            PRINT *, 'Define any of these options and rerun '
            CALL EXIT(status)
        ELSE
            IF (opts1%useroptions%maxeval < 0) THEN
                opts1%useroptions%maxeval = INT(1e12, KIND = C_LONG_LONG)
            END IF

            IF (opts1%useroptions%maxtime < 0) THEN
                opts1%useroptions%maxtime = INT(1e12, KIND = C_LONG_LONG)
            END IF
        END IF

        ! CARGA DE PARAMETROS POR DEFECTO 
        default1 = ssm_default()
        opts1 = ssm_optset(default1, opts1)
        
        !CALL print_problems(problem1)
        !CALL print_opts(opts1)
            
    END SUBROUTINE problem_specifications
        
    
    
! ----------------------------------------------------------------------
! SUBROUTINE DESTROY_PROBLEM
! ----------------------------------------------------------------------
    SUBROUTINE destroy_problem(problem1)
        IMPLICIT NONE
        TYPE(problem), INTENT(INOUT) :: problem1

        if (ALLOCATED(problem1%XL) ) then
            DEALLOCATE(problem1%XL)
        end if
        if (ALLOCATED(problem1%XU) ) then
            DEALLOCATE(problem1%XU)
        end if
        if (ALLOCATED(problem1%CL) ) then
            DEALLOCATE(problem1%CL)
        end if
        if (ALLOCATED(problem1%CU) ) then
            DEALLOCATE(problem1%CU)
        end if
        if (ALLOCATED(problem1%X0) ) then
            DEALLOCATE(problem1%X0)
        end if
        if (ALLOCATED(problem1%F0) ) then
            DEALLOCATE(problem1%F0)
        end if
        
        if (ALLOCATED(problem1%vtr)) DEALLOCATE(problem1%vtr)

    END SUBROUTINE destroy_problem
    

! ----------------------------------------------------------------------
! SUBROUTINE INIT_COMMON_VARS
! ----------------------------------------------------------------------
SUBROUTINE init_common_vars(exp1, common_vars, opts1)
    IMPLICIT NONE
    TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: INF
    TYPE(C_PTR), INTENT(INOUT) :: exp1
    TYPE(opts), INTENT(INOUT) :: opts1
    
    CALL setdblmax(INF)
    common_vars%inititer = 0
    common_vars%init_comm  = 0
    common_vars%parallel = 1
    common_vars%iter = 0
    common_vars%stopOptimization = 0
    common_vars%nfuneval = 0
    common_vars%fin = 0
    common_vars%nreset = 0
    common_vars%last_send = INF
    common_vars%lastevals = 0
    common_vars%index_cooperative = -1
    common_vars%ncounter_slave_recv_sol = 0
    common_vars%ncounter_slave_send_sol = 0
    common_vars%pending_adaptation = 0
    common_vars%out_solver = 0
    if (opts1%globaloptions%combination .eq. 1) then
       common_vars%nrand = common_vars%nvar
    else if (opts1%globaloptions%combination .eq. 2) then
       common_vars%nrand = 1
    end if

#ifdef MPI2
    CALL setNPROC(exp1, common_vars%NPROC)
    CALL chargeid(exp1,common_vars%idp)
#endif

END SUBROUTINE    


! ----------------------------------------------------------------------
! SUBROUTINE INIT_COMMON_VARS
! ----------------------------------------------------------------------
SUBROUTINE init_time_vars(time)
    IMPLICIT NONE
    TYPE(time_ess), INTENT(INOUT) :: time

    time%timeparallel = 0.0
    time%localsolvertime =  0.0

END SUBROUTINE
    
! ----------------------------------------------------------------------
! SUBROUTINE DESTROY_OPTS
! ----------------------------------------------------------------------
    SUBROUTINE destroy_opts(opts1)
        IMPLICIT NONE
        TYPE(opts), INTENT(INOUT) :: opts1

        if (ALLOCATED(opts1%useroptions%log_var) ) then
            DEALLOCATE(opts1%useroptions%log_var)
        end if
        if (ALLOCATED(opts1%localoptions%finish) ) then
            DEALLOCATE(opts1%localoptions%finish)
        end if

    END SUBROUTINE destroy_opts
    
    
    
! ----------------------------------------------------------------------
! FUNCTION SSM_DEFAULT
! ----------------------------------------------------------------------
    TYPE(opts) FUNCTION ssm_default()
        IMPLICIT NONE
        TYPE(opts) :: default1
        
        default1%empty = 0
        
        default1%useroptions%empty = 0
        default1%useroptions%nstuck_solution = 20
        default1%useroptions%maxeval = 1000
        default1%useroptions%maxtime = 60
        default1%useroptions%iterprint = 1
        default1%useroptions%plot = 0
        default1%useroptions%weight = 1d6
        default1%useroptions%tolc = REAL(1d-5,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        default1%useroptions%prob_bound = REAL(0.5d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        default1%useroptions%strategy = 0
        default1%useroptions%inter_save = 0
        default1%useroptions%init_point = 0

        default1%globaloptions%empty = 0
        default1%globaloptions%dim_refset = -1
        default1%globaloptions%ndiverse = -1
        default1%globaloptions%initiate = 1
        default1%globaloptions%combination = 1
        default1%globaloptions%regenerate = 3
        default1%globaloptions%delete = 'standard'
        default1%globaloptions%intens = 10
        default1%globaloptions%tolf = REAL(1e-4,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        default1%globaloptions%diverse_criteria = 1
        default1%globaloptions%tolx = REAL(1e-3,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        default1%globaloptions%n_stuck = 0

        default1%localoptions%empty = 1
        default1%localoptions%solver = 'dhc'
        default1%localoptions%tol = 2
        default1%localoptions%iterprint = 0
        default1%localoptions%n1 = -1
        default1%localoptions%n2 = -1
        default1%localoptions%balance = REAL(-1d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        

        
        default1%localoptions%bestx = 0
        default1%localoptions%merit_filter = 1
        default1%localoptions%distance_filter = 1
        default1%localoptions%thfactor = REAL(0.2,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        default1%localoptions%maxdistfactor = REAL(0.2,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        default1%localoptions%wait_maxdist_limit = 20
        default1%localoptions%wait_th_limit = 20
        default1%localoptions%threshold_local = 1d-1
        default1%localoptions%extrap%empty = 1

        ssm_default = default1
    END FUNCTION ssm_default
        
    
    
! ----------------------------------------------------------------------
! SUBROUTINE CREATE_NEW_OPTS
! ----------------------------------------------------------------------
    SUBROUTINE create_new_opts(opts1) 
        IMPLICIT NONE
        TYPE(opts), INTENT(INOUT) :: opts1
        
        opts1%empty = 1
        
        opts1%useroptions%empty = 1
        opts1%useroptions%nstuck_solution = -1
        opts1%useroptions%maxeval = -1
        opts1%useroptions%maxtime = -1
        opts1%useroptions%iterprint = -1
        opts1%useroptions%plot = -1
        opts1%useroptions%weight = -1
        opts1%useroptions%tolc = -1
        opts1%useroptions%prob_bound = -1
        opts1%useroptions%strategy = -1
        opts1%useroptions%inter_save = -1
        opts1%useroptions%init_point = -1

        opts1%globaloptions%empty = 1
        opts1%globaloptions%dim_refset = -1
        opts1%globaloptions%ndiverse = -1
        opts1%globaloptions%initiate = -1
        opts1%globaloptions%combination = -1
        opts1%globaloptions%regenerate = -1
        opts1%globaloptions%delete = " "
        opts1%globaloptions%intens = -1
        opts1%globaloptions%tolf = -1
        opts1%globaloptions%diverse_criteria = -1
        opts1%globaloptions%tolx = -1
        opts1%globaloptions%n_stuck = -1

        opts1%localoptions%empty = 1
        opts1%localoptions%solver = " "
        opts1%localoptions%tol = -1
        opts1%localoptions%iterprint = -1
        opts1%localoptions%n1 = -1
        opts1%localoptions%n2 = -1
        opts1%localoptions%balance = -1d0
        opts1%localoptions%bestx = -1
        opts1%localoptions%merit_filter = -1
        opts1%localoptions%distance_filter = -1
        opts1%localoptions%thfactor = -1
        opts1%localoptions%maxdistfactor = -1
        opts1%localoptions%wait_maxdist_limit = -1
        opts1%localoptions%wait_th_limit = -1
        opts1%localoptions%extrap%empty = 1
        opts1%localoptions%threshold_local = -1d0
    END SUBROUTINE create_new_opts
        
    
    
! ----------------------------------------------------------------------
! SUBROUTINE CREATE_NEW_PROBLEM
! ----------------------------------------------------------------------
    SUBROUTINE create_new_problem(problem1)
        IMPLICIT NONE
        TYPE(problem), INTENT(INOUT) :: problem1
        
        problem1%empty = 1        
        problem1%f = ''
        problem1%neq = -1
        problem1%int_var = -1
        problem1%bin_var = -1
        
    END SUBROUTINE create_new_problem
        
    
    
! ----------------------------------------------------------------------
! FUNCTION SSM_OPTSET
! ----------------------------------------------------------------------
    TYPE(opts) FUNCTION ssm_optset(default1,opts1)
        IMPLICIT NONE
        TYPE(opts), INTENT(INOUT):: default1, opts1
        INTEGER :: DIMEN(2)
        
        ! Aquí hai un código para recolocar o campo local da estructura de matlab  opts na última posición da estructura.
        ! En principio non vexo sentido en fortran
        IF (opts1%empty .eq. 0) THEN
            default1%empty= opts1%empty 
            if (opts1%useroptions%empty .eq. 0) then

                default1%useroptions%empty = opts1%useroptions%empty 
                
                if  (ALLOCATED(opts1%useroptions%log_var)) then
                    ALLOCATE(default1%useroptions%log_var(size(opts1%useroptions%log_var)))
                    default1%useroptions%log_var = opts1%useroptions%log_var
                end if
                if  (opts1%useroptions%nstuck_solution >= 0) then
                    default1%useroptions%nstuck_solution = opts1%useroptions%nstuck_solution
                end if                
                if  (opts1%useroptions%maxeval >= 0) then
                    default1%useroptions%maxeval = opts1%useroptions%maxeval
                end if
                if  (opts1%useroptions%maxtime >= 0) then
                    default1%useroptions%maxtime = opts1%useroptions%maxtime
                end if
                if  (opts1%useroptions%iterprint >= 0) then
                    default1%useroptions%iterprint = opts1%useroptions%iterprint
                end if 
                if  (opts1%useroptions%plot >= 0) then
                    default1%useroptions%plot = opts1%useroptions%plot
                end if   
                if  (opts1%useroptions%weight >= 0) then
                    default1%useroptions%weight = opts1%useroptions%weight
                end if     
                if  (opts1%useroptions%tolc >= 0) then
                    default1%useroptions%tolc = opts1%useroptions%tolc
                end if           
                if  (opts1%useroptions%prob_bound >= 0) then
                    default1%useroptions%prob_bound = opts1%useroptions%prob_bound
                end if      
                if  (opts1%useroptions%strategy >= 0) then
                    default1%useroptions%strategy = opts1%useroptions%strategy
                end if   
                if  (opts1%useroptions%inter_save >= 0) then
                    default1%useroptions%inter_save = opts1%useroptions%inter_save
                end if    
                if  (opts1%useroptions%init_point >= 0) then
                    default1%useroptions%init_point = opts1%useroptions%init_point
                end if   
            else
                default1%useroptions%empty = opts1%useroptions%empty
            end if    
            
            if (opts1%globaloptions%empty .eq. 0) then
                default1%globaloptions%empty = opts1%globaloptions%empty 
                if  (opts1%globaloptions%dim_refset >= 0) then
                    default1%globaloptions%dim_refset = opts1%globaloptions%dim_refset
                end if
                if  (opts1%globaloptions%ndiverse >= 0) then
                    default1%globaloptions%ndiverse = opts1%globaloptions%ndiverse
                end if     
                if  (opts1%globaloptions%initiate >= 0) then
                    default1%globaloptions%initiate = opts1%globaloptions%initiate
                end if    
                if  (opts1%globaloptions%combination > 0) then
                    default1%globaloptions%combination = opts1%globaloptions%combination
                end if  
!                if  (opts1%globaloptions%regenerate >= 0) then
!                    default1%globaloptions%regenerate = opts1%globaloptions%regenerate
!                end if 
                if  (opts1%globaloptions%delete .NE. '') then
                    default1%globaloptions%delete = opts1%globaloptions%delete
                end if    
                if  (opts1%globaloptions%intens >= 0) then
                    default1%globaloptions%intens = opts1%globaloptions%intens
                end if    
                if  (opts1%globaloptions%tolf >= 0) then
                    default1%globaloptions%tolf = opts1%globaloptions%tolf
                end if     
                if  (opts1%globaloptions%diverse_criteria >= 0) then
                    default1%globaloptions%diverse_criteria = opts1%globaloptions%diverse_criteria
                end if   
                if  (opts1%globaloptions%tolx >= 0) then
                    default1%globaloptions%tolx = opts1%globaloptions%tolx
                end if          
                if  (opts1%globaloptions%n_stuck >= 0) then
                    default1%globaloptions%n_stuck = opts1%globaloptions%n_stuck
                end if   
            else
                default1%globaloptions%empty = opts1%globaloptions%empty
            end if

            if (opts1%localoptions%empty .eq. 0) then
                default1%localoptions%empty = opts1%localoptions%empty 
                if  (opts1%localoptions%solver .NE. '') then
                    default1%localoptions%solver = opts1%localoptions%solver
                end if  
                if  (opts1%localoptions%tol >= 0) then
                    default1%localoptions%tol = opts1%localoptions%tol
                end if    
                if  (opts1%localoptions%iterprint >= 0) then
                    default1%localoptions%iterprint = opts1%localoptions%iterprint
                end if   
                if  (opts1%localoptions%n1 .GE. 0) then
                    default1%localoptions%n1 = opts1%localoptions%n1
                else
                    opts1%localoptions%empty = 1
                end if         
                if  (opts1%localoptions%n2 .GE. -1) then
                    default1%localoptions%n2 = opts1%localoptions%n2
                end if     
                if  (opts1%localoptions%balance  >= 0) then
                    default1%localoptions%balance = opts1%localoptions%balance
                end if     
                if  (opts1%localoptions%bestx  >= 0) then
                    default1%localoptions%bestx = opts1%localoptions%bestx
                end if       
                if  (opts1%localoptions%merit_filter  >= 0) then
                    default1%localoptions%merit_filter = opts1%localoptions%merit_filter
                end if        
                if  (opts1%localoptions%distance_filter  >= 0) then
                    default1%localoptions%distance_filter = opts1%localoptions%distance_filter
                end if    
                if  (opts1%localoptions%thfactor  >= 0) then
                    default1%localoptions%thfactor = opts1%localoptions%thfactor
                end if   
                if  (opts1%localoptions%maxdistfactor  >= 0) then
                    default1%localoptions%maxdistfactor = opts1%localoptions%maxdistfactor
                end if   
                if  (opts1%localoptions%wait_maxdist_limit  >= 0) then
                    default1%localoptions%wait_maxdist_limit = opts1%localoptions%wait_maxdist_limit
                end if   
                if  (opts1%localoptions%wait_th_limit  >= 0) then
                    default1%localoptions%wait_th_limit = opts1%localoptions%wait_th_limit
                end if      
                !if  (opts1%localoptions%threshold_local  >= 0) then
                !    default1%localoptions%threshold_local = opts1%localoptions%threshold_local
                !end if                  
                if (ALLOCATED(opts1%localoptions%finish) ) then
                    ALLOCATE(default1%localoptions%finish(1) )
                    default1%localoptions%finish = opts1%localoptions%finish
                end if
                if (opts1%localoptions%extrap%empty .eq. 0) then
                    if (default1%localoptions%solver .EQ. 'n2fb') then
                        if (ALLOCATED(opts1%localoptions%extrap%texp)) then
                            ALLOCATE(default1%localoptions%extrap%texp(size(default1%localoptions%extrap%texp )))
                            default1%localoptions%extrap%texp = opts1%localoptions%extrap%texp
                        end if
                        DIMEN = shape(opts1%localoptions%extrap%yexp)
                        if (ALLOCATED(opts1%localoptions%extrap%yexp)) then
                            ALLOCATE(default1%localoptions%extrap%yexp(DIMEN(1),DIMEN(2)))
                            default1%localoptions%extrap%yexp = opts1%localoptions%extrap%yexp
                        end if
                    end if
                    if (default1%localoptions%solver .EQ. 'solnp') then
                        default1%localoptions%extrap%k1 = opts1%localoptions%extrap%k1
                        default1%localoptions%extrap%k2 = opts1%localoptions%extrap%k2
                        default1%localoptions%extrap%k3 = opts1%localoptions%extrap%k3
                        default1%localoptions%extrap%k4 = opts1%localoptions%extrap%k4
                    end if
                end if
            else
                default1%localoptions%empty = opts1%localoptions%empty
            end if
        else
            default1%empty = opts1%empty 
        END IF   
        
        ssm_optset = default1
    END FUNCTION ssm_optset
        
    
    
! ----------------------------------------------------------------------
! SUBROUTINE PRINT_OPTS
! ----------------------------------------------------------------------
    SUBROUTINE print_opts(opts1) 
        IMPLICIT NONE
        TYPE(opts), INTENT(IN) :: opts1
        
        PRINT *, "OPTS OPTIONS ::"
        PRINT *, " opts1%useroptions%maxeval =", opts1%useroptions%maxeval
        PRINT *, " opts1%useroptions%maxtime =", opts1%useroptions%maxtime
        PRINT *, " opts1%useroptions%iterprint =", opts1%useroptions%iterprint
        PRINT *, " opts1%useroptions%plot =", opts1%useroptions%plot
        PRINT *, " opts1%useroptions%weight =", opts1%useroptions%weight 
        PRINT *, " opts1%useroptions%tolc =", opts1%useroptions%tolc 
        PRINT *, " opts1%useroptions%prob_bound =", opts1%useroptions%prob_bound 
        PRINT *, " opts1%useroptions%strategy =", opts1%useroptions%strategy 
        PRINT *, " opts1%useroptions%inter_save =", opts1%useroptions%inter_save

        PRINT *, " opts1%globaloptions%empty =", opts1%globaloptions%empty
        PRINT *, " opts1%globaloptions%dim_refset =", opts1%globaloptions%dim_refset
        PRINT *, " opts1%globaloptions%ndiverse =", opts1%globaloptions%ndiverse
        PRINT *, " opts1%globaloptions%initiate =", opts1%globaloptions%initiate 
        PRINT *, " opts1%globaloptions%combination =", opts1%globaloptions%combination 
        PRINT *, " opts1%globaloptions%regenerate =", opts1%globaloptions%regenerate
        PRINT *, " opts1%globaloptions%delete =", opts1%globaloptions%delete
        PRINT *, " opts1%globaloptions%intens =", opts1%globaloptions%intens
        PRINT *, " opts1%globaloptions%tolf =", opts1%globaloptions%tolf 
        PRINT *, " opts1%globaloptions%diverse_criteria =", opts1%globaloptions%diverse_criteria 
        PRINT *, " opts1%globaloptions%tolx =", opts1%globaloptions%tolx 
        PRINT *, " opts1%globaloptions%n_stuck =", opts1%globaloptions%n_stuck 

        PRINT *, " opts1%localoptions%empty =", opts1%localoptions%empty 
        PRINT *, " opts1%localoptions%solver =", opts1%localoptions%solver 
        PRINT *, " opts1%localoptions%tol =", opts1%localoptions%tol 
        PRINT *, " opts1%localoptions%iterprint =", opts1%localoptions%iterprint 
        PRINT *, " opts1%localoptions%n1 =", opts1%localoptions%n1 
        PRINT *, " opts1%localoptions%n2 =", opts1%localoptions%n2 
        PRINT *, " opts1%localoptions%balance =", opts1%localoptions%balance 
        PRINT *, " opts1%localoptions%bestx =", opts1%localoptions%bestx 
        PRINT *, " opts1%localoptions%merit_filter =", opts1%localoptions%merit_filter
        PRINT *, " opts1%localoptions%distance_filter =", opts1%localoptions%distance_filter 
        PRINT *, " opts1%localoptions%thfactor =", opts1%localoptions%thfactor 
        PRINT *, " opts1%localoptions%maxdistfactor =", opts1%localoptions%maxdistfactor 
        PRINT *, " opts1%localoptions%wait_maxdist_limit =", opts1%localoptions%wait_maxdist_limit 
        PRINT *, " opts1%localoptions%wait_th_limit =", opts1%localoptions%wait_th_limit
        if (ALLOCATED(opts1%localoptions%finish) ) then
            PRINT *, " opts1%localoptions%finish =", opts1%localoptions%finish
        end if
        if (opts1%localoptions%extrap%empty .eq. 0) then
            if (opts1%localoptions%solver .EQ. 'n2fb') then
                PRINT *, " opts1%localoptions%extrap%texp =", opts1%localoptions%extrap%texp
                PRINT *, " opts1%localoptions%extrap%yexp =", opts1%localoptions%extrap%yexp
            end if
            if (opts1%localoptions%solver .EQ. 'solnp') then
                PRINT *, " opts1%localoptions%extrap%k1 =", opts1%localoptions%extrap%k1
                PRINT *, " opts1%localoptions%extrap%k2 =", opts1%localoptions%extrap%k2
                PRINT *, " opts1%localoptions%extrap%k3 =", opts1%localoptions%extrap%k3
                PRINT *, " opts1%localoptions%extrap%k4 =", opts1%localoptions%extrap%k4
            end if
        end if
        PRINT *, "END PRINT"
    END SUBROUTINE print_opts 
        
    
    
! ----------------------------------------------------------------------
! SUBROUTINE PRINT_PROBLEMS
! ----------------------------------------------------------------------
    SUBROUTINE print_problems(problem1)
        IMPLICIT NONE
        TYPE(problem), INTENT(IN) :: problem1
        
        PRINT *, "PROBLEM OPTIONS ::"
        PRINT *, " problem1%empty =",problem1%empty     
        PRINT *, " problem1%f =",problem1%f
        PRINT *, " problem1%neq =",problem1%neq
        PRINT *, " problem1%int_var =",problem1%int_var
        PRINT *, " problem1%bin_var =",problem1%bin_var 
        
    END SUBROUTINE print_problems
    
! ----------------------------------------------------------------------
! SUBROUTINES NORMALIZE_VECTOR
! ----------------------------------------------------------------------
!    SUBROUTINE normalize_vector (vector)
!        IMPLICIT NONE
!!        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT) :: vector
!        INTEGER :: i
        
!        do i=1,size(vector)
!            vector(i) = exp(vector(i))
!        end do
!        
!    END SUBROUTINE normalize_vector       
        

   
    SUBROUTINE create_combination_scheme_ess(opts1, common_vars, ess_comm, problem1)
     IMPLICIT NONE
     TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
     TYPE(ess_common_vars), INTENT(INOUT) :: ess_comm
     TYPE(opts), INTENT(INOUT) :: opts1
     TYPE(problem), INTENT(INOUT) :: problem1
     INTEGER :: auxsize, i
     INTEGER, DIMENSION(:), ALLOCATABLE :: ncomb1
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: ncomb2
     INTEGER, DIMENSION(:), ALLOCATABLE :: diff_index

        ALLOCATE(ncomb1(opts1%globaloptions%dim_refset))
        ncomb1 = (/(i, i = 1, opts1%globaloptions%dim_refset)/)
        call nchoosek_fpar(ncomb1, 2, ncomb2)
        ess_comm%MaxSubSet = (REAL(opts1%globaloptions%dim_refset)**2d0 - REAL(opts1%globaloptions%dim_refset) )/2d0
        ess_comm%MaxSubSet2 = 2d0 * ess_comm%MaxSubSet

     if (opts1%globaloptions%global_solver(1:3) .EQ. 'eSS') then
        
        CALL calc_indexes(ncomb1,ncomb2,ess_comm%index1,ess_comm%index2,ess_comm%index, &
                                                    diff_index,auxsize)
        CALL calc_ppp(auxsize,ess_comm%ppp,opts1%globaloptions%dim_refset,diff_index)
        CALL calc_hyper(common_vars%nvar, ess_comm%hyper_x_L, ess_comm%hyper_x_U, &
                                                    ess_comm%MaxSubSet, problem1%XU, problem1%XL )
     end if

        IF (ALLOCATED(diff_index)) DEALLOCATE(diff_index)
        IF (ALLOCATED(ncomb1)) DEALLOCATE(ncomb1)
        IF (ALLOCATED(ncomb2)) DEALLOCATE(ncomb2)

    END SUBROUTINE create_combination_scheme_ess
  
 
    SUBROUTINE deallocate_scheme_ess(ess_comm)
     IMPLICIT NONE
     TYPE(ess_common_vars), INTENT(INOUT) :: ess_comm
        if (ALLOCATED(ess_comm%ppp))       DEALLOCATE(ess_comm%ppp)
        if (ALLOCATED(ess_comm%hyper_x_L)) DEALLOCATE(ess_comm%hyper_x_L)
        if (ALLOCATED(ess_comm%hyper_x_U)) DEALLOCATE(ess_comm%hyper_x_U)
        if (ALLOCATED(ess_comm%index1))    DEALLOCATE(ess_comm%index1)
        if (ALLOCATED(ess_comm%index2))    DEALLOCATE(ess_comm%index2)
        if (ALLOCATED(ess_comm%index))     DEALLOCATE(ess_comm%index)

    END SUBROUTINE deallocate_scheme_ess 

! ----------------------------------------------------------------------
! SUBROUTINES log_vector
! ----------------------------------------------------------------------
    SUBROUTINE log_vector (vector)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT) :: vector
        INTEGER :: i
        
        do i=1,size(vector)
            if (vector(i) .eq. 0) then
                vector(i) = REAL(0.001d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
            end if
            vector(i) = log(vector(i))
        end do
        
    END SUBROUTINE log_vector  
        
    
    
! ----------------------------------------------------------------------
! SUBROUTINES ssm_beyond
! ----------------------------------------------------------------------
! z1 --> father in Refset
! z2 --> solution generate -- Child
    SUBROUTINE ssm_beyond ( exp1, z1, z2, z2_val, nfuneval, &
        fitnessfunction,nrand,nconst, outfunction, opts1, problem1, idopenmp)
        
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: z1, z2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: z2_val
        INTEGER, INTENT(IN) :: idopenmp
        TYPE(opts), INTENT(IN) :: opts1
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(outfuction), INTENT(INOUT) :: outfunction
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        INTEGER, INTENT(IN) :: nconst, nrand
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: denom
        INTEGER ::  continuar, n_improve, sizer
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: delta
        INTEGER, DIMENSION(:), ALLOCATABLE :: index_result, index_out_bound
        CHARACTER(len = 2) :: typesearch
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: aux_zv, xnew, xnew2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: zv, randommatrix
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: random, cputime3
        INTEGER(C_INT) :: checkcooperativemigrationcriteriacessinner
        TYPE(outfuction):: outf
        INTEGER :: nvar, mig
        
        
        nvar = size(problem1%XL)
        continuar = 1
        denom = REAL(1,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        n_improve = 1
        
        ALLOCATE(delta(size(z1)))
        ALLOCATE(zv(nvar,2))
        ALLOCATE(xnew(nvar))
        ALLOCATE(xnew2(nvar))
        ALLOCATE(randommatrix(nrand,1))
            
        mig=0
        do while ((continuar .eq. 1) .AND. (mig .NE. 1))

#ifdef MPI2 
          !  mig =  checkcooperativemigrationcriteriacessinner(exp1)
#endif 
            if (mig .NE. 1) then
            ! CREAMOS O RECTÁNGULO
            delta= ( z2-z1) / denom
            zv(:,1)=z2
            zv(:,2)=z2+delta
            
            
            ALLOCATE(aux_zv(size(z2)))
            
            aux_zv=zv(:,2)
            ! Comprobamos os bounds
            typesearch = "LT"
            ALLOCATE(index_result(nvar))
            index_result = compare_multiple_vector(aux_zv,problem1%XL, sizer,typesearch)   
            ALLOCATE(index_out_bound(size(index_result)))
            CALL index_of_ones(index_result,index_out_bound)
            CALL reajust_index(index_out_bound)
            DEALLOCATE(index_result)
            if (ALLOCATED(index_out_bound)) then
                random = RETURN_RANDOM_NUMBER(exp1)
                if (random>opts1%useroptions%prob_bound) then
                    zv(index_out_bound,2)=problem1%XL(index_out_bound)
                end if
                DEALLOCATE(index_out_bound)
            end if
            
            typesearch = "GT"
            ALLOCATE(index_result(size(problem1%XU)))
            index_result = compare_multiple_vector(aux_zv,problem1%XU, sizer,typesearch)   
            ALLOCATE(index_out_bound(size(index_result)))
            CALL index_of_ones(index_result,index_out_bound) 
            CALL reajust_index(index_out_bound)
            DEALLOCATE(index_result)
            if (ALLOCATED(index_out_bound)) then
                random = RETURN_RANDOM_NUMBER(exp1)
                if (random>opts1%useroptions%prob_bound) then
                    zv(index_out_bound,2)=problem1%XU(index_out_bound)
                end if
                DEALLOCATE(index_out_bound)
            end if      
            
            DEALLOCATE(aux_zv)
            
            CALL random_Matrix_SEED_10(exp1,randommatrix,nrand,1)
            
            xnew = zv(:,1)+(zv(:,2)-zv(:,1))*randommatrix(:,1)
            xnew2 = xnew
            outf = ssm_evalfc(exp1,fitnessfunction, xnew ,problem1, opts1, nconst,idopenmp)
            nfuneval = nfuneval + 1
            
            xnew = xnew2
            
            if (outf%include .eq. 1) then
                
                if (outf%value_penalty < z2_val) then
                    z1=z2
                    z2=xnew
                    z2_val=outf%value_penalty
                    n_improve = n_improve+1

                    if ( n_improve .eq. 2 ) then
                        ! vas haciendo el rectángulo más grande mientras no mejore
                        denom=denom/REAL(2.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
                        n_improve=0
                    end if
                    if (.not. ALLOCATED(outfunction%x) ) ALLOCATE(outfunction%x(size(outf%x)))
                    outfunction%x = outf%x
                    outfunction%value_penalty = outf%value_penalty
                    outfunction%value = outf%value
                    outfunction%pena = outf%pena
                    
                    if (nconst .GT. 0 ) then 
                        if (.not. ALLOCATED(outfunction%nlc) )   ALLOCATE(outfunction%nlc(size(outf%nlc)))
                        outfunction%nlc = outf%nlc
                    end if
                    outfunction%include = outf%include
                else
                    continuar = 0
                end if
            else
                continuar = 0
            end if
            
            end if
            
        end do
        
        DEALLOCATE(delta)
        if (ALLOCATED (zv) ) DEALLOCATE(zv)
        DEALLOCATE(xnew)
        DEALLOCATE(randommatrix)    
            
        
        
    END SUBROUTINE ssm_beyond
    
      
    
! ----------------------------------------------------------------------
! SUBROUTINES create_init_solutions
! ----------------------------------------------------------------------    
    SUBROUTINE create_init_solutions(exp1,opts1, solutions,nvar,xl_log,xu_log)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: solutions
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:),   ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: ssss
        TYPE(opts), INTENT(IN) :: opts1
        INTEGER, INTENT(IN) :: nvar
        INTEGER :: j,i, sizet, error, bench;
        INTEGER(C_INT) :: openhdf5solutions, getbench
        !INTEGER(C_INT) :: savehdf5solutions
        !We generate 5 solutions
        sizet = opts1%globaloptions%ndiverse + 5
        !if (ALLOCATED(solutions)) DEALLOCATE(solutions)
        
        ALLOCATE( solutions(nvar, sizet))
        i=1
        j=1
        solutions=0d0
       ! solutions = reshape((/ ((REAL(0.0d0 , KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), i = 1, nvar) &
       !         , j = 1, sizet) /), (/ nvar, sizet /))
        
        ALLOCATE(temp(nvar))
        
        !if (opts1 % useroptions % init_point .eq. 0) then
            do i = 1, 5
                CALL random_Vector_5(exp1,temp, nvar, xl_log, xu_log, i)
                if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0) then
                        CALL converttonormal2(temp, nvar, opts1%useroptions%log_var )
                end if
                solutions(:,i) = temp
            end do
            
            
            do i = 6, sizet
                CALL random_Vector_SEED(exp1,temp, nvar, xl_log, xu_log)
                if (ALLOCATED(opts1 % useroptions % log_var) .and. size(opts1%useroptions%log_var) .GT. 0 ) then
                        CALL converttonormal2(temp, nvar, opts1%useroptions%log_var )
                end if
                solutions(:,i) = temp
            end do
            
            !error =  savehdf5solutions(solutions, nvar, sizet)
        !else 
  
            !error = openhdf5solutions(exp1, solutions,nvar,sizet)
              
            
        !end if
        
            
        DEALLOCATE(temp)
        
    END SUBROUTINE create_init_solutions
    
    
! ----------------------------------------------------------------------
! SUBROUTINES create_init_solutions
! ----------------------------------------------------------------------    
    SUBROUTINE create_init_solutions_det(exp1,opts1, solutions,nvar,xl_log,xu_log, contador, matrixR)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:),   ALLOCATABLE, INTENT(IN) :: matrixR
        INTEGER, INTENT(INOUT) :: contador       
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: solutions
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:),   ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp
        TYPE(opts), INTENT(IN) :: opts1
        INTEGER, INTENT(IN) :: nvar
        INTEGER :: j,i, sizet
        !We generate 5 solutions
        sizet = opts1%globaloptions%ndiverse + 5
        if (ALLOCATED(solutions)) DEALLOCATE(solutions)
        ALLOCATE( solutions(nvar, sizet))
        i=1
        j=1
        solutions = reshape((/ ((REAL(0.0d0 , KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), i = 1, nvar) &
                , j = 1, sizet) /), (/ nvar, sizet /))
        
        ALLOCATE(temp(nvar))
        
        !if (opts1 % useroptions % init_point .eq. 0) then
            do i = 1, 5
                CALL random_det_Vector_5(temp, nvar, xl_log, xu_log, i, contador, matrixR)
                if (ALLOCATED(opts1 % useroptions % log_var) .and. size(opts1%useroptions%log_var) .GT. 0) then
                        CALL converttonormal2(temp, nvar, opts1%useroptions%log_var)
                end if
                !CALL all_positive(temp)
                solutions(:,i) = temp
            end do
            
            
            do i = 6, sizet
                CALL random_det_Vector_SEED(temp, nvar, xl_log, xu_log, contador, matrixR)
                
                if (ALLOCATED(opts1 % useroptions % log_var) .and. size(opts1%useroptions%log_var) .GT. 0 ) then
                        CALL converttonormal2(temp, nvar, opts1%useroptions%log_var )
                end if
                !CALL all_positive(temp)
                
                solutions(:,i) = temp
            end do
            

        !end if

            
        
        DEALLOCATE(temp)
        
    END SUBROUTINE create_init_solutions_det
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES init_inequality_constraints
! ----------------------------------------------------------------------      
    SUBROUTINE init_inequality_constraints(neq, c_L,c_U)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: neq
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) ::c_U, c_L
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp
        INTEGER :: i
        
        if (neq .GT. 0) then
            if (.not.allocated(c_L)) then
                ALLOCATE(c_L(neq))
                i=1
                c_L = (/ (REAL(0.0d0, KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), i = 1, neq) /)
            else
                ALLOCATE(temp(size(c_L) + neq)) 
                i=1
                temp = (/ (/ (REAL(0.0d0, KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
                                    i = 1, neq) /), (/ (c_L(i - neq), i = neq + 1, size(c_L) + neq) /) /)
                CALL MOVE_ALLOC(temp, c_L)
            end if
            if (.not.allocated(c_U)) then
                ALLOCATE(c_U(neq))
                i=1
                c_U = (/ (REAL(0.0d0, KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), i = 1, neq) /)
            else
                ALLOCATE(temp(size(c_U) + neq))
                i=1
                temp = (/ (/ (REAL(0.0d0, KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
                                    i = 1, neq) /), (/ (c_U(i - neq), i = neq + 1, size(c_U) + neq) /) /)
                CALL MOVE_ALLOC(temp, c_U)
            end if
        else
            neq = 0
        end if
        
    END SUBROUTINE init_inequality_constraints
         
    
       
! ----------------------------------------------------------------------
! SUBROUTINES calc_dim_refset
! ----------------------------------------------------------------------     
    SUBROUTINE calc_dim_refset(dim_refset, nvar, iterprint,idp, NPROC) 
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: dim_refset
        INTEGER, INTENT(IN) :: nvar, iterprint, idp, NPROC
        REAL(KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: nnn(2), nn(1),polim(3)
        
        if (dim_refset .EQ. - 1) then

            polim = (/  REAL(1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
                        REAL(-1.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), &
                        -10.0d0 * REAL( nvar/1, KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) /)
                        
            CALL SQUARE_ROOTS_POLINOM(polim, nnn)
            nn = maxval(nnn, MASK = nnn .GT. 0)
            dim_refset = CEILING(nn(1))
            if (mod(dim_refset, 2) .EQ. 1) then
                dim_refset = dim_refset + 1
            end if
#ifdef MPI2
            !if ((iterprint .EQ. 1) .AND. (idp .EQ. 0)) then
            !    print *, "Refset size automatically calculated:", dim_refset, nvar
            !end if            
#else
            if (iterprint .EQ. 1) then
                print *, "Refset size automatically calculated:", dim_refset, nvar
            end if
#endif
        else 
            if (iterprint .EQ. 1) then
                print *, "Refset size :", dim_refset
            end if            
        end if
    END SUBROUTINE calc_dim_refset
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES calc_ndiverse
! ----------------------------------------------------------------------      
    SUBROUTINE calc_ndiverse(ndiverse, nvar, iterprint,idp)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: ndiverse
        INTEGER, INTENT(IN) :: nvar, iterprint,idp
        
        if (ndiverse .EQ. - 1) then
            ndiverse = 10 * nvar
#ifdef MPI2
        !    if ((iterprint .EQ. 1) .AND. (idp .EQ. 0)) then
        !        print *, "Number of diverse solutions automatically calculated::", ndiverse
        !    end if
#else            
            if (iterprint .EQ. 1) then
                print *, "Number of diverse solutions automatically calculated::", ndiverse
            end if
#endif

        end if
    END SUBROUTINE calc_ndiverse
        
! ----------------------------------------------------------------------
! SUBROUTINES calc_ndiverse_hyper
! ----------------------------------------------------------------------
    SUBROUTINE calc_ndiverse_hyper(ndiverse, nvar, iterprint,idp, dim_refset)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: ndiverse
        INTEGER, INTENT(IN) :: nvar, iterprint, idp, dim_refset

        if (ndiverse .EQ. - 1) then
            ndiverse = 10 * nvar
#ifndef MPI2
            if (iterprint .EQ. 1) then
                print *, "Number of diverse solutions automatically calculated::", ndiverse
            end if
#endif
        end if

        if (ndiverse .LT. dim_refset) then
            ndiverse = dim_refset
        end if
    END SUBROUTINE calc_ndiverse_hyper
    
    
! ----------------------------------------------------------------------
! SUBROUTINES destroy_refsettype
! ----------------------------------------------------------------------      
    SUBROUTINE destroy_refsettype(refset) 
        IMPLICIT NONE
        TYPE(Refsettype), INTENT(INOUT) :: refset
        
        if (ALLOCATED(refset%x ))       DEALLOCATE(refset%x )
        if (ALLOCATED(refset%f ))       DEALLOCATE(refset%f )      
        if (ALLOCATED(refset%fpen ))    DEALLOCATE(refset%fpen )
        if (ALLOCATED(refset%penalty )) DEALLOCATE(refset%penalty )
        if (ALLOCATED(refset%nlc ))     DEALLOCATE(refset%nlc )
        if (ALLOCATED(refset%parent_index) ) DEALLOCATE(refset%parent_index )
        if (ALLOCATED(refset%cooperative) ) DEALLOCATE(refset%cooperative)
        
    END SUBROUTINE destroy_refsettype
        
    
    
! ----------------------------------------------------------------------
! SUBROUTINES destroy_refsettype
! ----------------------------------------------------------------------      
    SUBROUTINE destroy_resultsf(results) 
        IMPLICIT NONE
        TYPE(resultsf), INTENT(INOUT) :: results
        
        if (ALLOCATED(results%f)) DEALLOCATE(results%f)
        if (ALLOCATED(results%x)) DEALLOCATE(results%x)
        if (ALLOCATED(results%time)) DEALLOCATE(results%time)
        if (ALLOCATED(results%neval )) DEALLOCATE(results%neval )
        if (ALLOCATED(results%fbest)) DEALLOCATE(results%fbest)
        if (ALLOCATED(results%xbest )) DEALLOCATE(results%xbest )
        CALL destroy_refsettype(results%refsett)
        
    END SUBROUTINE destroy_resultsf    
      
   
    
    
   
! ----------------------------------------------------------------------
! SUBROUTINES evaluate_solutions_set
! ----------------------------------------------------------------------   
   SUBROUTINE  evaluate_solutions_set(exp1,fitnessfunction,solutionset,problem1,opts1, solutions, nvar, nfuneval,nconst,&
        xl_log, xu_log)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        TYPE(Refsettype), INTENT(INOUT) :: solutionset
        TYPE(Refsettype) :: solutionset_aux
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(opts), INTENT(IN) :: opts1
        TYPE(outfuction) :: outfunct
        INTEGER, INTENT(IN) :: nconst
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: solutions
        INTEGER, INTENT(IN)  :: nvar
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER :: l_f_0, l_x_0, DIMEN(2), tempnum, counter, i, sizet
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: newsolution, entervalues
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp2, temp3
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: NAN
        INTEGER(C_INT) :: error, setNAN, continueeval
        INTEGER :: counter_eval            
          
        error = setNAN(NAN)
        
        if (ALLOCATED(problem1%F0)) then
            l_f_0 = size(problem1%F0)
        else
            l_f_0 = 0
        end if

        if (ALLOCATED(problem1%X0)) then
            DIMEN = shape(problem1%X0)
            l_x_0 = DIMEN(2)
        else
            l_x_0 = 0
        end if

        ! Put x_0 without their corresponding problem1%F0 values at the beginning of solutions
        if ( (l_x_0 - l_f_0) > 0) then
            DIMEN = shape(problem1%X0)
            
            ALLOCATE(temp2(DIMEN(1), DIMEN(2)- l_f_0))
            temp2 = problem1%X0(:,(l_f_0 + 1):DIMEN(2))
            CALL FUSION_MATRIX(temp2, solutions)
            tempnum = l_x_0 - l_f_0
            !DEALLOCATE(temp2)
        else
            tempnum = 0
        end if
         

        sizet = opts1%globaloptions%ndiverse + 5 + tempnum
        ALLOCATE(solutionset%x(nvar, sizet))
        solutionset%x = solutionset%x * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) 
        if (nconst .GT. 0 ) then 
            ALLOCATE(solutionset%nlc(nconst, sizet))
            solutionset%nlc = solutionset%nlc * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) 
        end if 
        ALLOCATE(solutionset%f(sizet))
        solutionset%f = solutionset%f * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(solutionset%fpen(sizet))
        solutionset%fpen = solutionset%fpen * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(solutionset%penalty(sizet))
        solutionset%penalty = solutionset%penalty * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        


        counter = 1
        !Evaluate the set of solutions
   
        do i = 1, sizet
            continueeval = 0
            ALLOCATE(newsolution(nvar))
            newsolution(:) = solutions(:,i)
            counter_eval = 0
            do while ((continueeval .eq. 0) .AND. (counter_eval .LT. 3))
                outfunct = ssm_evalfc(exp1, fitnessfunction, newsolution, problem1, opts1, nconst, 0)
                solutions(:,i) = outfunct%x
                nfuneval = nfuneval + 1

                if (outfunct%include .eq. 1) then
                    solutionset%x(:,counter) = outfunct%x
                    solutionset%f(counter) = outfunct%value
                    solutionset%fpen(counter) = outfunct%value_penalty
                    solutionset%penalty(counter) = outfunct%pena
                    if (nconst .GT. 0) then 
                        solutionset%nlc(:,counter) = outfunct%nlc
                    end if
                    counter = counter + 1
                    continueeval = 1
                else if ((outfunct%include .eq. 0) .AND. ( (counter-1) .LT. opts1%globaloptions%dim_refset )) then
                    counter_eval = counter_eval + 1
                    CALL random_Vector_SEED(exp1,newsolution, nvar, xl_log, xu_log)   
                    if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0 )  &
                           CALL converttonormal2(newsolution, nvar, opts1%useroptions%log_var)
                else    
                    continueeval = 1
                end if
                call destroy_outfuction(outfunct)
            end do
            if (ALLOCATED(newsolution)) DEALLOCATE(newsolution)
        end do

        !resize solutionset
        if ( (counter-1) .LT. sizet )  then
                sizet = counter - 1
                ALLOCATE(solutionset_aux%x(nvar, sizet))
                solutionset_aux%x(:,1:(counter-1))=solutionset%x(:,1:(counter-1))
                if (nconst .GT. 0 ) then
                        ALLOCATE(solutionset_aux%nlc(nconst, sizet))
                        solutionset_aux%nlc(:,1:(counter-1))=solutionset%nlc(:,1:(counter-1))
                end if
                ALLOCATE(solutionset_aux%f(sizet))
                solutionset_aux%f(1:(counter-1))=solutionset%f(1:(counter-1))
                ALLOCATE(solutionset_aux%fpen(sizet))
                solutionset_aux%fpen(1:(counter-1))=solutionset%fpen(1:(counter-1))
                ALLOCATE(solutionset_aux%penalty(sizet))
                solutionset_aux%penalty(1:(counter-1))=solutionset%penalty(1:(counter-1))
              
                DEALLOCATE(solutionset%x) 
                ALLOCATE(solutionset%x(nvar, sizet))
                if (nconst .GT. 0 ) then
                        DEALLOCATE(solutionset%nlc)
                        ALLOCATE(solutionset%nlc(nconst , sizet))
                end if
                DEALLOCATE(solutionset%f)
                ALLOCATE(solutionset%f(sizet))
                DEALLOCATE(solutionset%fpen)
                ALLOCATE(solutionset%fpen(sizet))
                DEALLOCATE(solutionset%penalty)
                ALLOCATE(solutionset%penalty(sizet))

                solutionset%x = solutionset_aux%x
                if (nconst .GT. 0 ) then
                        solutionset%nlc=solutionset_aux%nlc
                end if
                solutionset%f=solutionset_aux%f
                solutionset%fpen=solutionset_aux%fpen
                solutionset%penalty=solutionset_aux%penalty
                
                DEALLOCATE(solutionset_aux%x)
                if (nconst .GT. 0 ) then
                        DEALLOCATE(solutionset_aux%nlc)
                end if
                DEALLOCATE(solutionset_aux%f)
                DEALLOCATE(solutionset_aux%fpen)
                DEALLOCATE(solutionset_aux%penalty)
                 
        end if

        !Add points with f_0 values. We assume feasiblity
        if (l_f_0 > 0) then
            DIMEN = shape(problem1%X0)
            ALLOCATE(temp2(DIMEN(1),l_f_0))
            temp2 = problem1%X0(:,1:l_f_0)
            call FUSION_MATRIX(temp2, solutionset%x)
            call FUSION_VECTOR(problem1%F0, solutionset%f)
            call FUSION_VECTOR(problem1%F0, solutionset%fpen)
            ALLOCATE(entervalues(l_f_0))
            entervalues = entervalues * 0
            call FUSION_VECTOR(entervalues, solutionset%penalty)
            DEALLOCATE(entervalues)
            
            if (nconst .GT. 0) then
                ALLOCATE(temp3(nconst,l_f_0))
                temp3 = 1.0 * NAN
                CALL FUSION_MATRIX(temp3,solutionset%nlc)
            end if
        end if    
        
   END SUBROUTINE evaluate_solutions_set
   
  
! ----------------------------------------------------------------------
! SUBROUTINES create_refset
! ----------------------------------------------------------------------    
   SUBROUTINE create_refset ( exp1, solutionset, refset, nvar, opts1, nconst, dim_refset )
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(Refsettype), INTENT(INOUT) :: solutionset, refset
        INTEGER, INTENT(IN) :: dim_refset
        INTEGER, INTENT(IN) :: nvar, nconst
        TYPE(opts), INTENT(IN) :: opts1        
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp
        INTEGER :: first_members, last_members, DIMEN(2)
        INTEGER, DIMENSION(:), ALLOCATABLE :: indvect, ooo
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp2,temp3
        
        ALLOCATE(indvect(size(solutionset%fpen)))
        ALLOCATE(temp(size(solutionset%fpen)))
        temp = solutionset%fpen
        call QsortC_ind(temp, indvect)
        DEALLOCATE(temp)

        first_members = CEILING(REAL(dim_refset/2, KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
        last_members = dim_refset  - first_members

        ! Initialize Refset
        ALLOCATE(refset%x(nvar,first_members))
        refset%x = refset%x * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%f(first_members))
        refset%f = refset%f * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%fpen(first_members))
        refset%fpen = refset%fpen *REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%penalty(first_members))
        refset%penalty = refset%penalty * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        if (nconst .GT. 0) then
            ALLOCATE(refset%nlc(nconst,first_members))
            refset%nlc = refset%nlc * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))              
        end if
               
        ! Create the initial sorted Refset
        refset%x = solutionset%x(:,indvect(1:first_members))
        refset%f = solutionset%f(indvect(1:first_members))
        refset%fpen = solutionset%fpen(indvect(1:first_members))
        refset%penalty = solutionset%penalty(indvect(1:first_members))
        if (nconst .GT. 0) then
            refset%nlc = solutionset%nlc(:,indvect(1:first_members))
        end if

        ! The rest of Refset members are chosen randomly
        CALL delete_values_vector_INT(indvect, first_members)

        ALLOCATE (ooo(size(indvect)))
        
        CALL randperm_int(exp1,ooo)
        
        DIMEN = shape(solutionset%x)
        ALLOCATE(temp2(DIMEN(1), last_members))
        temp2 = solutionset%x(:,indvect(ooo(1:last_members)))
        CALL FUSION_MATRIX(refset%x, temp2)
        CALL MOVE_ALLOC(temp2, refset%x)


        ALLOCATE(temp(last_members))
        temp = solutionset%f(indvect(ooo(1:last_members)))
        CALL FUSION_VECTOR(refset%f, temp)
        CALL MOVE_ALLOC(temp, refset%f)

        ALLOCATE(temp(last_members))
        temp = solutionset%fpen(indvect(ooo(1:last_members)))
        CALL FUSION_VECTOR(refset%fpen, temp)
        CALL MOVE_ALLOC(temp, refset%fpen)
        
        ALLOCATE(temp(last_members))
        temp = solutionset%penalty(indvect(ooo(1:last_members)))
        CALL FUSION_VECTOR(refset%penalty, temp)
        CALL MOVE_ALLOC(temp, refset%penalty)
       
        if (nconst .GT. 0) then
           ALLOCATE(temp3(nconst, last_members))
           temp3 = solutionset%nlc(:,indvect(ooo(1:last_members)))
           CALL FUSION_MATRIX(refset%nlc, temp3)
           CALL MOVE_ALLOC(temp3, refset%nlc)
        end if
        
        if (ALLOCATED(ooo)) DEALLOCATE(ooo)
        
        ALLOCATE(refset%cooperative( opts1%globaloptions%dim_refset ))
        refset%cooperative = 0        
       
        
   END SUBROUTINE create_refset
        
   
   ! ----------------------------------------------------------------------
! SUBROUTINES create_refset_empty
! ----------------------------------------------------------------------    
   SUBROUTINE create_refset_empty ( refset, nvar, opts1, nconst, size1 )
        IMPLICIT NONE
        TYPE(Refsettype), INTENT(INOUT) :: refset
        INTEGER, INTENT(IN) :: nvar, nconst, size1
        TYPE(opts), INTENT(IN) :: opts1        

        ! Initialize Refset
        ALLOCATE(refset%x(nvar,size1))
        refset%x = refset%x * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%f(size1))
        refset%f = refset%f * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%fpen(size1))
        refset%fpen = refset%fpen *REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(refset%penalty(size1))
        refset%penalty = refset%penalty * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        if (nconst .GT. 0) then
            ALLOCATE(refset%nlc(nconst,size1))
            refset%nlc = refset%nlc * REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))              
        end if
        ALLOCATE(refset%cooperative(size1))
        refset%cooperative = 0
   END SUBROUTINE create_refset_empty
        
   
   
   
   
! ----------------------------------------------------------------------
! SUBROUTINES create_child
! ----------------------------------------------------------------------    
   SUBROUTINE create_childset ( childset, nvar, MaxSubSet2, nconst )   
       IMPLICIT NONE
       TYPE(Refsettype), INTENT(INOUT) :: childset
       INTEGER, INTENT(IN) :: nvar, nconst, MaxSubSet2
       
            ALLOCATE(childset%x(nvar,MaxSubSet2))
            ALLOCATE(childset%f(MaxSubSet2) )
            ALLOCATE(childset%fpen(MaxSubSet2) )
            ALLOCATE(childset%penalty(MaxSubSet2) )
            if (nconst .gt. 0) then
		ALLOCATE(childset%nlc(nconst,MaxSubSet2))
	    end if
            ALLOCATE(childset%parent_index(MaxSubSet2) )
            
            childset%x = 0d0
            childset%f = 0d0
            childset%fpen = 0d0
            childset%penalty = 0d0
            if (nconst .gt. 0) then
		childset%nlc = 0d0
	    end if
            childset%parent_index = 0
            
       
   END SUBROUTINE
! ----------------------------------------------------------------------
! SUBROUTINES check_best_value
! ---------------------------------------------------------------------- 
   SUBROUTINE check_best_value(refset, opts1, xbest, fbest )
        IMPLICIT NONE
        TYPE(Refsettype), INTENT(INOUT) :: refset
        TYPE(opts), INTENT(INOUT) :: opts1 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: xbest, fbest
        CHARACTER(len = 2) :: typesearch
        INTEGER, DIMENSION(:), ALLOCATABLE :: indexfind
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp
        INTEGER, DIMENSION(:), ALLOCATABLE :: iiim
        INTEGER :: iii

        typesearch = "LT"
        
        CALL find_in_vector(refset%penalty, opts1%useroptions%tolc, indexfind, typesearch)
        if (ALLOCATED(indexfind)) then
            CALL reajust_index(indexfind)
            ALLOCATE(temp(size(indexfind)))
            temp = refset%fpen(indexfind)
            if (.not. allocated(fbest)) allocate(fbest(1))
            fbest = minval(temp)
            ALLOCATE(iiim(1))
            iiim = minloc(temp)
            iii = iiim(1)
            xbest = refset%x(:,(indexfind(iii)))
            DEALLOCATE(iiim)
            DEALLOCATE(temp)
        end if
        
        IF (ALLOCATED(iiim)) DEALLOCATE(iiim)
        IF (ALLOCATED(indexfind)) DEALLOCATE(indexfind)
    END SUBROUTINE check_best_value
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES create_results
! ----------------------------------------------------------------------     
    SUBROUTINE create_results(results, opts1, xbest, fbest, nfuneval, cputime1, nvar)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: xbest, fbest
        TYPE(opts), INTENT(IN) :: opts1 
        INTEGER(C_LONG), INTENT(IN) :: nfuneval
        INTEGER, INTENT(IN):: nvar
        REAL(C_DOUBLE), INTENT(IN) :: cputime1
        TYPE(resultsf), INTENT(INOUT) :: results
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: INF
        
        CALL setdblmax(INF)
                
        ! Initialize Refset for restuls struct 
        ALLOCATE(results%refsett%x(nvar,opts1%globaloptions%dim_refset))
        results%refsett%x =  REAL ( 0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        
        ALLOCATE(results%refsett%f(opts1%globaloptions%dim_refset))
        results%refsett%f =  REAL ( 0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        
        ALLOCATE(results%refsett%fpen(opts1%globaloptions%dim_refset))
        results%refsett%fpen =  REAL ( 0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        
        ALLOCATE(results%refsett%penalty(opts1%globaloptions%dim_refset))
        results%refsett%penalty =  REAL ( 0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  
       
        ALLOCATE(results%refsett%cooperative(opts1%globaloptions%dim_refset))
        results%refsett%cooperative =  REAL ( 0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
 
        ALLOCATE(results%f(size(fbest)))
        results%f = fbest
        
        ALLOCATE(results%x(size(xbest),1))
        results%x(:,1) = xbest
        
        ALLOCATE(results%time(1))
        results%time(1) = cputime1
        
        ALLOCATE(results%neval(1))
        results%neval(1) = nfuneval
        
        results%timevtr=INF

        
        
    END SUBROUTINE create_results
    
    
! ----------------------------------------------------------------------
! SUBROUTINES sort_refset
! ----------------------------------------------------------------------     
    SUBROUTINE sort_refset(refset, opts1, indvect)   
        IMPLICIT NONE
        TYPE(Refsettype), INTENT(INOUT) :: refset
        TYPE(opts), INTENT(IN) :: opts1 
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: indvect
        
        if (ALLOCATED(indvect)) then
            DEALLOCATE(indvect)
        end if
        ALLOCATE(indvect(opts1%globaloptions%dim_refset))
        
        call QsortC_ind(refset%fpen, indvect)
        refset%f = refset%f(indvect)
        refset%penalty = refset%penalty(indvect)
        refset%x = refset%x(:,indvect)
        if (ALLOCATED(refset%nlc)) refset%nlc = refset%nlc(:,indvect)
        if (ALLOCATED(refset%cooperative)) refset%cooperative = refset%cooperative(indvect)
        if (ALLOCATED(refset%parent_index)) refset%parent_index = refset%parent_index(indvect)

    END SUBROUTINE sort_refset
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES create_candidate
! ----------------------------------------------------------------------     
    SUBROUTINE create_candidate(candidateset, refset)
        IMPLICIT NONE
        TYPE(Refsettype), INTENT(IN) :: refset
        TYPE(Refsettype), INTENT(INOUT) :: candidateset
        INTEGER :: DIMEN(2), DIMEN2(2)
        
        
        DIMEN = shape(refset%x)
        ALLOCATE(candidateset%x(DIMEN(1), DIMEN(2)))
        ALLOCATE(candidateset%f(size(refset%f)))
        ALLOCATE(candidateset%fpen(size(refset%fpen)))
        ALLOCATE(candidateset%penalty(size(refset%penalty)))
        if (ALLOCATED(refset%nlc)) then 
                DIMEN2=shape(refset%nlc)
                ALLOCATE(candidateset%nlc( DIMEN2(1)  , DIMEN(2)))
        end if

        candidateset%x = refset%x
        candidateset%f = refset%f
        candidateset%fpen = refset%fpen
        if (ALLOCATED(refset%nlc)) candidateset%nlc = refset%nlc
        candidateset%penalty = refset%penalty
        
    END SUBROUTINE create_candidate
        
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES check_duplicated_replace
! Eliminamos os elementos duplicados do REFSET, ou elementos que esten 
! moi próximos entre si.
! ----------------------------------------------------------------------     
    SUBROUTINE check_duplicated_replace (exp1, fitnessfunction, problem1, refset, opts1, parents_index1, parents_index2, &
                        xl_log, xu_log, refset_change, members_update, index1, index2, nfuneval, nvar, nconst )
            IMPLICIT NONE
            TYPE(problem), INTENT(IN) :: problem1
            TYPE(C_PTR), INTENT(INOUT) :: exp1
            TYPE(Refsettype), INTENT(INOUT) :: refset
            TYPE(opts), INTENT(IN) :: opts1 
            TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: parents_index1,&
                                                                                                                 parents_index2
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change            
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: index1, index2
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: members_update
            INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
            INTEGER, INTENT(IN) :: nvar, nconst
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
            INTEGER :: continuar, auxsize, DIMEN(2)
            INTEGER, DIMENSION(:), ALLOCATABLE :: indvect
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE  :: &
                                                                    parents_index1_values_penalty, parents_index2_values_penalty
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: denom22, AAA, temp2
            INTEGER, DIMENSION(:), ALLOCATABLE :: BBB
            CHARACTER(len = 2) :: typesearch
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: new_refset_member
            INTEGER, DIMENSION(:), ALLOCATABLE :: index_result, index_refset_out
            INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_result2
            TYPE(outfuction) :: outfunct
            INTEGER :: counter_eval
           
            continuar = 0
            do while (continuar .eq. 0)
                !Sort Refset
                !This is mandatory because the variable index is also sorted
                CALL sort_refset(refset,opts1, indvect)
                refset_change = refset_change(indvect)
                DIMEN = shape(refset%x)
                if (.not.ALLOCATED(parents_index1)) ALLOCATE(parents_index1(DIMEN(1), size(index1)))
                parents_index1 = refset%x(:,index1)
                if (.not.ALLOCATED(parents_index2)) ALLOCATE(parents_index2(DIMEN(1), size(index2)))
                parents_index2 = refset%x(:,index2)
                if (.not.ALLOCATED(parents_index1_values_penalty)) &
                ALLOCATE(parents_index1_values_penalty(size(index1)))
                parents_index1_values_penalty = refset%fpen(index1)
                
                if (.not.ALLOCATED(parents_index2_values_penalty)) &
                ALLOCATE(parents_index2_values_penalty(size(index2)))
                parents_index2_values_penalty = refset%fpen(index2)

                ALLOCATE(denom22(DIMEN(1),size(index1)))
                denom22 = max_matrix(abs_matrix(parents_index1), abs_matrix(parents_index2))

                call replace_matrix(denom22, REAL(1d0, KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
                ALLOCATE(temp2(DIMEN(1),size(index1)))
                temp2 = (parents_index1 - parents_index2)  / denom22
                
                ALLOCATE(AAA(DIMEN(1),size(index1)))
                AAA = abs_matrix(temp2)
                
                DEALLOCATE(temp2)
                typesearch = "LT"
                DIMEN = shape(AAA)
                auxsize = 0
                ALLOCATE(index_result2(DIMEN(1),DIMEN(2)))
                index_result2 = compare_matrix(AAA, REAL(1.0d-3, KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), typesearch)
                
                ALLOCATE(index_result(size(index1)))
                CALL row_one(index_result2, index_result, auxsize)
                ALLOCATE(BBB(size(index_result)))
                CALL index_of_ones(index_result, BBB)
                CALL reajust_index(BBB)
                DEALLOCATE(index_result2)
                DEALLOCATE(index_result)
                
                if (ALLOCATED(BBB)) then
                    ALLOCATE( index_result(size(BBB)))
                    index_result = index2(BBB)
                    ALLOCATE(index_refset_out(1))
                    index_refset_out = MAXVAL(index_result)
                    DEALLOCATE(index_result)
                    outfunct%include = 0
                    counter_eval = 0
                    do while ((outfunct%include .eq. 0) .AND. (counter_eval .LT. 3 ))
                        ALLOCATE(new_refset_member(nvar))
                        CALL random_Vector_SEED(exp1,new_refset_member, nvar, xl_log, xu_log)
                        
                        if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0 )  &
                                CALL converttonormal2(new_refset_member, nvar, opts1%useroptions%log_var)
                        
                        outfunct = ssm_evalfc(   exp1, fitnessfunction, new_refset_member, problem1, opts1, nconst, 0)
                        nfuneval = nfuneval + 1
                        if (outfunct%include .eq. 1) then
                            refset%x(:,index_refset_out(1)) = outfunct%x
                            refset%f(index_refset_out(1)) = outfunct%value
                            refset%fpen(index_refset_out(1)) = outfunct%value_penalty
                            if (ALLOCATED(refset%nlc)) refset%nlc(:,index_refset_out(1)) = outfunct%nlc
                            refset%penalty(index_refset_out(1)) = outfunct%pena
                            members_update(index_refset_out(1)) = 0
                        else
                            counter_eval = counter_eval + 1
                        end if
                        DEALLOCATE(new_refset_member)
                    end do
                    DEALLOCATE(index_refset_out)
                    call destroy_outfuction(outfunct)
                    
                else
                    continuar = 1
                end if

                !print *, "end compare_matrix AAA"
                DEALLOCATE(denom22)
                if ( ALLOCATED(AAA)) DEALLOCATE(AAA)
                if ( ALLOCATED(BBB)) DEALLOCATE(BBB)

            if (ALLOCATED(indvect)) DEALLOCATE(indvect)
            if (ALLOCATED(parents_index1_values_penalty)) DEALLOCATE(parents_index1_values_penalty)
            if (ALLOCATED(parents_index2_values_penalty)) DEALLOCATE(parents_index2_values_penalty)                
            end do ! Fin bucle interior
            
            
    END SUBROUTINE check_duplicated_replace        
          
    

    
! ----------------------------------------------------------------------
! SUBROUTINES check_vector_in_hyper_bounds
! ----------------------------------------------------------------------        
    SUBROUTINE check_vector_in_hyper_bounds(exp1,opts1, v, hyper_x_L, hyper_x_U )
            IMPLICIT NONE
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: hyper_x_L, hyper_x_U
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
            TYPE(C_PTR), INTENT(INOUT) :: exp1
            TYPE(opts), INTENT(IN) :: opts1 
            CHARACTER(len = 2) :: typesearch
            INTEGER :: DIMEN(2)
            INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_result2, index_result, index_result3
            INTEGER :: auxsize 
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp1            
                        
            DIMEN = shape(v)
            ALLOCATE(index_result(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result2(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result3(DIMEN(1), DIMEN(2)))
            ALLOCATE(temp(DIMEN(1), DIMEN(2)))
            
            typesearch = "LT"
            auxsize = 0
            index_result = compare_multiple_matrix(v, hyper_x_L, auxsize, typesearch)

            
            if (auxsize .GT. 0) then
                ALLOCATE(temp1(auxsize))
                call  random_Vector_SEED_10(exp1,temp1, auxsize)      
                call  set_one_matrix(temp,temp1,index_result)
                DEALLOCATE(temp1)                
                
                typesearch = "GT"
                index_result2 = compare_matrix(temp, opts1%useroptions%prob_bound, typesearch)
                
                typesearch = "EB"
                auxsize = 0
                index_result3 = compare_multiple_matrix_int(index_result, index_result2, auxsize, typesearch)

                if ( auxsize .GT. 0) then
                    call return_values_matrix_index(v, hyper_x_L, index_result3)
                end if
            end if
            DEALLOCATE(index_result)
            DEALLOCATE(index_result2)
            DEALLOCATE(index_result3)
            DEALLOCATE(temp)
            
            
            
            DIMEN = shape(v)
            ALLOCATE(index_result(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result2(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result3(DIMEN(1), DIMEN(2)))
            ALLOCATE(temp(DIMEN(1), DIMEN(2)))
            
            typesearch = "GT"
            auxsize= 0

            index_result = compare_multiple_matrix(v, hyper_x_U, auxsize, typesearch)

            if (auxsize .GT. 0) then
                ALLOCATE(temp1(auxsize))
                call  random_Vector_SEED_10(exp1,temp1, auxsize)      
                call  set_one_matrix(temp,temp1,index_result)
                DEALLOCATE(temp1)                     
                typesearch = "GT"
                index_result2 = compare_matrix(temp, opts1%useroptions%prob_bound, typesearch)
                
                typesearch = "EB"
                auxsize = 0
                index_result3 = compare_multiple_matrix_int(index_result, index_result2, auxsize, typesearch)

                if ( auxsize .GT. 0) then
                    call return_values_matrix_index(v, hyper_x_U, index_result3)
                end if
            end if
            

            
            DEALLOCATE(index_result)
            DEALLOCATE(index_result2)
            DEALLOCATE(index_result3)
            DEALLOCATE(temp)            

        
    END SUBROUTINE check_vector_in_hyper_bounds
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES check_vector_in_hyper_bounds
! ----------------------------------------------------------------------        
    SUBROUTINE check_vector_in_hyper_bounds_det(exp1,opts1, v, hyper_x_L, hyper_x_U, contador, matrixR )
            IMPLICIT NONE
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: hyper_x_L, hyper_x_U
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
            TYPE(C_PTR), INTENT(INOUT) :: exp1
            TYPE(opts), INTENT(IN) :: opts1 
            CHARACTER(len = 2) :: typesearch
            INTEGER :: DIMEN(2)
            INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_result2, index_result, index_result3
            INTEGER :: auxsize 
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp1
            INTEGER, INTENT(INOUT) :: contador  
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: matrixR
            REAL(C_DOUBLE) :: cputime4, cputime5

                    
            DIMEN = shape(v)
            ALLOCATE(index_result(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result2(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result3(DIMEN(1), DIMEN(2)))
            ALLOCATE(temp(DIMEN(1), DIMEN(2)))
            
            typesearch = "LT"
            auxsize = 0
            call cpu_time(cputime4) 
            index_result = compare_multiple_matrix(v, hyper_x_L, auxsize, typesearch)
            call cpu_time(cputime5) 
            
            !print *, "  time evals1.1", cputime5-cputime4   
            if (auxsize .GT. 0) then
                ALLOCATE(temp1(auxsize))
                call  random_det_Vector_SEED_10(temp1, auxsize, contador, matrixR)      
                call  set_one_matrix(temp,temp1,index_result)
                DEALLOCATE(temp1)
                
                typesearch = "GT"
                !call cpu_time(cputime4) 
                index_result2 = compare_matrix(temp, opts1%useroptions%prob_bound, typesearch)
                !call cpu_time(cputime5) 
                
                !print *, "  time evals1.2", cputime5-cputime4   
                typesearch = "EB"
                !call cpu_time(cputime4) 
                auxsize = 0
                index_result3 = compare_multiple_matrix_int(index_result, index_result2, auxsize, typesearch)
                !call cpu_time(cputime5) 
                !print *, "  time evals1.3", cputime5-cputime4   
                if ( auxsize .GT. 0) then
                    call return_values_matrix_index(v, hyper_x_L, index_result3)
                end if
                
            end if
            DEALLOCATE(index_result)
            DEALLOCATE(index_result2)
            DEALLOCATE(index_result3)
            DEALLOCATE(temp)
            
            
            
            DIMEN = shape(v)
            ALLOCATE(index_result(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result2(DIMEN(1), DIMEN(2)))
            ALLOCATE(index_result3(DIMEN(1), DIMEN(2)))
            ALLOCATE(temp(DIMEN(1), DIMEN(2)))
            
            typesearch = "GT"
            auxsize= 0
            !call cpu_time(cputime4) 
            index_result = compare_multiple_matrix(v, hyper_x_U, auxsize, typesearch)
            !call cpu_time(cputime5) 
            !print *, "  time evals2.1", cputime5-cputime4
            if (auxsize .GT. 0) then
                ALLOCATE(temp1(auxsize))
                call  random_det_Vector_SEED_10(temp1, auxsize, contador, matrixR)      
                call  set_one_matrix(temp,temp1,index_result)
                DEALLOCATE(temp1)
                
                typesearch = "GT"
                !call cpu_time(cputime4) 
                index_result2 = compare_matrix(temp, opts1%useroptions%prob_bound, typesearch)
                call cpu_time(cputime5) 
                !print *, "  time evals2.2", cputime5-cputime4
                typesearch = "EB"
                auxsize = 0
                !call cpu_time(cputime4) 
                index_result3 = compare_multiple_matrix_int(index_result, index_result2, auxsize, typesearch)
                !call cpu_time(cputime5) 
                !print *, "  time evals2.3", cputime5-cputime4
                if ( auxsize .GT. 0) then
                    call return_values_matrix_index(v, hyper_x_U, index_result3)
                end if
            end if
            

            
            DEALLOCATE(index_result)
            DEALLOCATE(index_result2)
            DEALLOCATE(index_result3)
            DEALLOCATE(temp)            

        
    END SUBROUTINE check_vector_in_hyper_bounds_det    
    
    

    
! ----------------------------------------------------------------------
! SUBROUTINES generate_new_comb_matrix
! Esta funcion crea novas solucion a partir de combinacion e outras operación
! entre outros membros do REFSET 
! ----------------------------------------------------------------------     
    SUBROUTINE generate_new_comb_matrix(exp1, new_comb, ppp, MaxSubSet, nrand, v1,v2,v3)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: MaxSubSet
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ppp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  new_comb
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: v1, v2, v3
        INTEGER, INTENT(IN) :: nrand
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: array1, array2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: new_comb1, new_comb2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: random
        
        if (nrand > 1) then
                !print *, " - generate_new_comb_matrix -- new_comb = ", MaxSubSet * 2
                ALLOCATE(array1(nrand, CEILING(MaxSubSet)))
                call random_Matrix_SEED_10(exp1,array1, nrand, CEILING(MaxSubSet))
                ALLOCATE( new_comb1(nrand,size(ppp)))
                new_comb1 = (v1 + (v2 - v1) * array1)
                DEALLOCATE(array1)
                
                ALLOCATE(array2(nrand, CEILING(MaxSubSet)))
                call random_Matrix_SEED_10(exp1,array2, nrand, CEILING(MaxSubSet))
                ALLOCATE( new_comb2(nrand,size(ppp)))

                new_comb2 = ( v2 + (v3 - v2) * array2)              
                DEALLOCATE(array2)
            else
                random = RETURN_RANDOM_NUMBER(exp1)
                new_comb1 = (v1 + (v2 - v1) * random) 
                random = RETURN_RANDOM_NUMBER(exp1)
                new_comb2 = (v2 + (v3 - v2) * random) 
            end if

            ALLOCATE(new_comb(nrand,CEILING(MaxSubSet) * 2))
            new_comb(:,1:CEILING(MaxSubSet)) = new_comb1(:,1:CEILING(MaxSubSet))
            new_comb(:,(CEILING(MaxSubSet) + 1):CEILING(MaxSubSet) * 2) = new_comb2(:,1:CEILING(MaxSubSet))
                                                    
            if (ALLOCATED(new_comb1))  DEALLOCATE(new_comb1)
            if (ALLOCATED(new_comb2))  DEALLOCATE(new_comb2)
            
    END SUBROUTINE  generate_new_comb_matrix
    
            
 ! ----------------------------------------------------------------------
! SUBROUTINES generate_new_comb_matrix
! ----------------------------------------------------------------------     
    SUBROUTINE generate_new_comb_matrix_det(exp1, new_comb, ppp, MaxSubSet, nrand, v1,v2,v3, contador, matrixR)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: MaxSubSet
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ppp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  new_comb
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: v1, v2, v3
        INTEGER, INTENT(IN) :: nrand
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: array1, array2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: new_comb1, new_comb2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: random
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: matrixR
        INTEGER, INTENT(INOUT) :: contador        

            if (nrand > 1) then
                !print *, " - generate_new_comb_matrix -- new_comb = ", MaxSubSet * 2
                ALLOCATE(array1(nrand, CEILING(MaxSubSet)))
                call random_det_Matrix_SEED_10(array1, nrand, CEILING(MaxSubSet), contador, matrixR)
                ALLOCATE( new_comb1(nrand,size(ppp)))
                new_comb1 = (v1 + (v2 - v1) * array1)
                DEALLOCATE(array1)
                
                ALLOCATE(array2(nrand, CEILING(MaxSubSet)))
                call random_det_Matrix_SEED_10(array2, nrand, CEILING(MaxSubSet), contador, matrixR)
                ALLOCATE( new_comb2(nrand,size(ppp)))

                new_comb2 = ( v2 + (v3 - v2) * array2)              
                DEALLOCATE(array2)
            else
                CALL random_det_NUMBER(random, contador, matrixR)
                new_comb1 = (v1 + (v2 - v1) * random) 
                CALL random_det_NUMBER(random, contador, matrixR)
                new_comb2 = (v2 + (v3 - v2) * random) 
            end if

            ALLOCATE(new_comb(nrand,CEILING(MaxSubSet) * 2))
            new_comb(:,1:CEILING(MaxSubSet)) = new_comb1(:,1:CEILING(MaxSubSet))
            new_comb(:,(CEILING(MaxSubSet) + 1):CEILING(MaxSubSet) * 2) = new_comb2(:,1:CEILING(MaxSubSet))
                                                    
            if (ALLOCATED(new_comb1))  DEALLOCATE(new_comb1)
            if (ALLOCATED(new_comb2))  DEALLOCATE(new_comb2)
            
    END SUBROUTINE  generate_new_comb_matrix_det   
    
! ----------------------------------------------------------------------
! SUBROUTINES update_candidateset_with_new_comb
! Avaliase a soluci�ns xeradas, e se melloran a solucion pai, sustituese
! ----------------------------------------------------------------------      
    SUBROUTINE update_candidateset_with_new_comb( exp1,opts1, fitnessfunction,problem1,new_comb,candidateset,childset, &
            candidate_update, members_update,MaxSubSet2,nrand,nconst,nfuneval,counter, index1 )
            IMPLICIT NONE
            TYPE(problem), INTENT(IN) :: problem1
            TYPE(C_PTR), INTENT(INOUT) :: exp1
            TYPE(Refsettype), INTENT(INOUT) :: candidateset, childset
            TYPE(opts), INTENT(IN) :: opts1 
            TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj        
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: members_update, candidate_update
            INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
            INTEGER, INTENT(INOUT) :: counter
            INTEGER, INTENT(IN) :: nconst, nrand
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: MaxSubSet2
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  new_comb        
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: index1
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: starttime
            INTEGER(C_INT) :: mig, checkcooperativemigrationcriteriacessinner
            INTEGER :: i
            TYPE(outfuction) :: outfunct 
            INTEGER :: DIMEN(2)
                    
       
            counter = 1 
            DIMEN = SHAPE(new_comb)
            mig = 0
            do i = 1, DIMEN(2)
                if (mig .NE. 1) then
                ! Evaluate
                ALLOCATE(temp(nrand))
                temp = new_comb(:,i)

                outfunct = ssm_evalfc(exp1, fitnessfunction, temp, problem1, opts1, nconst, 0)
                nfuneval = nfuneval + 1
                DEALLOCATE(temp)
                if (outfunct%include .eq. 1) then
                    childset%x(:,counter)=outfunct%x
                    childset%f(counter)= outfunct%value
                    childset%fpen(counter) = outfunct%value_penalty
                    childset%penalty(counter) = outfunct%pena
                    if (ALLOCATED(outfunct%nlc) ) childset%nlc(:,counter) = outfunct%nlc
                    childset%parent_index(counter) = index1(i)
                    counter = counter + 1
                    if (outfunct%value_penalty .lt. candidateset%fpen(index1(i))) then
                            candidateset%x(:,index1(i)) = outfunct%x
                            candidateset%fpen(index1(i)) = outfunct%value_penalty
                            candidateset%f(index1(i)) = outfunct%value
                            candidateset%penalty(index1(i)) = outfunct%pena
                            if (ALLOCATED(outfunct%nlc)) candidateset%nlc(:,index1(i)) = outfunct%nlc
                            candidate_update(index1(i)) = 1
                            members_update(index1(i)) = 0
                    end if
                end if
                call destroy_outfuction(outfunct)
                end if
            end do

            CALL reajustchild(childset, counter, nrand, nconst)
            
    END SUBROUTINE update_candidateset_with_new_comb
            
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES check_stopping_criteria
! ----------------------------------------------------------------------        
    FUNCTION check_stopping_criteria (problem1, opts1, cputime, stopOptimization, fbest, nfuneval)  RESULT(fin)
        IMPLICIT NONE
        REAL(C_DOUBLE), INTENT(IN) :: cputime
        INTEGER, INTENT(IN) :: stopoptimization
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        TYPE(opts), INTENT(IN) :: opts1 
        TYPE(problem), INTENT(IN) :: problem1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: fbest
        INTEGER :: fin, i, j
        fin = 0
            ! Check the stop criterion
            if (nfuneval >= opts1%useroptions%maxeval) then
                fin = 1
            else if (cputime >= opts1%useroptions%maxtime) then
                fin = 2
            else if (stopOptimization .eq. 1) then
                fin = 4
            elseif (ALLOCATED(problem1%vtr)) then
                do i = 1, size(problem1%vtr)
                    do j = 1, size(fbest)
                        if (fbest(j) .LE. problem1%vtr(i)) then
                            fin = 3
                            EXIT
                        end if
                    end do
                end do
            end if
            
    END FUNCTION check_stopping_criteria
    
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES select_better_solution
! ----------------------------------------------------------------------      
    SUBROUTINE select_better_solution(results, opts1, refset, xbest, fbest, use_bestx, nfuneval, time, fin, exp1) 
        IMPLICIT NONE
        TYPE(opts), INTENT(IN) :: opts1 
        TYPE(resultsf), INTENT(INOUT) :: results
        TYPE(Refsettype), INTENT(INOUT) :: refset
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: xbest, fbest
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: time
        INTEGER :: clock_stop
        INTEGER, INTENT(IN) :: fin
        INTEGER, INTENT(INOUT) :: use_bestx
        INTEGER, DIMENSION(:), ALLOCATABLE :: index_result
        INTEGER, DIMENSION(:), ALLOCATABLE :: gg, hh
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE ::new_fbest,temp 
        INTEGER(C_LONG), INTENT(IN) :: nfuneval
        CHARACTER(len = 2) :: typesearch
        INTEGER :: auxsize
        TYPE(C_PTR), INTENT(INOUT) :: exp1

            ! Extraemos los indices de 
            ALLOCATE(index_result(size(refset%penalty)))
            typesearch = "LT"
            index_result = compare_vector(refset%penalty, opts1%useroptions%tolc, typesearch, auxsize)
            ALLOCATE(gg(auxsize))
            call index_of_ones(index_result, gg)
            CALL reajust_index(gg)
            DEALLOCATE(index_result)
            ! indice de los valores de penalty que son menores que la tolerancia de las restricciones
            if (ALLOCATED(gg)) then
                ALLOCATE(temp(size(gg)))
                temp(1: size(gg)) = refset%fpen(gg)
                if (.not. ALLOCATED(new_fbest)) ALLOCATE(new_fbest(1))
                if (.not. ALLOCATED(hh)) ALLOCATE(hh(1))                
                new_fbest = minval(temp)
                hh = minloc(temp)
                DEALLOCATE(temp)
                if (new_fbest(1) < fbest(1)) then
                    xbest=refset%x(:,gg(hh(1)))
                    fbest=new_fbest
                    use_bestx=1
#ifdef MPI2
                    CALL printsolution(exp1, xbest, fbest)
#endif
                    if (opts1%useroptions%inter_save .eq. 1) then  
                        !if (.not. ALLOCATED(results%fbest) ) ALLOCATE(results%fbest(size(fbest)))
                        !results%fbest(1)=fbest(1);
                        !if (.not. ALLOCATED(results%xbest) ) ALLOCATE(results%xbest(size(xbest)))
                        !results%xbest=xbest
                        !results%numeval=nfuneval
                        !results%end_crit=fin
                        !results%cputime=time
                        !results%refsett%x=refset%x                       
                        !results%refsett%f=refset%f                  
                        !results%refsett%fpen=refset%fpen  
                        !if (ALLOCATED(refset%nlc)) results%refsett%nlc=refset%nlc             
                        !results%refsett%penalty=refset%penalty  
                    end if
                end if
                if (ALLOCATED(new_fbest)) DEALLOCATE(new_fbest)
            end if
            IF (ALLOCATED(gg)) DEALLOCATE(gg)
            IF (ALLOCATED(hh)) DEALLOCATE(hh)
            
     END SUBROUTINE select_better_solution
 
     
     
     
! ----------------------------------------------------------------------
! SUBROUTINES update_refset_change
! ----------------------------------------------------------------------      
     SUBROUTINE update_refset_change(members_update, refset_change) 
         IMPLICIT NONE
         INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change
         INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: members_update
         INTEGER, DIMENSION(:), ALLOCATABLE :: members_update_non_zero_index, members_update_zero_index
         CHARACTER(len = 2) :: typesearch
         
         typesearch = "GE"
         ALLOCATE(members_update_non_zero_index(size(members_update)))
         members_update_non_zero_index = 0
         CALL find_in_vector_int(members_update, 1 , members_update_non_zero_index, typesearch)
         if (ALLOCATED(members_update_non_zero_index)) then
              CALL reajust_index(members_update_non_zero_index)
              if (ALLOCATED(members_update_non_zero_index)) &
                refset_change(members_update_non_zero_index) = refset_change(members_update_non_zero_index) + 1
         end if
            
         typesearch = "EQ"
         ALLOCATE(members_update_zero_index(size(members_update)))
         members_update_zero_index = members_update_zero_index * 0
         CALL find_in_vector_int(members_update, 0 , members_update_zero_index, typesearch)
         if (ALLOCATED(members_update_zero_index)) then
              CALL reajust_index(members_update_zero_index)  
              if (ALLOCATED(members_update_zero_index)) &
                refset_change(members_update_zero_index) = refset_change(members_update_zero_index) * 0
         end if
         
         
         if (ALLOCATED(members_update_non_zero_index))  DEALLOCATE(members_update_non_zero_index)
         if (ALLOCATED(members_update_zero_index))  DEALLOCATE(members_update_zero_index)         
            
    END SUBROUTINE update_refset_change
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES remove_possible_stuck_members_in_refset
! Esta funcion elimina solucion que estan estacadas e non melloran en X
! solucions
! ----------------------------------------------------------------------      
    SUBROUTINE remove_possible_stuck_members_in_refset(problem1, opts1, exp1, fitnessfunction, nconst, &
                        refset, refset_change,nfuneval, nvar, xl_log, xu_log )
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change
        TYPE(opts), INTENT(INOUT) :: opts1 
        TYPE(Refsettype), INTENT(INOUT) :: refset
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        INTEGER, INTENT(IN) :: nconst, nvar 
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        CHARACTER(len = 2) :: typesearch
        INTEGER, DIMENSION(:), ALLOCATABLE ::  index_members_to_change
        INTEGER :: to_replace, replaced, replaced_index
        TYPE(outfuction) :: outfunct
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: new_refset_member
        
        typesearch = "GT"
        ALLOCATE(index_members_to_change(size(refset_change)))
        index_members_to_change = 0

        
        CALL find_in_vector_int(refset_change, opts1 % useroptions % nstuck_solution, &
        index_members_to_change, typesearch)
        if (ALLOCATED(index_members_to_change)) then
            CALL reajust_index(index_members_to_change)
            if (ALLOCATED(index_members_to_change)) &
            refset_change(index_members_to_change) = 0
        end if
            
            if (ALLOCATED(index_members_to_change) ) then
                to_replace = size(index_members_to_change)
                replaced = 0
                replaced_index = 0
                do while (replaced < to_replace )
                        ALLOCATE(new_refset_member(nvar))
                        CALL random_Vector_SEED(exp1,new_refset_member, nvar, xl_log, xu_log)
                        if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0)  &
                                CALL converttonormal2(new_refset_member, nvar, opts1%useroptions%log_var)
                        outfunct = ssm_evalfc(exp1, fitnessfunction, new_refset_member, problem1, opts1, nconst, 0)
                        nfuneval = nfuneval + 1
                        if (outfunct%include .eq. 1) then
                            refset%x(:, index_members_to_change(replaced_index+1)) = outfunct%x
                            refset%f(index_members_to_change(replaced_index+1)) = outfunct%value
                            refset%fpen(index_members_to_change(replaced_index+1)) = outfunct%value_penalty
                            if (ALLOCATED(refset%nlc))  refset%nlc(:, index_members_to_change(replaced_index+1))=outfunct%nlc
                            refset%penalty(index_members_to_change(replaced_index+1)) = outfunct%pena
                            refset_change(index_members_to_change(replaced_index+1)) = 0;
                            replaced_index=replaced_index +1
                            replaced=replaced +1
                        end if
                        call destroy_outfuction(outfunct)
                        DEALLOCATE(new_refset_member)
                end do
                
            end if
            
            if (ALLOCATED(index_members_to_change))  DEALLOCATE(index_members_to_change)
     END SUBROUTINE remove_possible_stuck_members_in_refset       
     
     

     
     
     
! ----------------------------------------------------------------------
! SUBROUTINES remove_possible_stuck_members_in_refset
! ----------------------------------------------------------------------      
    SUBROUTINE remove_possible_stuck_members_in_refset_det(problem1, opts1, exp1, fitnessfunction, nconst, &
                        refset, refset_change,nfuneval, nvar, xl_log, xu_log, contador, matrixR )
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change
        TYPE(opts), INTENT(INOUT) :: opts1 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: matrixR
        INTEGER, INTENT(INOUT) :: contador
        TYPE(Refsettype), INTENT(INOUT) :: refset
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        INTEGER, INTENT(IN) :: nconst, nvar 
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        CHARACTER(len = 2) :: typesearch
        INTEGER, DIMENSION(:), ALLOCATABLE ::  index_members_to_change
        INTEGER :: to_replace, replaced,replaced_index
        TYPE(outfuction) :: outfunct
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: new_refset_member
        
        typesearch = "GT"
        ALLOCATE(index_members_to_change(size(refset_change)))
        index_members_to_change = 0

        
        CALL find_in_vector_int(refset_change, opts1 % useroptions % nstuck_solution, &
        index_members_to_change, typesearch)
        if (ALLOCATED(index_members_to_change)) then
            CALL reajust_index(index_members_to_change)
            if (ALLOCATED(index_members_to_change)) &
            refset_change(index_members_to_change) = 0
        end if
            
            if (ALLOCATED(index_members_to_change) ) then
                to_replace = size(index_members_to_change)
                replaced = 0
                replaced_index = 0
                 !print *, " - remove_possible_stuck_members_in_refset -- to_remplace= ", &
                 !   to_replace, "solutions"
                do while (replaced < to_replace )
                    !print *, "stuck"
                        ALLOCATE(new_refset_member(nvar))
                        
                        CALL random_det_Vector_SEED(new_refset_member, nvar, xl_log, xu_log, contador, matrixR)
                        if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0)  &
                                CALL converttonormal2(new_refset_member, nvar, opts1%useroptions%log_var)
                        outfunct = ssm_evalfc(exp1, fitnessfunction, new_refset_member, problem1, opts1, nconst, 0)
                        nfuneval = nfuneval + 1
                        if (outfunct%include .eq. 1) then
                            refset%x(:, index_members_to_change(replaced_index+1)) = outfunct%x
                            refset%f(index_members_to_change(replaced_index+1)) = outfunct%value
                            print *," Antes- ", refset%fpen(index_members_to_change(replaced_index+1)), " despois", &
                                            outfunct%value_penalty
                            refset%fpen(index_members_to_change(replaced_index+1)) = outfunct%value_penalty
                            if (ALLOCATED(refset%nlc))  refset%nlc(:, index_members_to_change(replaced_index+1))=outfunct%nlc
                            refset%penalty(index_members_to_change(replaced_index+1)) = outfunct%pena
                            refset_change(index_members_to_change(replaced_index+1)) = 0
                            replaced_index = replaced_index + 1
                        end if
                        call destroy_outfuction(outfunct)
                        DEALLOCATE(new_refset_member)
                        replaced=replaced +1
                end do
                
            end if
            
            if (ALLOCATED(index_members_to_change))  DEALLOCATE(index_members_to_change)
    END SUBROUTINE remove_possible_stuck_members_in_refset_det       
     
     
     
! ----------------------------------------------------------------------
! SUBROUTINES check_the_improvement_of_fbest
! ----------------------------------------------------------------------                               
     SUBROUTINE check_the_improvement_of_fbest(results, opts1, fbest, xbest, fbest_lastiter, nreset, cputime2, fin, nfuneval)
        IMPLICIT NONE
        TYPE(opts), INTENT(IN) :: opts1 
        TYPE(resultsf), INTENT(INOUT) :: results
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: fbest, fbest_lastiter,&
                                                                                                           xbest
        INTEGER, INTENT(INOUT) :: nreset, fin
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        REAL(C_DOUBLE), INTENT(IN) :: cputime2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE ::  temp
        REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: timeR
        INTEGER(C_LONG), DIMENSION(:), ALLOCATABLE :: index_result
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: auxiliarvalue
        INTEGER :: DIMEN(2)
            
        auxiliarvalue = REAL(0.0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
                
        if (fbest(1)<fbest_lastiter(1)) auxiliarvalue=1
                
        auxiliarvalue = auxiliarvalue/fbest_lastiter(1)
            
        if (auxiliarvalue > opts1%useroptions%tolc ) then
            results%f = fbest
            results%x(:,1) = xbest
            results%time(1) = cputime2
            results%neval(1) = nfuneval
            
               !ALLOCATE(temp(size(results%f) + 1))
               !temp(1:size(results%f)) = results%f
               !temp(   size(results%f)+1) = fbest(1)
               !CALL MOVE_ALLOC(temp, results%f )
                
               !DIMEN = shape(results%x)
               !ALLOCATE(temp2(DIMEN(1),DIMEN(2)+1))
               !temp2(:,1:DIMEN(2)) = results%x
               !temp2(:,DIMEN(2)+1) = xbest
               !CALL MOVE_ALLOC(temp2,results%x)
                
               !ALLOCATE(timeR(size(results%time) + 1))
               !timeR(1:size(results%time)) = results%time
               !timeR( size(results%time) +1) = cputime2
               !CALL MOVE_ALLOC(timeR, results%time )
               
               !ALLOCATE(index_result(size(results%neval) + 1))
               !index_result(1:size(results%neval)) = results%neval
               !index_result( size(results%neval)+1) = nfuneval
               !CALL MOVE_ALLOC(index_result, results%neval )
               
               nreset = 0
               
               !if opts1%useroptions%plot
         else 
                if (opts1%globaloptions%n_stuck > 0) then
                    nreset=nreset+1
                    if (nreset .eq. opts1%globaloptions%n_stuck) then
                        fin = 5
                    end if
                end if
         end if
     
    END SUBROUTINE check_the_improvement_of_fbest
          
    

! ----------------------------------------------------------------------
! SUBROUTINES apply_beyond_to_members_to_update
! ----------------------------------------------------------------------   
    SUBROUTINE apply_beyond_to_members_to_update( exp1, opts1, fitnessfunction, problem1, members_to_update,  nfuneval, &
        refset, candidateset, nconst, nrand, nvar)
        IMPLICIT NONE
        TYPE(opts), INTENT(IN) :: opts1 
        TYPE(Refsettype), INTENT(INOUT) :: refset, candidateset
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj 
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  ::  members_to_update
        INTEGER, INTENT(IN) :: nconst, nvar , nrand
        INTEGER :: i
        TYPE(outfuction) :: outfunct
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp1, temp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: auxiliarvalue
        
        if (ALLOCATED(members_to_update)) then
                do i = 1, size(members_to_update)
                    ALLOCATE(temp(nvar))
                    temp = refset % x(:,members_to_update(i))
                    ALLOCATE(temp1(nvar))
                    
                    temp1 = candidateset % x(:,members_to_update(i))
                    auxiliarvalue = candidateset % fpen(members_to_update(i))

                    ! ssm_beyond
                    CALL ssm_beyond(exp1, temp, temp1, auxiliarvalue, nfuneval, &
                          fitnessfunction, nrand, nconst, outfunct, opts1, problem1, 0)
                    DEALLOCATE(temp)
                    DEALLOCATE(temp1)

                    if (ALLOCATED(outfunct%x) ) then  
                        refset % x(:,members_to_update(i)) = outfunct % x
                        refset % f(members_to_update(i)) = outfunct % value
                        refset % fpen(members_to_update(i)) = outfunct % value_penalty
                        refset % penalty(members_to_update(i)) = outfunct % pena
                        if (nconst .GT. 0) then
                                refset % nlc(:,members_to_update(i)) = outfunct % nlc
                        end if
                    else    
                        refset % x(:,members_to_update(i)) = candidateset % x(:,members_to_update(i))
                        refset % f(members_to_update(i)) = candidateset % f(members_to_update(i))
                        refset % fpen(members_to_update(i)) = candidateset % fpen(members_to_update(i))
                        refset % penalty(members_to_update(i)) = candidateset % penalty(members_to_update(i))
                        if (nconst .GT. 0) then
                                refset%nlc(:,members_to_update(i)) = candidateset%nlc(:,members_to_update(i))
                        end if
                    end if
                    
                    call destroy_outfuction(outfunct)
                end do
            end if
            
    END SUBROUTINE apply_beyond_to_members_to_update
    
    
! ----------------------------------------------------------------------
! SUBROUTINES reajustchild
! ----------------------------------------------------------------------      
    SUBROUTINE reajustchild(childset, cont, nvar, nconst)
       IMPLICIT NONE
       TYPE(Refsettype), INTENT(INOUT) :: childset
       TYPE(Refsettype) :: childset_aux
       INTEGER, INTENT(IN) :: nvar, nconst, cont
       INTEGER :: counter 
       
       counter = cont-1
       
       if (counter .lt. size(childset%fpen)) then
        CALL create_childset(childset_aux, nvar, counter, nconst)
        childset_aux % x = childset % x(:, 1:counter)
        childset_aux % f = childset % f(1:counter)
        childset_aux % fpen = childset % fpen(1:counter)
        childset_aux % penalty = childset % penalty(1:counter)
        if (nconst .gt. 0) childset_aux % nlc = childset % nlc(:, 1:counter)
        childset_aux % parent_index = childset % parent_index(1:counter)
        CALL MOVE_ALLOC(childset_aux % x, childset % x)
        CALL MOVE_ALLOC(childset_aux % f, childset % f)
        CALL MOVE_ALLOC(childset_aux % fpen, childset % fpen)
        CALL MOVE_ALLOC(childset_aux % penalty, childset % penalty)
        CALL MOVE_ALLOC(childset_aux % parent_index, childset % parent_index)
        if (nconst .gt. 0) CALL MOVE_ALLOC(childset_aux % nlc, childset % nlc)

        end if
    END SUBROUTINE reajustchild



    SUBROUTINE check_duplicated_replace_DE (exp1, fitnessfunction,problem1,refset, opts1, dim_refset,&
                                                                      xl_log,xu_log,members_update, nfuneval, nvar, nconst )
        IMPLICIT NONE
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(Refsettype), INTENT(INOUT) :: refset
        TYPE(opts), INTENT(IN) :: opts1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction
        INTEGER, INTENT(INOUT) :: dim_refset
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: members_update
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER, INTENT(IN) :: nvar, nconst
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        INTEGER :: continuar, auxsize, DIMEN(2)
        INTEGER, DIMENSION(:), ALLOCATABLE :: indvect
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE :: denom22, AAA, temp2
        INTEGER :: BBB
        CHARACTER(len = 2) :: typesearch
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE :: new_refset_member
        INTEGER, DIMENSION(:), ALLOCATABLE :: index_result
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_result2
        TYPE(outfuction) :: outfunct
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE :: vi, vj
        INTEGER :: i,j
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: difference
        INTEGER :: counter_eval

        continuar = 0
        DIMEN = shape(refset%x)
        ALLOCATE(vi(DIMEN(1),1))
        ALLOCATE(vj(DIMEN(1),1))

        do i=1, dim_refset
          do j=1, dim_refset
             vi(:,1) = refset%x(:,i)
             if ( j .LT. i ) then
                vj(:,1) = refset%x(:,j)

                ALLOCATE(denom22(DIMEN(1),1))
                denom22 = max_matrix(abs_matrix(vi),abs_matrix(vj))

                call replace_matrix(denom22,REAL(1d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
                ALLOCATE(temp2(DIMEN(1),1))
                temp2 = (vi - vj)  / denom22

                ALLOCATE(AAA(DIMEN(1),1))
                AAA = abs_matrix(temp2)
                DEALLOCATE(temp2)
                typesearch = "LT"
                DIMEN = shape(AAA)
                auxsize = 0
                ALLOCATE(index_result2(DIMEN(1),1))

                !difference = ABS(refset%fpen(i) - refset%fpen(j))
                !if (difference .GE. 1.0d-6 ) then
                index_result2 = compare_matrix(AAA,REAL(1.0d-3,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), typesearch)
                !else
                ! index_result2 = compare_matrix(AAA,REAL(1.0d-1,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), typesearch)
                !end if
                ALLOCATE(index_result(1))
                CALL row_one(index_result2, index_result, auxsize)

                BBB = index_result(1)

                DEALLOCATE(index_result2)
                DEALLOCATE(index_result)

                counter_eval=0
                if (BBB .EQ. 1) then
                    outfunct%include = 0
                    do while ((outfunct%include .eq. 0) .AND. (counter_eval .LT. 3))
                        ALLOCATE(new_refset_member(nvar))
                        CALL random_Vector_SEED(exp1,new_refset_member,nvar,xl_log, xu_log)
                        if (ALLOCATED(opts1%useroptions%log_var) .and.size(opts1%useroptions%log_var) .GT. 0 )  &
                                CALL converttonormal2(new_refset_member,nvar,opts1%useroptions%log_var)

                        outfunct = ssm_evalfc( exp1, fitnessfunction,new_refset_member, problem1, opts1, nconst, 0)
                        nfuneval = nfuneval + 1
                        if (outfunct%include .eq. 1) then
                            refset%x(:,j) = outfunct%x
                            refset%f(j) = outfunct%value
                            refset%fpen(j) = outfunct%value_penalty
                            if (ALLOCATED(refset%nlc)) refset%nlc(:,j) = outfunct%nlc
                            refset%penalty(j) = outfunct%pena
                            members_update(j) = 0
                        else
                            counter_eval= counter_eval + 1
                        end if
                        DEALLOCATE(new_refset_member)
                    end do
                    call destroy_outfuction(outfunct)


                end if

                DEALLOCATE(denom22)
                if ( ALLOCATED(AAA)) DEALLOCATE(AAA)


              end if
          end do
        end do ! Fin bucle interior

            DEALLOCATE(vi)
            DEALLOCATE(vj)
    END SUBROUTINE check_duplicated_replace_de

! ----------------------------------------------------------------------
! SUBROUTINES serializerefset
! ----------------------------------------------------------------------      
    SUBROUTINE serializerefset(refset, nvar, tam, servector, rest)
        IMPLICIT NONE
        TYPE(Refsettype), INTENT(INOUT) :: refset
        INTEGER, INTENT(IN) :: nvar, tam      
        INTEGER, INTENT(INOUT) :: rest
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: servector
        INTEGER :: i, DIMEN(2)
        
        DO i=1, tam
            servector(1:nvar,i) =  refset%x(:,i)
            if (ALLOCATED(refset % nlc)) then
                DIMEN = shape(refset % nlc)
                servector((nvar + 1):(nvar + DIMEN(1)), i) = refset % nlc(:, i)
                servector(nvar + DIMEN(1) + 1, i) = refset % f(i)
                servector(nvar + DIMEN(1) + 2, i) = refset % penalty(i)
                servector(nvar + DIMEN(1) + 3, i) = refset % fpen(i)
                rest = nvar + 3
            else
                servector(nvar  + 1, i) = refset % f(i)
                servector(nvar  + 2, i) = refset % penalty(i)
                servector(nvar  + 3, i) = refset % fpen(i)    
                rest = 3
            end if

        END DO
            
            
    END SUBROUTINE    serializerefset
    
    
! ----------------------------------------------------------------------
! SUBROUTINES serializerefset
! ----------------------------------------------------------------------      
    SUBROUTINE deserializerefset(refset, nvar, tam, servector)
        IMPLICIT NONE
        TYPE(Refsettype), INTENT(INOUT) :: refset
        INTEGER, INTENT(IN) :: nvar, tam      
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: servector
        INTEGER :: i, DIMEN(2)
        
        DO i=1, tam
            refset%x(:,i) = servector(1:nvar,i) 
            if (ALLOCATED(refset%nlc)) then
                DIMEN=shape(refset%nlc)
                refset % nlc(:, i) =  servector((nvar + 1):(nvar + DIMEN(1)), i) 
                refset % f(i) =       servector(nvar + DIMEN(1) + 1, i)
                refset % penalty(i) = servector(nvar + DIMEN(1) + 2, i) 
                refset % fpen(i) =    servector(nvar + DIMEN(1) + 3, i) 
            else
                refset % f(i) = servector(nvar  + 1, i) 
                refset % penalty(i) = servector(nvar  + 2, i) 
                refset % fpen(i)   = servector(nvar  + 3, i)               
            end if

        END DO
            
            
    END SUBROUTINE    deserializerefset    
    
! ----------------------------------------------------------------------
! SUBROUTINES sizeseriallize
! ----------------------------------------------------------------------     
    SUBROUTINE sizeseriallize(nvar, nconst, size)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: nvar, nconst, size
       
        ! FORMAT: NCONST+NVARS+F+FPEN+PEN 
        if (nconst .GT. 0) then
            size = nconst+nvar+1+1+1
        else
            size = nvar+1+1+1
        endif
        
    END SUBROUTINE sizeseriallize
    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES checkbounds
! ----------------------------------------------------------------------  
    SUBROUTINE checkbounds(problem1,st)
        IMPLICIT NONE
        TYPE(problem), INTENT(IN) :: problem1
        INTEGER, INTENT(INOUT) :: st
        ! Check if bounds have the same dimension   CHECKBOUNDS
        if (size(problem1%XU) .NE. size(problem1%XL)) then
            PRINT *, "Upper and lower bounds have different dimension!!!!"
            PRINT *, "EXITING"
            CALL EXIT(st)
        end if
        ! ..   
    END SUBROUTINE checkbounds
    
! ----------------------------------------------------------------------
! SUBROUTINES initbestvars
! ----------------------------------------------------------------------  
    SUBROUTINE initbestvars(problem1,xbest,fbest,nvar) 
        IMPLICIT NONE
        TYPE(problem), INTENT(IN) :: problem1
        INTEGER, INTENT(IN) :: nvar
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: xbest, fbest
        INTEGER, DIMENSION(:), ALLOCATABLE :: iiim
        INTEGER :: iii
        
        !INIT BEST VARIABLES
        if (ALLOCATED(problem1%F0)) then
            if (ALLOCATED(fbest)) DEALLOCATE(fbest) 
            ALLOCATE(fbest(1))
            fbest = minval(problem1%F0)
            ALLOCATE(iiim(1))
            iiim = minloc(problem1%F0)
            iii = iiim(1)
            if (ALLOCATED(xbest)) DEALLOCATE(xbest) 
            ALLOCATE(xbest(nvar))
            xbest = problem1%X0(:,iii)
        else
            if (ALLOCATED(fbest)) DEALLOCATE(fbest) 
            ALLOCATE(fbest(1))
            fbest(1) = INFINITE
            if (ALLOCATED(xbest)) DEALLOCATE(xbest) 
            ALLOCATE(xbest(nvar))
            xbest = INFINITE
        end if
        
        IF (ALLOCATED(iiim)) DEALLOCATE(iiim)
    
    END SUBROUTINE initbestvars  
    
    
! ----------------------------------------------------------------------
! SUBROUTINES initintbinvars
! ----------------------------------------------------------------------  
    SUBROUTINE initintbinvars(problem1)
        IMPLICIT NONE
        TYPE(problem), INTENT(INOUT) :: problem1        
        
        if (problem1%int_var .EQ. - 1) then
            problem1%int_var = 0
        end if

        if (problem1%bin_var .EQ. - 1) then
            problem1%bin_var = 0
        end if
        ! ---
    END SUBROUTINE initintbinvars 
    
    
! ----------------------------------------------------------------------
! SUBROUTINES initlocaloptions
! ----------------------------------------------------------------------  
    SUBROUTINE initlocaloptions(opts1)    
        IMPLICIT NONE
        TYPE(opts), INTENT(INOUT) :: opts1         
            ! INIT LOCAL OPTIONS
        
        if (opts1 % localoptions % empty .eq. 0) then
           ! if (opts1 % localoptions % n1 .EQ. - 1) then
           !     opts1 % localoptions % n1 = 1
           ! end if

           ! if (opts1 % localoptions % n2 .EQ. - 1) then
           !     opts1 % localoptions % n2 = 10
           ! end if

            if (ALLOCATED(opts1 % localoptions % finish)) then
                opts1 % localoptions % finish(1) = opts1 % localoptions % solver
            end if
        end if
        
    END SUBROUTINE initlocaloptions 
        
    
! ----------------------------------------------------------------------
! SUBROUTINES initlocaloptions
! ----------------------------------------------------------------------  
    SUBROUTINE calcnconst(problem1,nconst) 
        IMPLICIT NONE
        TYPE(problem), INTENT(IN) :: problem1 
        INTEGER, INTENT(INOUT) :: nconst
        
        if (ALLOCATED(problem1%CU)) then
            nconst = size(problem1%CU)
        else
            nconst = 0
        end if
        
    END SUBROUTINE calcnconst 
    
    
! ----------------------------------------------------------------------
! SUBROUTINES check_output_obj_funct
! ----------------------------------------------------------------------     
    SUBROUTINE check_output_obj_funct(problem1, opts1,status, nconst)
        IMPLICIT NONE    
        TYPE(opts), INTENT(INOUT) :: opts1 
        TYPE(problem), INTENT(IN) :: problem1 
        INTEGER, INTENT(INOUT) :: status, nconst
                
        ! CHECK OUTPUT OBJECTIVE FUNCTION
        !Transform the objective function into a handle
        !fobj=str2func(problem.f);

        !Check if the objecttive function has 2 output arguments in constrained problems
        if (opts1 % localoptions % empty .eq. 0) then
            if (nconst .EQ. 1) then
                if (opts1 % localoptions % extrap % empty .eq. 0) then
                    PRINT *, "For constrained problems the objective function must have at least 2 output arguments \n"
                    PRINT *, "EXITING"
                    CALL EXIT(status)
                end if
            end if

            if (((problem1 % int_var + problem1 % bin_var) > 0) .AND. (opts1 % localoptions % solver .NE. '') &
                .AND. (opts1 % localoptions % solver .EQ. 'misqp')) then
                PRINT *, "For problems with integer and/or binary variables you must use MISQP as a local solver"
                PRINT *, "EXITING"
                CALL EXIT(status)
            end if

            ! --- COMPROBAMOS O ESTADO DOS PARÁMETROS DO LOCAL SOLVER N2FB, QUEDA PENDENTE
            if ((opts1 % localoptions % solver .EQ. 'nf2b') .OR. (opts1 % localoptions % solver .EQ. 'dnf2b')) then
                if (opts1 % localoptions % extrap % empty .eq. 0) then
                    PRINT *, "nf2d or dnf2b require 3 output arguments\n"
                    PRINT *, "EXITING"
                    CALL EXIT(status)
                else
                    !CALL random_Vector_SEED(randomV, nvar, problem1%XL, problem1%XU)
                end if
            end if
        end if
    
    END SUBROUTINE check_output_obj_funct
    
       
! ----------------------------------------------------------------------
! SUBROUTINES finalize_algorithm
! ----------------------------------------------------------------------       
    SUBROUTINE finalize_algorithm(exp1,problem1,results,opts1,fitnessfunction,refset,xbest,fbest,fin,nfuneval,NPROC, &
        cputime3, nconst, nvar,local_solver_var,idp)
        USE localsolver
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction       
        TYPE(resultsf), INTENT(INOUT) :: results
        TYPE(opts), INTENT(INOUT) :: opts1    
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: cputime3
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: xbest, fbest
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER(C_INT), INTENT(INOUT) :: fin, NPROC, nconst, nvar, idp
        TYPE(Refsettype), INTENT(INOUT) :: refset
        TYPE(problem), INTENT(INOUT) :: problem1 
        
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: X0
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: F0
        INTEGER(C_LONG) :: nevals
        INTEGER :: DIMEN(2)
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp2
        REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: timeR
        INTEGER(C_LONG), DIMENSION(:), ALLOCATABLE :: temp_int
        TYPE(outfuction) :: outfunct
        TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_var    

        
        if (fin .GT. 0) then
            if (opts1%localoptions%empty .eq. 0 ) then
               if (ALLOCATED(opts1%localoptions%finish) .and. (fin .LT. 4)) then
                    ALLOCATE(X0(nvar))
                    X0 = xbest
                    F0 = fbest(1)
                    
                    if  (opts1%useroptions%iterprint .eq. 1 ) then
                        WRITE (*,'(A10, I3, A33, A5)')    '[SLAVE ID=', idp,'] Final local refinement with  : ', &
                                opts1%localoptions%solver 
                        WRITE (*,'(A10, I3, A33, F25.10)')'(SLAVE ID=', idp,'] Initial point function value : ', F0
                    end if
                    nevals=0
                    CALL call_local_solver(problem1, exp1,opts1,fitnessfunction, X0, F0, NPROC,nevals,1,local_solver_var)
                    outfunct = ssm_evalfc(exp1,fitnessfunction, X0 ,problem1, opts1, nconst,0)
                    F0 = outfunct%value_penalty
                    nfuneval=nfuneval+nevals
                    if ( (outfunct%include .eq. 1 ) .and. (outfunct%pena .LE. opts1%useroptions%tolc ) ) then
                        if (  outfunct%value_penalty .LT. fbest(1) ) then
                            fbest(1) = F0
                            xbest=outfunct%value
                            if ( fbest(1) .LT. problem1%vtr(1)  ) then
                                fin = 3
                            end if
                        end if
                    end if
                
                    if ( opts1%useroptions%iterprint .eq. 1 ) then
                        WRITE (*,'(A10, I3, A33, F25.10)')'[SLAVE ID=', idp,'] Local solution function value: ', F0
                    end if
            
                end if
            end if
               
            ALLOCATE(temp(size(results%f) + 1))
            temp(1:size(results%f)) = results%f
            temp(   size(results%f)+1) = fbest(1)
            CALL MOVE_ALLOC(temp, results%f )
                
            DIMEN = shape(results%x)
            ALLOCATE(temp2(DIMEN(1),DIMEN(2)+1))
            temp2(:,1:DIMEN(2)) = results%x
            temp2(:,DIMEN(2)+1) = xbest
            CALL MOVE_ALLOC(temp2,results%x)
                
            ALLOCATE(timeR(size(results%time) + 1))
            timeR(1:size(results%time)) = results%time
            timeR( size(results%time) +1) = cputime3
            CALL MOVE_ALLOC(timeR, results%time )
               
            ALLOCATE(temp_int(size(results%neval) + 1))
            temp_int(1:size(results%neval)) = results%neval
            temp_int( size(results%neval)+1) = nfuneval
            CALL MOVE_ALLOC(temp_int, results%neval )
               
            if (.not. ALLOCATED(results%fbest) ) ALLOCATE(results%fbest(size(fbest)))
               
            results%fbest(1)=fbest(1)
            if (.not. ALLOCATED(results%xbest) ) ALLOCATE(results%xbest(size(xbest)))
            results%xbest = xbest
            results%numeval = nfuneval
            results%end_crit = fin
            results%cputime = cputime3
            results%refsett%x = refset%x
            results%refsett%f = refset%f
            results%refsett%fpen = refset%fpen
            if (ALLOCATED(results%refsett%nlc))    results%refsett%nlc =refset%nlc;                
            results%refsett%penalty = refset%penalty
            
        end if  
            
    END SUBROUTINE finalize_algorithm
    
    
! ----------------------------------------------------------------------
! SUBROUTINES calc_indexes
! ----------------------------------------------------------------------        
    SUBROUTINE calc_indexes(ncomb1,ncomb2,index1,index2,index,diff_index,auxsize)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ncomb1
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: ncomb2   
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: index1, index2,index, diff_index
        INTEGER, INTENT(INOUT) :: auxsize
        ! CALC INDEXES
        auxsize = choose(ncomb1)
        ALLOCATE(index1(auxsize))
        ALLOCATE(index2(auxsize))
        ALLOCATE(index(auxsize))
        ALLOCATE(diff_index(auxsize))
        
        index1 = ncomb2(:, 1)
        index2 = ncomb2(:, 2)
        index = ncomb2(:, 2)
        
        CALL FUSION_VECTOR_INT(index1, index)

        ! CALC PPP
        diff_index = (index2 - index1)
    
    END SUBROUTINE 
  
! ----------------------------------------------------------------------
! SUBROUTINES calc_ppp
! ----------------------------------------------------------------------    
    SUBROUTINE calc_ppp(auxsize,ppp,dim_refset,diff_index )
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: ppp
        INTEGER :: lb_p
        INTEGER, INTENT(INOUT) :: auxsize,dim_refset
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: st_p
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: diff_index
        lb_p = 1 !Minimum ppp
        st_p = REAL(0.75d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        !Defines maximum ppp. max(ppp)= lb_p+ st_p
        !Example: lb_p=1.25 and st_p=0.5 varies ppp from
        !%1.25 to 1.75

        ALLOCATE (ppp(auxsize))
        ppp = st_p * REAL((diff_index -1),KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/ &
           REAL((dim_refset - 2),KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) + &
           REAL(lb_p,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
    
    END SUBROUTINE calc_ppp
    
! ----------------------------------------------------------------------
! SUBROUTINES calc_hyper
! ----------------------------------------------------------------------    
    SUBROUTINE calc_hyper(nvar, hyper_x_L, hyper_x_U, MaxSubSet, XU, XL )    
        IMPLICIT NONE        
        ! CALC HYPER
        INTEGER, INTENT(INOUT):: nvar
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) ::MaxSubSet 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:),ALLOCATABLE, INTENT(INOUT):: hyper_x_L, hyper_x_U
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: XU, XL
        INTEGER :: i, j
        
        ALLOCATE(hyper_x_L(nvar,CEILING(MaxSubSet)))
        ALLOCATE(hyper_x_U(nvar,CEILING(MaxSubSet)))
        
        i=1
        j=1
        hyper_x_L = reshape((/ ((XL(i), i = 1,nvar), j = 1,  CEILING(MaxSubSet)) /),(/ nvar, CEILING(MaxSubSet) /))
        hyper_x_U = reshape((/ ((XU(i), i = 1,nvar), j = 1,  CEILING(MaxSubSet))/),(/ nvar, CEILING(MaxSubSet) /))

    END SUBROUTINE calc_hyper       



#ifdef MPI2
! ----------------------------------------------------------------------
! SUBROUTINES stoppingcriteriamodule
! ----------------------------------------------------------------------
    SUBROUTINE stoppingcriteriamodule(exp1,common_vars,opts1,results,fbest,xbest,time)
        IMPLICIT NONE
        TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
        TYPE(opts), INTENT(INOUT) :: opts1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(resultsf), INTENT(INOUT) :: results
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:),&
                        ALLOCATABLE, INTENT(INOUT) :: xbest,fbest
        REAL(C_DOUBLE) :: cputime1
        INTEGER :: dest, l, flag
        REAL(C_DOUBLE) :: calctimempi
        TYPE(time_ess), INTENT(INOUT) :: time
        INTEGER :: cooperativempitestess
                
        if (common_vars%fin .LT. 1) THEN
            do l = 1, common_vars%NPROC
                if (common_vars%fin .NE. 3) then
                    dest = l - 1
                    if (dest .NE. common_vars%idp) THEN
                        flag = cooperativempitestess(exp1, dest )
                        if (flag .EQ. 1) THEN
                           common_vars%fin = 3
                        end if
                    end if
                end if
            end do
        end if
        cputime1 = calctimeMPI(exp1,time%starttime)

        CALL asynchronousstoppingcriteriaess(exp1, common_vars%fin,common_vars%nfuneval,&
                fbest, cputime1,&
                opts1%useroptions%maxtime, common_vars%stopOptimization)
        CALL printcheckvtr(exp1, fbest(1), common_vars%nfuneval, cputime1, results%timevtr )

    END SUBROUTINE stoppingcriteriamodule
#endif



! ----------------------------------------------------------------------
! SUBROUTINES combine_scatter_search
! ----------------------------------------------------------------------
    SUBROUTINE combine_scatter_search(exp1,opts1,common_vars,ess_comm,pop,new_comb)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(populations), INTENT(INOUT) :: pop    
        TYPE(opts),INTENT(INOUT) :: opts1                                                                     
        TYPE(ess_common_vars), INTENT(INOUT) :: ess_comm
        TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:), ALLOCATABLE :: &
                                               factor, v1, v2, v3
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:), ALLOCATABLE, &
                                        INTENT(INOUT) :: new_comb
        INTEGER:: i, j
        !TYPE(problem), INTENT(INOUT) :: problem1
        !INTEGER :: i, j
        !INTEGER, INTENT(INOUT) :: counter
        
        ALLOCATE( factor(common_vars%nvar, size(ess_comm%ppp)))
        factor = reshape((/ ((ess_comm%ppp(j), i = 1, common_vars%nvar), j = 1,&
                                size(ess_comm%ppp)) /), (/common_vars%nvar,size(ess_comm%ppp) /))
        factor = factor * (pop%parents_index2 - pop%parents_index1)/ REAL(1.5d0,KIND = &
                                SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))

        ALLOCATE( v1(common_vars%nvar, size(ess_comm%ppp)))
        ALLOCATE( v2(common_vars%nvar, size(ess_comm%ppp)))
        ALLOCATE( v3(common_vars%nvar, size(ess_comm%ppp)))

        v1 = pop%parents_index1 - factor
        v2 = pop%parents_index2 - factor
        v3 = REAL(2.0d0,KIND =SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) &
                                        * pop%parents_index2  - pop%parents_index1 - factor

        CALL check_vector_in_hyper_bounds(exp1,opts1, v1, ess_comm%hyper_x_L,ess_comm%hyper_x_U )
        CALL check_vector_in_hyper_bounds(exp1,opts1, v3, ess_comm%hyper_x_L,ess_comm%hyper_x_U )
        CALL generate_new_comb_matrix(exp1, new_comb, ess_comm%ppp,ess_comm%MaxSubSet, &
                                common_vars%nrand, v1,v2,v3)

        if (ALLOCATED(factor))  DEALLOCATE(factor)
        if (ALLOCATED(v1))  DEALLOCATE(v1)
        if (ALLOCATED(v2))  DEALLOCATE(v2)
        if (ALLOCATED(v3))  DEALLOCATE(v3)

    END SUBROUTINE combine_scatter_search


    SUBROUTINE combine_differential_evolution(exp1,opts1,common_vars,pop,new_comb,bestx)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(populations), INTENT(INOUT) :: pop
        TYPE(opts),INTENT(INOUT) :: opts1
!        TYPE(ess_common_vars), INTENT(INOUT) :: ess_comm
        TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE, &
                                        INTENT(INOUT) :: new_comb
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE ::  bestx
        INTEGER :: i, j, point1, point2, point3, point4, point5
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::random

        ALLOCATE( new_comb(common_vars%nrand,opts1%globaloptions%dim_refset))

        !pop%parent_index1  
        DO i=1,opts1%globaloptions%dim_refset
                point1 = i
                point2 = i
                point3 = i
                point4 = i
                point5 = i
                do while (point1 .EQ. i)
                        point1 =CEILING( RETURN_RANDOM_NUMBER(exp1) * &
                               REAL(opts1%globaloptions%dim_refset,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
                end do
                do while ((point2 .EQ. i) .OR. (point2 .EQ. point1))
                        point2 =CEILING( RETURN_RANDOM_NUMBER(exp1) * &
                               REAL(opts1%globaloptions%dim_refset,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
                end do 
                do while ((point3 .EQ. i) .OR. (point3 .EQ. point2) .OR. (point3 .EQ. point1))
                        point3 =CEILING( RETURN_RANDOM_NUMBER(exp1) * &
                               REAL(opts1%globaloptions%dim_refset,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
                end do
                do while ((point4 .EQ. i) .OR. (point4 .EQ. point3) .OR. (point4 .EQ. point2) .OR. (point4 .EQ. point1))
                        point4 =CEILING( RETURN_RANDOM_NUMBER(exp1) * &
                               REAL(opts1%globaloptions%dim_refset,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
                end do
                do while ((point5 .EQ. i) .OR. (point5 .EQ. point4) .OR. (point5 .EQ. point3) .OR. &
                                                                        (point5 .EQ. point2) .OR. (point5 .EQ. point1))
                        point5 =CEILING( RETURN_RANDOM_NUMBER(exp1) * &
                               REAL(opts1%globaloptions%dim_refset,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)))
                end do

                DO j=1,common_vars%nvar
                   random = RETURN_RANDOM_NUMBER(exp1) 
                   if (random .LT. opts1%globaloptions%cr_de) then
                    if (opts1%globaloptions%mutation_de(1:6) .EQ. 'best/1') then
                      ! vi = x1 + F(x2 - x3 + x4 - x5)
                      new_comb(j,i) = bestx(j) + opts1%globaloptions%f_de * &
                             (pop%refset%x(j,point2) - pop%refset%x(j,point3))

                    end if
                    if (opts1%globaloptions%mutation_de(1:6) .EQ. 'best/2') then
                      ! vi = x1 + F(x2 - x3 + x4 - x5)
                      new_comb(j,i) = bestx(j) + opts1%globaloptions%f_de * &
                             (pop%refset%x(j,point2) + pop%refset%x(j,point3) - &
                               pop%refset%x(j,point4) - pop%refset%x(j,point5))
                           
                    end if
                    if (opts1%globaloptions%mutation_de(1:6) .EQ. 'rand/2') then
                      ! vi = best + F(x2 - x3 + x4 - x5)
                      new_comb(j,i) = pop%refset%x(j,point1) + opts1%globaloptions%f_de * &
                              (pop%refset%x(j,point2) + pop%refset%x(j,point3) - &
                              pop%refset%x(j,point4) - pop%refset%x(j,point5))
                    end if
                   end if
                END DO

        END DO
    END SUBROUTINE combine_differential_evolution






! ----------------------------------------------------------------------
! SUBROUTINES eval_solution
! ----------------------------------------------------------------------
    SUBROUTINE eval_solution(exp1,opts1,common_vars,SIZE_CHILDSET,pop,&
                                        fitnessfunction,problem1,openmp,counter,new_comb,index11)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(populations), INTENT(INOUT) :: pop
        TYPE(opts),INTENT(INOUT) :: opts1
        INTEGER, INTENT(IN) :: SIZE_CHILDSET
        TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
        INTEGER, INTENT(IN) :: openmp
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: index11
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE,&
                                                               INTENT(INOUT) :: new_comb
        INTEGER(C_INT) :: openmp_pos
        INTEGER(C_INT) :: getopenmpoption
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        TYPE(problem), INTENT(INOUT) :: problem1
        INTEGER, INTENT(INOUT) :: counter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: CHILDSIZE2

        if (openmp .EQ. 1) then  
               CHILDSIZE2 = REAL(SIZE_CHILDSET,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
               openmp_pos = getopenmpoption(exp1)
               if (openmp_pos .EQ. 1) then
#ifdef OPENMP
                       CALL update_candidateset_with_new_comb_parallel(exp1, opts1,&
                          fitnessfunction,problem1,new_comb,&
                          pop%candidateset,pop%childset,pop%candidate_update,pop%members_update,&
                          SIZE_CHILDSET,common_vars%nrand,common_vars%nconst,&
                          common_vars%nfuneval,counter,index11)
#endif
               else
                       CALL update_candidateset_with_new_comb(exp1,opts1,fitnessfunction,problem1,&
                                new_comb,pop%candidateset,pop%childset,pop%candidate_update,&
                                pop%members_update,CHILDSIZE2,common_vars%nrand,&
                                common_vars%nconst,common_vars%nfuneval,counter,index11 )

               end if
        else

              CALL update_candidateset_with_new_comb(exp1,opts1,fitnessfunction,problem1,new_comb,&
                                pop%candidateset,pop%childset,pop%candidate_update,pop%members_update,&
                                CHILDSIZE2,common_vars%nrand,common_vars%nconst,&
                                common_vars%nfuneval,counter, index11 )
        end if
                        
    END SUBROUTINE eval_solution



#ifdef OPENMP
! ----------------------------------------------------------------------
! SUBROUTINES evaluate_solutions_set_parallel
! ----------------------------------------------------------------------   
   SUBROUTINE evaluate_solutions_set_parallel(exp1,fitnessfunction,solutionset,problem1,opts1,solutions, nvar, nfuneval,nconst, &
                     xl_log, xu_log)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        TYPE(Refsettype), INTENT(INOUT) :: solutionset
        TYPE(Refsettype) :: solutionset_aux
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(opts), INTENT(IN) :: opts1
        INTEGER, INTENT(IN) :: nconst
        REAL(C_DOUBLE) ::  cputime1, cputime2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: solutions
        INTEGER, INTENT(IN)  :: nvar
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER :: l_f_0, l_x_0, DIMEN(2), tempnum, counter, counter2, i, sizet, neval
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: newsolution
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: entervalues
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: temp2, temp3
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: NAN, calctimeMPI
        INTEGER(C_INT) :: error, setNAN, continueeval
        INTEGER :: OMP_NUM_THREADS, IDPROC, clocal
        TYPE(outfuction), DIMENSION(:), ALLOCATABLE ::  outfunct1_VECTOR
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: xl_log, xu_log
        REAL :: SIZE_MIN
        INTEGER ::        counter_eval
        IDPROC = 0
        error = setNAN(NAN)
        
        if (ALLOCATED(problem1%F0)) then
            l_f_0 = size(problem1%F0)
        else
            l_f_0 = 0
        end if

        if (ALLOCATED(problem1%X0)) then
            DIMEN = shape(problem1%X0)
            l_x_0 = DIMEN(2)
        else
            l_x_0 = 0
        end if

        ! Put x_0 without their corresponding problem1%F0 values at the
        ! beginning of solutions
        !print *, "l_x_0", l_x_0, "-","l_f_0", l_f_0
        if ( (l_x_0  - l_f_0) > 0) then
            DIMEN = shape(problem1%X0)
            if (ALLOCATED(temp2)) DEALLOCATE(temp2)
            ALLOCATE(temp2(DIMEN(1), DIMEN(2)- l_f_0))
            temp2 = problem1%X0(:,(l_f_0 + 1):DIMEN(2))
            CALL FUSION_MATRIX(temp2, solutions)
            tempnum = l_x_0 - l_f_0
            !DEALLOCATE(temp2)
        else
            tempnum = 0
        end if
         
        sizet = opts1%globaloptions%ndiverse + 5 + tempnum

        counter = 0
        neval = 0 
        IDPROC = -1
        clocal = 0
        ALLOCATE(outfunct1_VECTOR(sizet))
        
        OMP_NUM_THREADS = omp_get_max_threads()
           
        SIZE_MIN = CEILING(REAL(opts1%globaloptions%dim_refset)/REAL(OMP_NUM_THREADS))
 !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) FIRSTPRIVATE(IDPROC,clocal) PRIVATE(newsolution,i,continueeval,counter_eval)
 !  REDUCTION(+:neval,counter) DEFAULT(SHARED) 
        do i = 1, sizet
            continueeval = 0
            if (IDPROC .EQ. -1 ) IDPROC = OMP_GET_THREAD_NUM()
            if ( .not. ALLOCATED(newsolution)) ALLOCATE(newsolution(nvar))
            newsolution = solutions(:,i)
            counter_eval=0
            do while ((continueeval .EQ. 0) .AND. (counter_eval .LT. 3))
                outfunct1_VECTOR(i) = ssm_evalfc(exp1, fitnessfunction, newsolution, problem1, opts1, nconst, IDPROC)
                neval = neval + 1

                if (outfunct1_VECTOR(i)%include .eq. 1) then
                        counter = counter + 1
                        clocal = clocal + 1
                        continueeval = 1
                else if ( (outfunct1_VECTOR(i)%include .eq. 0) .AND. &
                  ( REAL(clocal) .LT. SIZE_MIN)) then
                    counter_eval = counter_eval + 1
                    CALL random_Vector_SEED(exp1,newsolution, nvar, xl_log, xu_log)   
                    if (ALLOCATED(opts1%useroptions%log_var) .and. size(opts1%useroptions%log_var) .GT. 0 )  &
                           CALL converttonormal2(newsolution, nvar, opts1%useroptions%log_var)
                    call destroy_outfuction(outfunct1_VECTOR(i))
                else
                    continueeval = 1
                end if
            end do
            if (ALLOCATED(newsolution)) DEALLOCATE(newsolution)
        end do
 !$OMP END PARALLEL DO
        
        ALLOCATE(solutionset%x(nvar, counter))
        solutionset%x =  REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        if (nconst .GT. 0 ) then 
            ALLOCATE(solutionset%nlc(nconst , counter))
            solutionset%nlc = REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) 
        end if
        ALLOCATE(solutionset%f(counter))
        solutionset%f =  REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(solutionset%fpen(counter))
        solutionset%fpen = REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(solutionset%penalty(counter))
        solutionset%penalty = REAL(0.0d0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        
        counter2 = 1
        do i = 1, counter
            if (outfunct1_VECTOR(i)%include .eq. 1) then
                solutionset%x(:,counter2) = outfunct1_VECTOR(i)%x
                solutionset%f(counter2) = outfunct1_VECTOR(i)%value
                solutionset%fpen(counter2) = outfunct1_VECTOR(i)%value_penalty
                solutionset%penalty(counter2) = outfunct1_VECTOR(i)%pena
                if (nconst .GT. 0) then 
                    solutionset%nlc(:,counter2) = outfunct1_VECTOR(i)%nlc
                end if
                counter2 = counter2 + 1
            end if
            call destroy_outfuction(outfunct1_VECTOR(i))
        end do
        DEALLOCATE(outfunct1_VECTOR)
        nfuneval = nfuneval + neval

        !resize solutionset
        counter = counter2-1
        if ( counter .LT. sizet )  then
                sizet = counter
                ALLOCATE(solutionset_aux%x(nvar, sizet))
                solutionset_aux%x(:,1:counter)=solutionset%x(:,1:counter)
                if (nconst .GT. 0 ) then
                        ALLOCATE(solutionset_aux%nlc(nconst , sizet))
                        solutionset_aux%nlc(:,1:counter)=solutionset%nlc(:,1:counter)
                end if
                ALLOCATE(solutionset_aux%f(sizet))
                solutionset_aux%f(1:counter)=solutionset%f(1:counter)
                ALLOCATE(solutionset_aux%fpen(sizet))
                solutionset_aux%fpen(1:counter)=solutionset%fpen(1:counter)
                ALLOCATE(solutionset_aux%penalty(sizet))
                solutionset_aux%penalty(1:counter)=solutionset%penalty(1:counter)

                DEALLOCATE(solutionset%x)
                ALLOCATE(solutionset%x(nvar, sizet))
                if (nconst .GT. 0 ) then
                        DEALLOCATE(solutionset%nlc)
                        ALLOCATE(solutionset%nlc(nconst , sizet))
                end if
                DEALLOCATE(solutionset%f)
                ALLOCATE(solutionset%f(sizet))
                DEALLOCATE(solutionset%fpen)
                ALLOCATE(solutionset%fpen(sizet))
                DEALLOCATE(solutionset%penalty)
                ALLOCATE(solutionset%penalty(sizet))

                solutionset%x = solutionset_aux%x
                if (nconst .GT. 0 ) then
                        solutionset%nlc=solutionset_aux%nlc
                end if
                solutionset%f=solutionset_aux%f
                solutionset%fpen=solutionset_aux%fpen
                solutionset%penalty=solutionset_aux%penalty

                DEALLOCATE(solutionset_aux%x)
                if (nconst .GT. 0 ) then
                        DEALLOCATE(solutionset_aux%nlc)
                end if
                DEALLOCATE(solutionset_aux%f)
                DEALLOCATE(solutionset_aux%fpen)
                DEALLOCATE(solutionset_aux%penalty)

        end if
            

        !Add points with f_0 values. We assume feasiblity
        if (l_f_0 .GT. 0) then
            DIMEN = shape(problem1%X0)
            if (ALLOCATED(temp2)) DEALLOCATE(temp2)
            ALLOCATE(temp2(DIMEN(1),l_f_0))
            temp2 = problem1%X0(:,1:l_f_0)
            call FUSION_MATRIX(temp2, solutionset%x)
            call FUSION_VECTOR(problem1%F0, solutionset%f)
            call FUSION_VECTOR(problem1%F0, solutionset%fpen)
            ALLOCATE(entervalues(l_f_0))
            entervalues = entervalues * 0d0
            call FUSION_VECTOR(entervalues, solutionset%penalty)
            DEALLOCATE(entervalues)
            
            if (ALLOCATED(solutionset%nlc)) then
                ALLOCATE(temp3(nconst,l_f_0))
                temp3=0
                CALL FUSION_MATRIX(temp3,solutionset%nlc)
            end if
        end if    
        
        
   END SUBROUTINE evaluate_solutions_set_parallel  
#endif    


#ifdef OPENMP    
! ----------------------------------------------------------------------
! SUBROUTINES apply_beyond_to_members_to_update_parallel
! ----------------------------------------------------------------------   
    SUBROUTINE apply_beyond_to_members_to_update_parallel( exp1, opts1, fitnessfunction, problem1, members_to_update,  nfuneval, &
        refset, candidateset, nconst, nrand, nvar, timeparallel )
        IMPLICIT NONE
        TYPE(opts), INTENT(IN) :: opts1 
        TYPE(Refsettype), INTENT(INOUT) :: refset, candidateset
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj 
        INTEGER(C_LONG), INTENT(INOUT) :: nfuneval
        INTEGER :: clock_rate, clock_start, clock_stop        
        INTEGER(C_LONG) :: neval
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  ::  members_to_update
        INTEGER, INTENT(IN) :: nconst, nvar , nrand
        REAL(C_DOUBLE), INTENT(INOUT) :: timeparallel
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp1, temp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: auxiliarvalue
        INTEGER :: OMP_NUM_THREADS, IDPROC
        TYPE(outfuction), DIMENSION(:), ALLOCATABLE :: outfunct1
        INTEGER :: i 
        IDPROC = 0
!$OMP PARALLEL SHARED(OMP_NUM_THREADS)           
           OMP_NUM_THREADS = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL    
           

               
       ALLOCATE(outfunct1(OMP_NUM_THREADS)) 
       ALLOCATE(temp(nvar))
       ALLOCATE(temp1(nvar))     

            
       if (ALLOCATED(members_to_update)) then
                ! exploramos los miembros a actualizar (los que prometen más resultados)
                ! y aplicamos ls técnica go-beyond en ellos
            
 !$OMP PARALLEL PRIVATE(OMP_NUM_THREADS,IDPROC) DEFAULT(SHARED)
                IDPROC = OMP_GET_THREAD_NUM()
                neval = 0
 !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(temp,temp1,i,auxiliarvalue) REDUCTION(+:neval)
                do i = 1, size(members_to_update)

                    temp = refset%x(:,members_to_update(i))
                    temp1 = candidateset%x(:,members_to_update(i))
                    auxiliarvalue = candidateset%fpen(members_to_update(i))

                    ! ssm_beyond
                    CALL ssm_beyond(exp1, temp, temp1, auxiliarvalue, neval, &
                    fitnessfunction, nrand, nconst, outfunct1(IDPROC+1), opts1, problem1, IDPROC)
                          
                             
                    if (ALLOCATED(outfunct1(IDPROC+1)%x) ) then  
                        refset%x(:,members_to_update(i)) = outfunct1(IDPROC+1)%x
                        refset%f(members_to_update(i)) = outfunct1(IDPROC+1)%value
                        refset%fpen(members_to_update(i)) = outfunct1(IDPROC+1)%value_penalty
                        refset%penalty(members_to_update(i)) = outfunct1(IDPROC+1)%pena
                        if (nconst .GT. 0) refset%nlc(:,members_to_update(i)) = outfunct1(IDPROC+1) % nlc
                        !end if
                    else    
                        refset%x(:,members_to_update(i)) = candidateset%x(:,members_to_update(i))
                        refset%f(members_to_update(i)) = candidateset%f(members_to_update(i))
                        refset%fpen(members_to_update(i)) = candidateset%fpen(members_to_update(i))
                        refset%penalty(members_to_update(i)) = candidateset%penalty(members_to_update(i))
                        if (nconst .GT. 0) refset%nlc(:,members_to_update(i)) = candidateset%nlc(:,members_to_update(i))
                    end if
                    
                end do
 !$OMP END DO

 !$OMP END PARALLEL   
        
       nfuneval = nfuneval + neval
       end if
            
            DEALLOCATE(temp)
            DEALLOCATE(temp1)
            call destroy_outfuction(outfunct1(IDPROC+1))
            if (ALLOCATED(outfunct1)) DEALLOCATE(outfunct1)    
            
    END SUBROUTINE apply_beyond_to_members_to_update_parallel
#endif    
    
#ifdef OPENMP 
#ifdef MPI2
! ----------------------------------------------------------------------
! SUBROUTINES update_candidateset_with_new_comb_parallel
! ----------------------------------------------------------------------   
    SUBROUTINE update_candidateset_with_new_comb_parallel( exp1,opts1, fitnessfunction,problem1,new_comb,candidateset,childset, &
            candidate_update, members_update,MaxSubSet2,nvar,nconst,nfuneval,counter, index3 )
            IMPLICIT NONE
            TYPE(problem), INTENT(IN) :: problem1
            TYPE(C_PTR),  INTENT(INOUT) :: exp1
            TYPE(Refsettype), INTENT(INOUT) :: candidateset, childset
            TYPE(opts), INTENT(IN) :: opts1 
            TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj        
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: members_update, candidate_update
            INTEGER(C_LONG), INTENT(INOUT) ::  nfuneval
            INTEGER, INTENT(INOUT) :: counter
            INTEGER :: neval    
            INTEGER, INTENT(IN) :: nconst, nvar
            INTEGER, INTENT(IN) :: MaxSubSet2
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  new_comb        
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: index3
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: newsolution
            INTEGER :: i, OMP_NUM_THREADS, IDPROC, ii
            TYPE(outfuction), DIMENSION(:), ALLOCATABLE :: outfunct1, outfunct1_VECTOR
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: cputime1,cputime2,calctimeMPI
            INTEGER :: DIMEN(2) 
            IDPROC = 0
            
       
            
            OMP_NUM_THREADS = omp_get_max_threads()

            ALLOCATE(outfunct1(OMP_NUM_THREADS)) 
            ALLOCATE(outfunct1_VECTOR((MaxSubSet2)))

            neval = 0
            DIMEN=shape(new_comb)


 !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(newsolution,i,IDPROC)
 !REDUCTION(+:neval)
            do i = 1, DIMEN(2)
                ! Evaluate
                IDPROC = OMP_GET_THREAD_NUM()
                ALLOCATE(newsolution(nvar))
                newsolution = new_comb(:,i)
                
                
                outfunct1(IDPROC+1) = ssm_evalfc(exp1, fitnessfunction, newsolution, problem1, opts1, nconst,IDPROC)
                neval = neval + 1 
                outfunct1_VECTOR(i) = outfunct1(IDPROC+1)

                call destroy_outfuction(outfunct1(IDPROC+1))
                if (ALLOCATED(newsolution))  DEALLOCATE(newsolution)                
            end do
 !$OMP END PARALLEL DO
            
            if (ALLOCATED(outfunct1)) DEALLOCATE(outfunct1)
            counter = 1
            do i = 1, DIMEN(2)
                if (outfunct1_VECTOR(i)%include .eq. 1) then
                    childset%x(:,counter)=outfunct1_VECTOR(i)%x
                    childset%f(counter)= outfunct1_VECTOR(i)%value
                    childset%fpen(counter) = outfunct1_VECTOR(i)%value_penalty
                    childset%penalty(counter) = outfunct1_VECTOR(i)%pena
                    if (ALLOCATED(outfunct1_VECTOR(i)%nlc) ) childset%nlc(:,counter) = outfunct1_VECTOR(i)%nlc
                    childset%parent_index(counter) = index3(i)
                    counter = counter + 1
                   if (outfunct1_VECTOR(i)%value_penalty .LT. candidateset%fpen(index3(i))) then
                            candidateset%x(:,index3(i)) = outfunct1_VECTOR(i)%x
                            candidateset%fpen(index3(i)) = outfunct1_VECTOR(i)%value_penalty
                            candidateset%f(index3(i)) = outfunct1_VECTOR(i)%value
                            candidateset%penalty(index3(i)) = outfunct1_VECTOR(i)%pena
                            if (ALLOCATED(outfunct1_VECTOR(i)%nlc)) candidateset%nlc(:,index3(i)) = outfunct1_VECTOR(i)%nlc
                            candidate_update(index3(i)) = 1
                            members_update(index3(i)) = 0
                    end if

                end if
                call destroy_outfuction(outfunct1_VECTOR(i))              
            end do            
            
            CALL reajustchild(childset, counter, nvar, nconst)
            if (ALLOCATED(outfunct1_VECTOR)) DEALLOCATE(outfunct1_VECTOR)
            

            nfuneval = nfuneval + neval
            
            
                
    END SUBROUTINE update_candidateset_with_new_comb_parallel   
#endif
#endif

#ifdef MPI2
    SUBROUTINE local_solver_method1(exp1,opts1,common_vars,local_solver_var,time,pop,&
                                    problem1,results,fitnessfunction,xbest,fbest,fbest_lastiter) 
      IMPLICIT NONE
#ifdef MPI2
      REAL(C_DOUBLE) :: calctimempi
#endif
      REAL(C_DOUBLE) :: cputime1
      TYPE(algorithm_common_vars),INTENT(INOUT) :: common_vars
      TYPE(problem),INTENT(INOUT) :: problem1
      TYPE(local_solver_help_vars),INTENT(INOUT) :: local_solver_var
      TYPE(opts), INTENT(INOUT) :: opts1  
      TYPE(time_ess), INTENT(INOUT) :: time
      TYPE(C_PTR), INTENT(INOUT) :: exp1
      TYPE(populations), INTENT(INOUT) :: pop    
      TYPE(resultsf), INTENT(INOUT) :: results
      TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
      REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE&
                                                  ,INTENT(INOUT) :: xbest,fbest,fbest_lastiter

      if (( opts1%localoptions%empty .eq. 0 )  .and. (common_vars%fin .eq. 0)) then
        if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. (opts1%localoptions%bestx .gt. 0) .and. &
                   (local_solver_var%n_minimo .ge. local_solver_var%n_critico) ) then
            if ( local_solver_var%use_bestx .gt. 0 ) then
                cputime1 = calctimeMPI(exp1,time%starttime)
                CALL printgant(exp1,cputime1,1)
                CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                            pop%childset,xbest,fbest,common_vars%nfuneval,fbest_lastiter,pop%refset,pop%refset_change,1)
                cputime1 = calctimeMPI(exp1,time%starttime)
                CALL printgant(exp1,cputime1,2)
                local_solver_var%n_minimo=0
            end if
        else
            if ( ( LEN(opts1%localoptions%solver) .gt. 0 ) .and. (local_solver_var%n_minimo .ge. &
                               local_solver_var%n_critico) ) then
                cputime1 = calctimeMPI(exp1,time%starttime)
                CALL printgant(exp1,cputime1,1)
                CALL ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                            pop%childset,xbest,fbest,common_vars%nfuneval,fbest_lastiter,pop%refset,pop%refset_change,1)
                cputime1 = calctimeMPI(exp1,time%starttime)
                CALL printgant(exp1,cputime1,2)
                local_solver_var%n_minimo=0
            end if
        end if
     end if

   END SUBROUTINE local_solver_method1
#endif

   SUBROUTINE calc_init_solt_and_refset_global(opts1,exp1,pop,common_vars,xl_log,xu_log,&
                                                                   problem1,fitnessfunction,openth)
        IMPLICIT NONE
        TYPE(opts),  INTENT(INOUT) :: opts1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(populations), INTENT(INOUT) :: pop
        TYPE(algorithm_common_vars), INTENT(INOUT) :: common_vars
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:), ALLOCATABLE,&
                                         INTENT(INOUT):: xl_log, xu_log
        TYPE(problem), INTENT(INOUT) :: problem1
        INTEGER, INTENT(IN) :: openth
        INTEGER(C_INT) :: getopenmpoption, openmp_pos
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
                
        if (common_vars%idp .NE. 0) then
            CALL create_init_solutions(exp1,opts1,pop%solutions,common_vars%nvar,xl_log,xu_log)
            if ((problem1%int_var .GT. 0) .OR. (problem1%bin_var .GT. 1)) then
                CALL ssm_round_int( pop%solutions, problem1%int_var + problem1%bin_var, &
                         problem1%XL,problem1%XU)
            end if
            if ( openth .EQ. 1 ) then
               openmp_pos = getopenmpoption(exp1)
               if (openmp_pos .EQ. 1) then
#ifdef OPENMP
                 CALL evaluate_solutions_set_parallel(exp1,fitnessfunction,&
                        pop%solutionset, problem1, opts1, pop%solutions, &
                        common_vars%nvar, common_vars%nfuneval, &
                        common_vars%nconst, xl_log, xu_log)
#endif
               else
                 CALL evaluate_solutions_set(exp1,fitnessfunction, pop%solutionset,problem1,&
                                 opts1,pop%solutions,common_vars%nvar,common_vars%nfuneval, &
                                 common_vars%nconst, xl_log, xu_log)
               end if
            else
               CALL evaluate_solutions_set(exp1,fitnessfunction,pop%solutionset,problem1,opts1, &
                       pop%solutions, common_vars%nvar,common_vars%nfuneval, common_vars%nconst,&
                       xl_log, xu_log)
            end if

            CALL create_refset (exp1, pop%solutionset, pop%refset,common_vars%nvar, opts1, common_vars%nconst, &
                opts1%globaloptions%dim_refset )

        else
            CALL create_refset_empty( pop%refset, common_vars%nvar, opts1,common_vars%nconst, &
                                                                    opts1%globaloptions%dim_refset)
        end if
        
    END SUBROUTINE calc_init_solt_and_refset_global



END MODULE scattersearchfunctions

