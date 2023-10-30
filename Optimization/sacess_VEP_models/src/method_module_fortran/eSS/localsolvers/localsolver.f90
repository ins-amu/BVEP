#define TIMELOCAL 1

MODULE localsolver
    USE iso_c_binding
    USE scattersearchtypes
    USE common_functions
    USE misqp_interface
    USE qsort_module
    USE funcevalinterface
    USE dhc_mod
CONTAINS
! ----------------------------------------------------------------------
! SUBROUTINES build_aux_local_search
! ----------------------------------------------------------------------     
    SUBROUTINE build_aux_local_search(problem1,local_solver_vars, opts1, exp1)
        IMPLICIT NONE    
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(problem), INTENT(INOUT) :: problem1
        TYPE(opts), INTENT(INOUT) :: opts1
        TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_vars   
        INTEGER :: evals, iter    
  
        CALL chargelsevalmax(exp1,local_solver_vars%evals_per_iter)

        ! Build auxiliary files in case local search is activated
        if (opts1%localoptions%empty .eq. 0) then
            if ( opts1%localoptions%solver(1:5) .EQ. 'misqp' ) then
               CALL ssm_aux_local_misqp (problem1,local_solver_vars)
            end if
        end if
       
        !Initialize variables local solver
        local_solver_vars%stage_1=1
        local_solver_vars%pass_misqp = 0
        local_solver_vars%stage_2=0
        local_solver_vars%n_minimo = 0;
        local_solver_vars%use_bestx = 0        
        if (opts1%localoptions%empty .eq. 0) then 
            local_solver_vars%n_critico = opts1%localoptions%n1
        end if
        
    END SUBROUTINE build_aux_local_search    

    SUBROUTINE  call_dhc( problem1,exp1,opts1,fitnessfunction, x0, fval, numeval, local_solver_vars )
        IMPLICIT NONE
        ! fobj - fitnessfunction
        ! weight - opts%user%weight
        ! iterprint - opts%local%iterprint
        ! tolc - opts%user%tolc
        TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_vars
        TYPE(problem), INTENT(INOUT) :: problem1
        TYPE(opts), INTENT(INOUT) :: opts1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        INTEGER(C_LONG), INTENT(INOUT) :: numeval
        INTEGER(C_LONG) :: eval
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) ::fval
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: x0        
        INTEGER :: budget
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: initsize
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: thres
        
        X0 = (X0 - problem1%XL)/(problem1%XU-problem1%XL)
        initsize = 1d-1
        budget = local_solver_vars%evals_per_iter

        if (opts1%localoptions%tol .EQ. 1) then
            thres = 1d-6
        else if (opts1%localoptions%tol .EQ. 2) then
            thres = 1d-8
        else if (opts1%localoptions%tol .EQ. 3) then    
            thres = 1d-10
        end if
        
        CALL dhc(problem1,exp1,opts1,fitnessfunction,x0,initsize,thres,budget, eval, fval)

        numeval = numeval + eval
    END SUBROUTINE call_dhc
    
 
    SUBROUTINE  call_misqp( problem1,exp1,opts1,fitnessfunction, x0, fval, numeval, exitflag, local_solver_vars )
        IMPLICIT NONE
        ! fobj - fitnessfunction
        ! weight - opts%user%weight
        ! iterprint - opts%local%iterprint
        ! tolc - opts%user%tolc
        TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_vars    
        TYPE(problem), INTENT(INOUT) :: problem1
        TYPE(opts), INTENT(INOUT) :: opts1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        INTEGER(C_LONG) :: nfunc
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) ::fval
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: x0 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: g
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: acc
        INTEGER(C_LONG), INTENT(INOUT) :: numeval
        INTEGER, INTENT(INOUT) :: exitflag        
        INTEGER :: ifail
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::  NAN, INF     
        INTEGER, DIMENSION(:), ALLOCATABLE :: RESULTG
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: valuecomp
        CHARACTER(len = 2) :: typesearch        
        !ncont = size(X0)- problem1%int_var - problem1%bin_var
        !n = ncont + problem1%int_var + problem1%bin_var
        !m = size(problem1%XU) + size(problem1%XL) + problem1%neq
        CALL setnan(NAN)
        CALL setdblmax(INF)        

        if (opts1%localoptions%tol .EQ. 1) then
                     acc = opts1%useroptions%tolc*100
        else if (opts1%localoptions%tol .EQ. 2) then
                     acc = opts1%useroptions%tolc
        else if (opts1%localoptions%tol .EQ. 3) then
                    acc = opts1%useroptions%tolc/100
        end if

        
        CALL RUN_MISQP( problem1,exp1,opts1,fitnessfunction,acc, x0, fval, nfunc, ifail, g, local_solver_vars  )
        
        
        if (ALLOCATED(g)) then
            typesearch = "GE"
            valuecomp = -1d0 * opts1%useroptions%tolc        
            CALL find_in_vector(g, valuecomp, RESULTG, typesearch)  
        
            if ((fval .NE. NAN) .AND. (fval .NE. INF ) .AND. (ALLOCATED(RESULTG))) then
                exitflag=1
            else
                exitflag=-1
            end if
        else
            if ((fval .NE. NAN) .AND. (fval .NE. INF ) ) then
                exitflag=1
            else
                exitflag=-1
            end if
        end if
        
        numeval = numeval + nfunc
        
    END SUBROUTINE call_misqp   

! ----------------------------------------------------------------------
! SUBROUTINE call_local_solver
! ----------------------------------------------------------------------
   SUBROUTINE call_local_solver(problem1, exp1,opts1,fitnessfunction, x0, fval,NPROC,neval, centercriteria,&
                                    local_solver_vars)
       IMPLICIT NONE
       TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_vars    
       TYPE(problem), INTENT(INOUT) :: problem1
       TYPE(opts), INTENT(INOUT) :: opts1
       TYPE(C_PTR), INTENT(INOUT) :: exp1
       INTEGER, INTENT(INOUT) ::  NPROC
       INTEGER, INTENT(IN) :: centercriteria
       INTEGER(C_LONG), INTENT(INOUT) :: neval
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) ::fval
       TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: x0
       INTEGER :: exitflag, opt
       INTEGER(C_LONG) :: iter
     
  
       if ( opts1%localoptions%solver(1:12) .EQ. 'nl2sol.dn2gb' ) then          
            iter = 100
            CALL callnl2sol(exp1,fitnessfunction,x0,neval,local_solver_vars%evals_per_iter,iter, NPROC, 0)
       else if ( opts1%localoptions%solver(1:12) .EQ. 'nl2sol.dn2fb' ) then
            iter = 100
            CALL callnl2sol(exp1,fitnessfunction,x0,neval,local_solver_vars%evals_per_iter,iter, NPROC, 1)
       else if ( opts1%localoptions%solver(1:5) .EQ. 'misqp' ) then
            CALL call_misqp( problem1,exp1,opts1,fitnessfunction, x0, fval, neval, exitflag, local_solver_vars )
       else if ( opts1%localoptions%solver(1:3) .EQ. 'dhc' ) then
            CALL call_dhc( problem1,exp1,opts1,fitnessfunction, x0, fval, neval, local_solver_vars )
       end if
       
   END SUBROUTINE call_local_solver

! ----------------------------------------------------------------------
! SUBROUTINE transform_norm
! ----------------------------------------------------------------------
   SUBROUTINE transform_norm ( oldmt, newmt, UB, LB)
       IMPLICIT NONE
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: oldmt
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  newmt
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: UB, LB
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: new, old
       INTEGER :: z, j, DIMEN(2) 
       
       DIMEN = shape(oldmt)
       ALLOCATE(new(DIMEN(1)))
       ALLOCATE(old(DIMEN(1)))
       
          do j=1,DIMEN(2)
            old = oldmt(:,j)
            do z=1,DIMEN(1)
                if ( old(z) .eq. LB(z) ) then 
                    new(z) = 0d0
                else if  ( (abs(LB(z)) + abs(UB(z))) .eq. 0 ) then
                    new(z) = 0d0
                else if  ( LB(z) .eq. UB(z) ) then
                    new(z) = 0d0
                else
                      new(z) = old(z) / ( UB(z) - LB(z) )
                end if
            end do
            newmt(:,j) = new
          end do
       DEALLOCATE(new)
       DEALLOCATE(old)
       
   END SUBROUTINE transform_norm





! ----------------------------------------------------------------------
! SUBROUTINE matrix_eucl_dist
! ----------------------------------------------------------------------
   SUBROUTINE matrix_eucl_dist ( a, b , res)
       IMPLICIT NONE
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: a, b
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: res
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: sum
       INTEGER :: i, j, z, DIMEN(2), DIMEN2(2)
       
       DIMEN = shape(res)
       DIMEN2 = shape(a)
       do i=1,  DIMEN(1)
           do j=1, DIMEN(2)
               sum = 0d0
               do z=1,DIMEN2(1)
                    sum =  sum + (a(z,j) - b(z,i))**2d0
               end do
               res(i,j) =  REAL(sqrt(dabs(sum)), KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))             
           end do
       end do
   END SUBROUTINE matrix_eucl_dist   
   
! ----------------------------------------------------------------------
! SUBROUTINE ssm_isdif2
! ----------------------------------------------------------------------
   SUBROUTINE ssm_isdif2(x,group,tol,flag,f,ind,ind2)
       IMPLICIT NONE
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: x
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: group
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: tol 
       INTEGER, INTENT(IN) :: flag
       INTEGER :: i, DIMEN(2)
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: f
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: eps
       INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ind, ind2
       INTEGER, DIMENSION(:), ALLOCATABLE :: ind_aux, aaa, indexfind
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), ALLOCATABLE, DIMENSION(:) :: num, denom, difference
       CHARACTER(len = 2) :: typesearch
       
       typesearch = "LT"
       eps = 2.2204d-016
       DIMEN = shape(group)
       f=0d0
       
       do i=1,DIMEN(2)
           ALLOCATE(num(size(x)))
           ALLOCATE(denom(size(x)))
           ALLOCATE(difference(size(x)))
           num = abs(x-group(:,i) )
           denom = min(abs(x),abs(group(:,i)))
           CALL find_in_vector(denom, eps, indexfind, typesearch)
           if (ALLOCATED(indexfind)) then
               denom(indexfind) = 1
               DEALLOCATE(indexfind)
           end if
           difference = num / denom
           typesearch = "GE"
           CALL find_in_vector(difference, tol, aaa, typesearch)


           if (flag .eq. 1 ) then
               if (.not. ALLOCATED(aaa)) then
                   f = f + 1d0
                   if ( .not. ALLOCATED( ind ) ) then
                       ALLOCATE(ind (1) )
                       ind(1) = i
                   else
                       ALLOCATE( ind_aux ( size(ind) + 1 ) )
                       ind_aux(size(ind) + 1) = i
                       CALL MOVE_ALLOC(ind_aux,ind)
                   end if
               end if
           else if ( flag .eq. 2) then
               if  ( ALLOCATED(aaa) .and. (size(aaa) .ne. size(x) ) ) then 
                       if ( .not. ALLOCATED( ind ) ) then
                            ALLOCATE(ind (1) )
                            ind(1) = i
                       else
                            ALLOCATE( ind_aux ( size(ind) + 1 ) )
                            ind_aux(size(ind) + 1) = i
                            CALL MOVE_ALLOC(ind_aux,ind)
                       end if
               else
                       if ( .not. ALLOCATED( ind2 ) ) then
                            ALLOCATE(ind2 (1) )
                            ind2(1) = i
                       else
                            ALLOCATE( ind_aux ( size(ind2) + 1 ) )
                            ind_aux(size(ind2) + 1) = i
                            CALL MOVE_ALLOC(ind_aux,ind2)
                       end if
               end if
           end if
           
           if (ALLOCATED( denom) )DEALLOCATE(denom)
           if (ALLOCATED( difference) )DEALLOCATE(difference)
           if (ALLOCATED( num) )DEALLOCATE(num)
           if (ALLOCATED( aaa) ) DEALLOCATE(aaa)
           
       end do
       
   END SUBROUTINE ssm_isdif2   

! ----------------------------------------------------------------------
! SUBROUTINE ssm_local_filters
! ----------------------------------------------------------------------       
   SUBROUTINE ssm_local_filters(exp1, problem1, fitnessfunction, opts1,local_solver_var, common_vars, results, time, &
                    childset,xbest,fbest,numeval,fbest_lastiter,refset,refset_change,centercriteria) 
        IMPLICIT NONE
        TYPE(algorithm_common_vars) :: common_vars
        TYPE(time_ess) :: time        
        TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_var
        TYPE(C_PTR), INTENT(INOUT) :: exp1 
        TYPE(resultsf), INTENT(INOUT) :: results
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
        TYPE(opts), INTENT(INOUT) :: opts1
        TYPE(Refsettype), INTENT(INOUT) :: childset, refset
        INTEGER :: clock_rate, clock_start, clock_stop 
        INTEGER, INTENT(IN) ::  centercriteria
        INTEGER(C_LONG), INTENT(INOUT) :: numeval
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),  DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: fbest_lastiter
        TYPE(problem), INTENT(INOUT) :: problem1
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: refset_change 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:,:),ALLOCATABLE :: new_initial_points, &
                                                                                        new_local_solutions
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:),ALLOCATABLE :: new_local_solutions_values, &
                                                                                                   new_initial_points_values
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: xbest,fbest
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: temp_fpen
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: minimo, X0, X, dismin, child_points, &
                                                                             dismin_save
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: local_init, child_norm, &
                                                                                local_init_norm, repmat, dist
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: F0, indexj, fval, f11
        INTEGER, DIMENSION(:), ALLOCATABLE :: I, indrefset, indmin, ind11, ind12
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ::  calctimeMPI, cputime1
        INTEGER :: ind,j,nconst, DIMEN(2), DIMEN2(2), adiccionar_local
        TYPE(outfuction) :: outfunct
        INTEGER :: MODE,  sizelocal
        INTEGER :: idp, par, getidp
        INTEGER(C_LONG) :: last_evals
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: MAX_DBL

        CALL setdblmax(MAX_DBL)
        par = 0
        idp = 0
       
        idp = getidp(exp1)
       
        ALLOCATE(X0(common_vars%nvar))
        ALLOCATE(X(common_vars%nvar))   
    if (opts1%localoptions%empty .eq. 0) then
          if (local_solver_var%stage_1 .eq. 1) then           
            ALLOCATE(minimo(1) )
            minimo = minval(childset%fpen)
            ALLOCATE(I(1))
            I = minloc(childset%fpen)
            ind = I(1)
            if (minimo(1) .LT. fbest(1) ) then
                X0 =  childset%x(:,ind)
                F0 =  minimo(1)
            else
                X0 = xbest
                F0 = fbest(1)
            end if
            DEALLOCATE(minimo)
            local_solver_var%n_critico = opts1%localoptions%n2
            
            if (opts1%localoptions%bestx .eq. 1) then
                local_solver_var%use_bestx=0
            end if
            
            
          else if (local_solver_var%stage_2 .eq. 1 ) then 
            if ((opts1%localoptions%bestx .eq. 1).AND.(local_solver_var%use_bestx .eq. 1)) then
                X0 = xbest
                F0 = fbest(1)
                local_solver_var%use_bestx = 0
            else
                    DIMEN  =  shape(local_solver_var%local_solutions)
                    DIMEN2 =  shape(local_solver_var%initial_points)
                    ALLOCATE( local_init(common_vars%nvar,DIMEN2(2) + DIMEN(2) ))
                    local_init(:, 1:DIMEN(2) ) = local_solver_var%local_solutions
                    local_init(:, (DIMEN(2)+1): (DIMEN2(2) + DIMEN(2)) ) = local_solver_var%initial_points
                    
                    DIMEN =  shape(local_init)
                    ALLOCATE( local_init_norm(common_vars%nvar, DIMEN(2)) )
                    DIMEN =  shape(childset%x)
                    ALLOCATE( child_norm(common_vars%nvar, DIMEN(2)) )
                    
                    CALL transform_norm(local_init, local_init_norm, problem1%XU, problem1%XL)
                    CALL transform_norm(childset%x, child_norm, problem1%XU, problem1%XL)
                    
                    DIMEN = shape(local_init_norm)
                    DIMEN2 = shape(child_norm)
                    
                    ALLOCATE(dist(DIMEN(2), DIMEN2(2)))
                    CALL matrix_eucl_dist(child_norm,local_init_norm, dist )
                           
                    DIMEN = shape(dist)
                    ALLOCATE(dismin( DIMEN(2) ))
                    ALLOCATE(dismin_save( DIMEN(2)) )
 
                    IF ( opts1%localoptions%balance .GE. 0 ) THEN
                        dismin = minval(dist,DIM=1)
                        dismin_save = dismin
                        dismin = (-1.0d0) * dismin
                        ALLOCATE(indmin(  size(dismin) ))
                        call QsortC_ind(dismin, indmin)
                    
                        ALLOCATE(temp_fpen( size(childset%fpen) ))
                        temp_fpen=childset%fpen
                        ALLOCATE(indrefset( size(childset%fpen) ))
                        call QsortC_ind(temp_fpen, indrefset)   
                  
 
                        ALLOCATE(child_points(size(childset%fpen)))
                        child_points=0d0                        
                        do j=1,size(childset%fpen)
                            indexj = REAL(j,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
                            child_points(indmin(j)) = &
                                REAL(child_points(indmin(j)),KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))  &
                                + opts1%localoptions%balance* indexj
                            child_points(indrefset(j)) = &
                                   REAL(child_points(indrefset(j)),KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) &
                                 + ( 1d0 - opts1%localoptions%balance)* indexj                       
                        end do
                        
                        DEALLOCATE(indrefset)             
                    
                        ALLOCATE(minimo(1) )
                        minimo = minval(child_points)
                        ALLOCATE(I(1))
                        I = minloc(child_points)
                    
                        X0 = childset%x(:,I(1) )
                        F0 = childset%fpen( I(1) )
                                                
                    ELSE 
                      !  PRINT *, "*****NO BALANCE"
                        ALLOCATE(temp_fpen( size(childset%fpen) ))
                        temp_fpen=childset%fpen
                        ALLOCATE(indrefset( size(childset%fpen) ))
                        call QsortC_ind(temp_fpen, indrefset)   
                        dismin = minval(dist,DIM=1)
                        ALLOCATE(I(1))
                        
                        do j=1,size(childset%fpen)
                            IF ( dismin(indrefset(j)) .GT. opts1%localoptions%threshold_local ) THEN
                                I(1) = indrefset(j)
                                X0 = childset%x(:, I(1) )
                                F0 = childset%fpen( I(1) ) 
                                EXIT
                            END IF      
                        end do  
                        
                        DEALLOCATE(indrefset)  

                    END IF
                    
                    DEALLOCATE(dismin_save)
                    DEALLOCATE(dismin)                    
                    if ( F0 .EQ. MAX_DBL ) then
                        X0 = xbest
                        F0 = fbest(1)
                    end if
                !else
                !    X0 = xbest
                !    F0 = fbest(1)
                !end if
            end if
        end if 
        
        if ( opts1%useroptions%iterprint .eq. 1 ) then
            if ( opts1%localoptions%solver == 'n2sol' ) then
                write (*,*)  idp, "Call local solver: Nl2sol - INPUT: Initial point function value: [", f0, "]"
            else if ( opts1%localoptions%solver == 'misqp' ) then
                write (*,*)  idp, "Call local solver: misqp - INPUT: Initial point function value: [", f0, "]"
            end if
#ifdef MPI2            
            CALL printlsinitlog( exp1, f0 )
#endif            
        end if
        
        X=X0
        
        CALL setentervalueess( exp1, fbest);
        CALL printverboselocaloutput(exp1, size(x0), X, f0, idp)
         
#ifdef TIMELOCAL     
#ifdef MPI2 
        cputime1 = calctimeMPI(exp1,time%starttime)
#else
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        CALL SYSTEM_CLOCK(count=clock_start)
#endif
#endif 
        
#ifndef MPI2
        WRITE (*, '(A33,F25.10)') 'Call to local solver. Input Fx:', f0        
#endif        
        last_evals = numeval 
        CALL call_local_solver(problem1,exp1,opts1,fitnessfunction, x, fval,common_vars%NPROC,numeval,centercriteria, &
                        local_solver_var)

        nconst = problem1%ineq
        outfunct = ssm_evalfc(exp1,fitnessfunction,x,problem1, opts1, nconst, 0)
        numeval = numeval + 1
        fval = outfunct%value
         
#ifdef TIMELOCAL    
#ifdef MPI2
        time%localsolvertime = calctimeMPI(exp1,time%starttime)
        time%localsolvertime = time%localsolvertime - cputime1
#else        
        CALL SYSTEM_CLOCK(count=clock_stop)
        time%localsolvertime = time%localsolvertime + REAL(clock_stop-clock_start)/REAL(clock_rate)
#endif
#endif    
        CALL printverboselocaloutput2(exp1, fval,idp)
        CALL improvelocalsolver(exp1,fval,f0 )
        
        if ( opts1%useroptions%iterprint .eq. 1 ) then
            last_evals = numeval - last_evals
#ifdef MPI2               
            CALL printlsendlog( exp1, fval, time%localsolvertime, last_evals )
#else 
            WRITE (*, '(A31,F25.10)') 'Stop local solver. Output Fx:', fval  
#endif            
        end if
        
        

        if (.not. ALLOCATED(local_solver_var%initial_points) ) then 
            ALLOCATE(local_solver_var%initial_points(common_vars%nvar,1) )
            local_solver_var%initial_points(:,1) = X0
            ALLOCATE(local_solver_var%initial_points_values(1))
            local_solver_var%initial_points_values(1) = f0
        else 
            if (size(local_solver_var%initial_points_values) .LT.  (opts1%globaloptions%dim_refset)*100) then
                DIMEN = shape(local_solver_var%initial_points)
                ALLOCATE(new_initial_points(DIMEN(1), DIMEN(2) +1 ) )
                new_initial_points(:, 1:DIMEN(2)) = local_solver_var%initial_points
                new_initial_points(:, DIMEN(2) +1) = X0
                CALL MOVE_ALLOC(new_initial_points, local_solver_var%initial_points)
                
                ALLOCATE(new_initial_points_values(size(local_solver_var%initial_points_values)+1 ) )
                new_initial_points_values(1:size(local_solver_var%initial_points_values)) = local_solver_var%initial_points_values
                new_initial_points_values(size(local_solver_var%initial_points_values)+1) = F0
                CALL MOVE_ALLOC(new_initial_points_values, local_solver_var%initial_points_values)
                
            else
                DIMEN = shape(local_solver_var%initial_points)
                
                ALLOCATE( new_initial_points(DIMEN(1),DIMEN(2)))
                new_initial_points(:,1:DIMEN(2)-1) = new_initial_points(:,2:DIMEN(2))
                new_initial_points(:,DIMEN(2)) = X0 
                CALL MOVE_ALLOC(new_initial_points, local_solver_var%initial_points)
                 

                ALLOCATE(new_initial_points_values(size(local_solver_var%initial_points_values)))
                new_initial_points_values(1:size(local_solver_var%initial_points_values)-1) = &
                                                                local_solver_var%initial_points_values(2:DIMEN(2))
                new_initial_points_values( size(local_solver_var%initial_points_values) ) = F0
                CALL MOVE_ALLOC(new_initial_points_values, local_solver_var%initial_points_values)

            end if
            sizelocal = size ( local_solver_var%initial_points_values )
        end if
        
        
        if (local_solver_var%stage_1 .eq. 1) then
            local_solver_var%stage_1 = 0
            local_solver_var%stage_2 = 1
        end if
        
        if ( (outfunct%include .eq. 1) .and. ( outfunct%pena .le. opts1%useroptions%tolc )) then
            if ( outfunct%value_penalty .lt. fbest(1) ) then 
                fbest(1) = outfunct%value_penalty 
                xbest = X
                if ( fbest(1) .lt. fbest_lastiter(1) ) then
                        if (.not. ALLOCATED(results%fbest) ) ALLOCATE(results%fbest(size(fbest)))
                        results%fbest(1)=fbest(1);
                        if (.not. ALLOCATED(results%xbest) ) ALLOCATE(results%xbest(size(xbest)))
                        results%xbest=xbest
                        results%numeval=numeval
                end if
            end if 
        end if
        
        
        if (outfunct%include .eq. 1) then
            if (ALLOCATED(local_solver_var%local_solutions)) then
                adiccionar_local=1
                CALL ssm_isdif2(outfunct%x,local_solver_var%local_solutions,1d-2,1,f11,ind11,ind12)
                sizelocal = size(local_solver_var%local_solutions_values)
                if (f11 .ge. 1d0 ) then 
                    adiccionar_local = 0
                end if
            else 
                adiccionar_local=1
            end if
          
            if ( adiccionar_local .eq. 1 ) then 
                 
                if ( outfunct%value_penalty .LT. refset%fpen(childset%parent_index(I(1))) ) then
#ifdef MPI2                       
                    CALL printlocalsolverinsert(exp1,outfunct%value_penalty,local_solver_var%local_solutions_values,&
                                 sizelocal, 1, outfunct%value_penalty, refset%fpen(childset%parent_index(I(1))))
#endif                                 
                    refset%x(:,childset%parent_index(I(1) )) = outfunct%x
                    refset%f(childset%parent_index(I(1) )) = outfunct%value
                    refset%fpen(childset%parent_index(I(1) )) = outfunct%value_penalty
                    refset%penalty(childset%parent_index(I(1) )) = outfunct%pena
                    if (ALLOCATED (refset%nlc) ) refset%nlc(:,childset%parent_index(I(1) )) = outfunct%nlc
                    refset_change (childset%parent_index(I(1) )) = 0
                else
#ifdef MPI2                       
                    CALL printlocalsolverinsert(exp1, outfunct%value_penalty, local_solver_var%local_solutions_values, &
                                sizelocal, 0, outfunct%value_penalty, refset%fpen(childset%parent_index(I(1))))
#endif                                
                end if
                
                if (.not. ALLOCATED(local_solver_var%local_solutions)) then
                    ALLOCATE(local_solver_var%local_solutions(common_vars%nvar,1))
                    local_solver_var%local_solutions(:,1) = outfunct%x
                    ALLOCATE(local_solver_var%local_solutions_values(1) )
                    local_solver_var%local_solutions_values(1) = outfunct%value_penalty
                else
                    if (size(local_solver_var%local_solutions_values) .LT. &
                                                                (opts1%globaloptions%dim_refset*100) ) then                        
                        DIMEN = shape( local_solver_var%local_solutions )
                        ALLOCATE( new_local_solutions(DIMEN(1),DIMEN(2)+1))
                        new_local_solutions(:,1:DIMEN(2)) = local_solver_var%local_solutions
                        new_local_solutions(:,DIMEN(2)+1) = outfunct%x
                        CALL MOVE_ALLOC(new_local_solutions, local_solver_var%local_solutions)
                    
                        ALLOCATE(new_local_solutions_values(size(local_solver_var%local_solutions_values)+1))
                        new_local_solutions_values(1:size(local_solver_var%local_solutions_values)) = &
                                                                                local_solver_var%local_solutions_values
                        new_local_solutions_values(size(local_solver_var%local_solutions_values)+1) = outfunct%value_penalty
                        CALL MOVE_ALLOC(new_local_solutions_values, local_solver_var%local_solutions_values)
 
                    else
                        DIMEN = shape( local_solver_var%local_solutions )
                        ALLOCATE( new_local_solutions(DIMEN(1),DIMEN(2)))
                        new_local_solutions(:,1:(DIMEN(2)-1)) = local_solver_var%local_solutions(:,2:DIMEN(2))
                        new_local_solutions(:,DIMEN(2)) = outfunct%x
                        CALL MOVE_ALLOC(new_local_solutions, local_solver_var%local_solutions)

                        ALLOCATE(new_local_solutions_values(size(local_solver_var%local_solutions_values)))
                        new_local_solutions_values(1:(size(local_solver_var%local_solutions_values)-1))= &
                                                                local_solver_var%local_solutions_values(2:DIMEN(2))
                        new_local_solutions_values(size(local_solver_var%local_solutions_values))= outfunct%value_penalty
                        CALL MOVE_ALLOC(new_local_solutions_values,local_solver_var%local_solutions_values)
                          
                    end if
                end if
            else
#ifdef MPI2                   
                CALL printlocalsolverinsert(exp1,outfunct%value_penalty, local_solver_var%local_solutions_values, &
                                sizelocal, 0, outfunct%value_penalty, refset%fpen(childset%parent_index(I(1))))
#endif                                
            end if
        end if
       
        
        if (ALLOCATED(minimo) ) DEALLOCATE(minimo)
        if (ALLOCATED(repmat) ) DEALLOCATE(repmat)
        if (ALLOCATED(local_init) ) DEALLOCATE(local_init)
        if (ALLOCATED(child_norm) ) DEALLOCATE(child_norm)
        if (ALLOCATED(local_init_norm) ) DEALLOCATE( local_init_norm)
        if (ALLOCATED(I) ) DEALLOCATE(I)
        if (ALLOCATED(X0) ) DEALLOCATE(X0)
        if (ALLOCATED(X) ) DEALLOCATE(X)
        if (ALLOCATED(new_initial_points) ) DEALLOCATE(new_initial_points)
        if (ALLOCATED(new_local_solutions) ) DEALLOCATE(new_local_solutions)
        if (ALLOCATED(dismin) ) DEALLOCATE(dismin)
        if (ALLOCATED(indmin) ) DEALLOCATE(indmin)
        if (ALLOCATED(dist) ) DEALLOCATE(dist)
        if (ALLOCATED(child_points) ) DEALLOCATE(child_points)
        if (ALLOCATED(ind11) ) DEALLOCATE(ind11)
        if (ALLOCATED(ind12) ) DEALLOCATE(ind12)
        if (ALLOCATED(temp_fpen) ) DEALLOCATE(temp_fpen)
        CALL destroy_outfuction(outfunct)
    
        if (opts1%localoptions%n2 .LE. 0)  then 
            opts1%localoptions%empty=1
        end if
    
    end if
    
    END SUBROUTINE ssm_local_filters
    



      
    
! ----------------------------------------------------------------------
! SUBROUTINE local_solver
! ----------------------------------------------------------------------    
    SUBROUTINE destroy_local_solver_vars(local_solver_var)
        TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_var
        
        IF (ALLOCATED(local_solver_var%initial_points)) DEALLOCATE(local_solver_var%initial_points)
        IF (ALLOCATED(local_solver_var%initial_points_values)) DEALLOCATE(local_solver_var%initial_points_values)
        IF (ALLOCATED(local_solver_var%local_solutions)) DEALLOCATE(local_solver_var%local_solutions)
        IF (ALLOCATED(local_solver_var%local_solutions_values)) DEALLOCATE(local_solver_var%local_solutions_values)        
    END SUBROUTINE destroy_local_solver_vars
    
    
    
END MODULE localsolver
