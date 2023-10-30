MODULE misqp_interface
    USE iso_c_binding
    USE scattersearchtypes
    USE common_functions
    USE qsort_module
    USE funcevalinterface

#ifdef OPENMP
    USE omp_lib
#endif

!    USE MISQP
CONTAINS


SUBROUTINE evaluate_gradient ( problem1,exp1,opts1, x, f,  g, ncont, nfunc, m, df, dg, fobj_mysqp, pass_misqp, &
    neq, fitnessfunction  )
    IMPLICIT NONE
    TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction
    TYPE(C_FUNPTR), INTENT(INOUT) :: fobj_mysqp
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(C_PTR), INTENT(INOUT) :: exp1
    TYPE(opts), INTENT(INOUT) :: opts1
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: x, g, df
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)  :: dg
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE  :: x_aux
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT)  :: f
    INTEGER, INTENT(INOUT) :: ncont, m, neq, pass_misqp
    INTEGER(C_LONG), INTENT(INOUT)  :: nfunc
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: eps, feps, epsa, epsi
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: callfitnessfunctionfortranopenmp2
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: wa1
    INTEGER :: index1, i, j

    if ( pass_misqp .EQ. 1 ) then
        ! USAR fobj_mysqp
    else
        ALLOCATE(x_aux(size(x)))
        eps = REAL(1.0d-7,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        index1=ncont
        feps=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        ALLOCATE(wa1(m))
        wa1=0
        do i=1,index1
            epsa = eps * MAX( REAL(1d-5,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ,ABS(x(i)))
            epsi =REAL(1d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/epsa
            x_aux(i) = x(i)
            x(i) = x(i) + epsa
            feps = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, x, wa1, 0)
            nfunc = nfunc  + 1
            df(i)=epsi*(feps-f)
            do j=1,m
                dg(j,i) = epsi *(wa1(j)-g(j))
            end do
            x(i) = x_aux(i)
        end do

        DEALLOCATE(wa1)
        DEALLOCATE(x_aux)
    end if


END SUBROUTINE evaluate_gradient

#ifdef OPENMP
SUBROUTINE evaluate_gradient_parallel ( problem1,exp1,opts1, x, f,  g, ncont, nfunc, m, df, dg, fobj_mysqp, pass_misqp, &
    neq, fitnessfunction  )
    IMPLICIT NONE
    TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction
    TYPE(C_FUNPTR), INTENT(INOUT) :: fobj_mysqp
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(C_PTR), INTENT(INOUT) :: exp1
    TYPE(opts), INTENT(INOUT) :: opts1
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: x, g, df
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)  :: dg
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE  :: x_aux
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT)  :: f
    INTEGER, INTENT(INOUT) :: ncont, m, neq, pass_misqp
    INTEGER(C_LONG), INTENT(INOUT)  :: nfunc
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: eps, feps, epsa, epsi
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: callfitnessfunctionfortranopenmp2
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: wa1
    INTEGER :: index1, i, j, IDPROC

    if ( pass_misqp .EQ. 1 ) then
        ! USAR fobj_mysqp
    else
        eps = REAL(1.0d-7,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        index1=ncont
        feps=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))

  !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)  FIRSTPRIVATE(IDPROC) PRIVATE(i,epsa,epsi,feps,j,x_aux,wa1)
  !  REDUCTION(+:nfunc) DEFAULT(SHARED)
        do i=1,index1
            ALLOCATE(wa1(m))
            ALLOCATE(x_aux(size(x)))
            IDPROC = OMP_GET_THREAD_NUM()
            wa1=0
            epsa = eps * MAX(REAL(1d-5,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ,ABS(x(i)))
            epsi=REAL(1d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))/epsa
            x_aux=x
            x_aux(i) = x(i) + epsa
            feps = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, x_aux, wa1, IDPROC)
            nfunc = nfunc  + 1
           
            df(i)=epsi*(feps-f)
            do j=1,m
                dg(j,i) = epsi *(wa1(j)-g(j))
            end do
            DEALLOCATE(wa1)
            DEALLOCATE(x_aux)
        end do
 !OMP END PARALLEL DO
    end if

END SUBROUTINE evaluate_gradient_parallel
#endif

SUBROUTINE fobj_misqp(fitnessfunction, exp1, x, g, neq, local_solver_vars, problem1)
   TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_vars
   TYPE(C_PTR), INTENT(INOUT) :: exp1
   TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction   
   INTEGER, INTENT(IN) :: neq
   REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
   REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: x
   TYPE(problem), INTENT(INOUT) :: problem1 

!global n_fun_eval ccll ccuu n_upper n_lower

!c=[];

!if ~isempty(ccll) | ~isempty(ccuu)
!    [fx,ggg] = feval(fobj,x,varargin{:});
!else
!    [fx] = feval(fobj,x,varargin{:});
!    ggg=[];
!end
!n_fun_eval = n_fun_eval + 1;


!ninequ=length(n_upper);
!for i=1:ninequ
!    c=[c ccuu(n_upper(i))-ggg(n_upper(i))];
!end

!nineql=length(n_lower);
!for j=1:nineql
!    c=[c ggg(n_lower(j))-ccll(n_lower(j))];
!end

!cx=[ggg(1:neq) c];   

END SUBROUTINE fobj_misqp

! callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, x, g_aux, 0)

SUBROUTINE ssm_aux_local_misqp (problem1,local_solver_vars)
    IMPLICIT NONE
    TYPE(problem), INTENT(INOUT) :: problem1
    TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_vars
    INTEGER, DIMENSION(:), ALLOCATABLE :: N_UPPER, N_LOWER, N_UPPER_AUX, N_LOWER_AUX
    INTEGER, DIMENSION(:), ALLOCATABLE :: CLINF, CUINF   
    INTEGER :: i, j, contador, delete
    CHARACTER(len = 2) :: typesearch
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: CCUU
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: CCLL    
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: C_INF, C_MINF

    CALL SETINF(C_INF)
    CALL SETMINF(C_MINF)
! SSM_AUX_LOCAL CODE
    IF (ALLOCATED(problem1%CU) ) THEN
        ALLOCATE(CCUU(size(problem1%CU)))
        ALLOCATE(CCLL(size(problem1%CL)))
        CCUU = problem1%CU
        CCLL = problem1%CL
           
        ALLOCATE(N_UPPER(size(problem1%CU)))
        ALLOCATE(N_LOWER(size(problem1%CL)))
        N_UPPER =  (/ (i, i = 1, size(problem1%CU)) /) 
        N_LOWER = N_UPPER
           
        typesearch = "EQ"
        CALL find_in_vector(problem1%CU, C_INF, CUINF, typesearch)  
        if (ALLOCATED(CUINF)) then
            CALL reajust_index(CUINF)
            ALLOCATE(N_UPPER_AUX( SIZE(N_UPPER) - SIZE(CUINF)  ) )
            contador = 1
            DO i=1, SIZE(N_UPPER)
                delete = 0
                DO j=1, SIZE(CUINF)
                    IF (i .EQ. CUINF(j) ) THEN
                        delete = 1
                        EXIT
                    END IF
                END DO
                IF (delete .NE. 1 ) THEN
                    N_UPPER_AUX(contador) = N_UPPER(i)
                    contador = contador + 1
                END IF
            END DO
            CALL MOVE_ALLOC(N_UPPER_AUX, N_UPPER)
            DEALLOCATE(CUINF)
        end if
        
        CALL find_in_vector(problem1%CL, C_MINF, CLINF, typesearch)
        if (ALLOCATED(CLINF)) then
            CALL reajust_index(CLINF)
           ALLOCATE(N_LOWER_AUX( SIZE(N_LOWER) - SIZE(CLINF)  ) )
            contador = 1
            DO i=1, SIZE(N_LOWER)
                delete = 0
                DO j=1, SIZE(CLINF)
                    IF (i .EQ. CLINF(j) ) THEN
                        delete = 1
                        EXIT
                    END IF
                END DO
                IF (delete .NE. 1 ) THEN
                    N_LOWER_AUX(contador) = N_LOWER(i)
                    contador = contador + 1
                END IF
            END DO
            CALL MOVE_ALLOC(N_LOWER_AUX, N_LOWER)     
            DEALLOCATE(CLINF)
        end if     
       
        if ( problem1%neq .GT. 0 ) then
                ALLOCATE(N_UPPER_AUX(size(N_UPPER)-problem1%neq))
                ALLOCATE(N_LOWER_AUX(size(N_LOWER)-problem1%neq)) 
                N_UPPER_AUX=N_UPPER(problem1%neq:size(N_UPPER))
                N_LOWER_AUX=N_LOWER(problem1%neq:size(N_LOWER))
                CALL MOVE_ALLOC(N_UPPER_AUX,local_solver_vars%N_UPPER)
                CALL MOVE_ALLOC(N_LOWER_AUX,local_solver_vars%N_LOWER)
        else 
                CALL MOVE_ALLOC(N_UPPER,local_solver_vars%N_UPPER) 
                CALL MOVE_ALLOC(N_LOWER,local_solver_vars%N_LOWER)
        end if
       
        
    END IF    
    
END SUBROUTINE ssm_aux_local_misqp


SUBROUTINE RUN_MISQP( problem1,exp1,opts1,fitnessfunction,acc, x0, fval, nfunc, ifail, g, local_solver_vars  )
       IMPLICIT NONE
        
       TYPE(problem), INTENT(INOUT) :: problem1
       TYPE(local_solver_help_vars), INTENT(INOUT) :: local_solver_vars    
       TYPE(opts), INTENT(INOUT) :: opts1
       TYPE(C_PTR), INTENT(INOUT) :: exp1
       TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
       TYPE(outfuction) :: outfunct

       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: acc, fval
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) ::g
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: x0
       INTEGER(C_LONG), INTENT(INOUT) :: nfunc
       INTEGER, INTENT(INOUT) ::ifail
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),DIMENSION(:), ALLOCATABLE :: x, g_aux, x_best, g_best

       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: U, DF
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: B, DG
       CHARACTER(len = 2) :: typesearch
       INTEGER :: i,j,delete, contador
       INTEGER :: n, m, me, mme, iprint, ngrad, ncont
              
       ! working arrays
       LOGICAL :: RESOPT
       INTEGER :: MAXIT, MAXPEN, MAXUND, MODE, NONMON, MAXNDE
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: accqp
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: f_ls, f, f_best
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: g_ls,x_ls,xl_ls,xu_ls,x_aux
       INTEGER :: lrw, liw, llw
       REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: RW
       LOGICAL, DIMENSION(:), ALLOCATABLE :: LW
       INTEGER, DIMENSION(:), ALLOCATABLE :: IW
       REAL(C_DOUBLE) :: callfitnessfunctionfortranopenmp2
       INTEGER :: exitloop, size_upper, size_lower, ite
       INTEGER :: idp
       INTEGER(C_INT) :: checkcooperativemigrationcriteriacessinner
       INTEGER :: mig, dest, cooperativempitestess
       INTEGER(C_INT) :: getopenmpoption, openmp_pos
       
       ncont = size(x0) - problem1%int_var - problem1%bin_var
       n = problem1%int_var + problem1%bin_var + ncont
     
       if (.not. ALLOCATED(local_solver_vars%N_UPPER)) then
            size_upper = 0
       else 
            size_upper = size(local_solver_vars%N_UPPER)
       end if
       
       if (.not. ALLOCATED(local_solver_vars%N_LOWER)) then 
            size_lower = 0
       else 
            size_lower = size(local_solver_vars%N_LOWER)
       end if
       
       m = size_upper+size_lower+problem1%neq
       me = problem1%neq
       ALLOCATE(x(size(x0)))
       x = x0
       
       
! RUN_MISQP
       f_ls = 0
       ALLOCATE(g_ls(m+me+1))
       g_ls(1:(m+me+1))=0
       ALLOCATE(g_aux(m))
       ALLOCATE(g(m)) 
       ALLOCATE(x_aux(n))
       ALLOCATE(x_best(n))
       ALLOCATE(g_best(m+me)) 
       ALLOCATE(x_ls(n+1))
       ALLOCATE(xl_ls(n+1))
       ALLOCATE(xu_ls(n+1))
       ALLOCATE(U (m+n+n) )
       ALLOCATE(B (n+1,n+1)  )
       ALLOCATE(DG (m+problem1%neq+1, n+1) )
       ALLOCATE(DF (n+1) )
       
       x_ls=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
       xl_ls=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
       xu_ls=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
       U=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
       B=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
       DG=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
       DF=REAL(0d0,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))

       ifail = 0
!   MAXIT :   Maximum number of iterations, where one iteration corresponds to
!             one evaluation of a set of gradients (e.g. 100).       
       maxit= 500
!   MAXPEN :  Maximum number of successive increments of the penalty parameter 
!             without success (e.g. 50).       
       maxpen=50
!   MAXUND :  Maximum number of successive iterations without improvements
!             of the iterate X (e.g. 10). 
       maxund=50
!   MAXNDE :  Maximum number of branch-and-bound steps for solving MIQP.       
       MAXNDE=10000
!   NONMON :  Maximum number of successive iterations, which are to be 
!             considered for the non-monotone trust region algorithm.       
       nonmon=2
       
       mode=0
       accqp = REAL(1d-12,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
       resopt=.TRUE.
       mme=m+me
       iprint = 0
       nfunc = 0
       ngrad = 0
       
       ! EVALUATION
       nfunc = nfunc + 1
       
       idp = 0
       f = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, x, g, idp)

       
       ! GRADIENT
       ngrad = ngrad + 1
#ifdef OPENMP       
       openmp_pos = getopenmpoption(exp1)
       if (openmp_pos .EQ. 1) then
       CALL evaluate_gradient_parallel( problem1, exp1, opts1, x, f,  g, ncont, nfunc, m, DF, DG, local_solver_vars%fobj_mysqp, &
                    local_solver_vars%pass_misqp, problem1%neq, fitnessfunction  )
       
       else 
       CALL evaluate_gradient( problem1, exp1, opts1, x, f,  g, ncont, nfunc, m, DF, DG, local_solver_vars%fobj_mysqp, &
                    local_solver_vars%pass_misqp, problem1%neq, fitnessfunction)
       end if
#else      
       CALL evaluate_gradient( problem1, exp1, opts1, x, f,  g, ncont, nfunc, m, DF, DG, local_solver_vars%fobj_mysqp, &
                    local_solver_vars%pass_misqp, problem1%neq, fitnessfunction)        
#endif
               
       x_best=x
       f_best=f
       g_best=g
      
       f_ls = f
       if (problem1%ineq .GT. 0) then
            g_ls(1:m) = g
       end if
     
       x_ls(1:n)= x
       xl_ls(1:n) = problem1%xl
       xu_ls(1:n) = problem1%xu
       exitloop = 0
       ite = 0
       
       DO WHILE ( (ifail .LE. 0) .AND. (exitloop .NE. 1))
           if (ifail .EQ. 0) then
                lrw = 22*n+11*m+6*me+2*nonmon+85
                liw=problem1%int_var+problem1%bin_var+16
                llw=15   
                IF (ALLOCATED( RW )) DEALLOCATE(RW)
                IF (ALLOCATED( IW )) DEALLOCATE(IW)
                IF (ALLOCATED( LW )) DEALLOCATE(LW)                      
                ALLOCATE(RW(lrw))
                ALLOCATE(IW(liw))
                ALLOCATE(LW(llw))
                RW = 0
                IW = 0
                LW = .FALSE.
                                 
           end if
           CALL MISQP(m,me,N,problem1%int_var,problem1%bin_var,x_ls,f_ls,g_ls,df,dg,U,xl_ls,xu_ls,B, &
                     acc,accqp,maxit,maxpen,maxund,resopt, &
                     nonmon,iprint,mode,ifail, &
                     RW,lrw,IW,liw,LW,llw,MAXNDE) 
           
           if (ifail .EQ. -1) then
               nfunc = nfunc+1
               x_aux=x_ls(1:n)
               g_aux=g_ls(1:m)
               f = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, x_aux, g_aux, 0)
               f_ls = f
               
               if (problem1%ineq .GT. 0) then
                        g_ls(1:m) = g(1:m)
               end if
               if (nfunc .GE. local_solver_vars%evals_per_iter) then 
                       ifail=0                       
                       x_ls(1:n)=x_best
                       f_ls=f_best
                       g_ls(1:m+me)=g_best 
               end if
           end if

           if (ifail .EQ. -2) then
               ngrad=ngrad+1
               g_aux=g_ls(1:m)
               x_aux=x_ls(1:n)
#ifdef OPENMP
               CALL evaluate_gradient_parallel( problem1, exp1, opts1, x_aux, f_ls,  g_aux,ncont, nfunc, m, DF, DG, &
                    local_solver_vars%fobj_mysqp, local_solver_vars%pass_misqp,problem1%neq, fitnessfunction  )
#else
               CALL evaluate_gradient( problem1, exp1, opts1, x_aux, f_ls,  g_aux,ncont, nfunc, m, DF, DG, &
                    local_solver_vars%fobj_mysqp, local_solver_vars%pass_misqp,problem1%neq, fitnessfunction  )

#endif
           end if
           
        
#ifdef MPI2
           dest = 0
           mig = cooperativempitestess(exp1, dest )
           if (mig .EQ. 1) exitloop = 1
#endif

           if (ifail .GE. 0) exitloop = 1
           ite = ite + 1

           if (f_ls .LT. f_best) then
                 x_best=x_ls(1:n)
                 f_best=f_ls
                 g_best=g_ls(1:m+me)                
           end if
       END DO
         
       x0=x_ls(1:n)
       fval=f_ls
       if (problem1%ineq .GT. 0) then
        g=g_ls(1:m)
       end if
       
       
       
       IF (ALLOCATED(U) ) DEALLOCATE(U)
       IF (ALLOCATED(B) ) DEALLOCATE(B)
       IF (ALLOCATED(DG) ) DEALLOCATE(DG)
       IF (ALLOCATED(DF) ) DEALLOCATE(DF)
       IF (ALLOCATED(g_ls) ) DEALLOCATE(g_ls)
       IF (ALLOCATED(x_ls) ) DEALLOCATE(x_ls)
       IF (ALLOCATED(g_best) ) DEALLOCATE(g_best)
       IF (ALLOCATED(x_best) ) DEALLOCATE(x_best)
       IF (ALLOCATED(xl_ls)) DEALLOCATE(xl_ls)
       IF (ALLOCATED(xu_ls)) DEALLOCATE(xu_ls)
       IF (ALLOCATED(x)) DEALLOCATE(x)
       IF (ALLOCATED(g_aux)) DEALLOCATE(g_aux)
       IF (ALLOCATED(x_aux)) DEALLOCATE(x_aux)
END SUBROUTINE RUN_MISQP
    
    
    
    
    
END MODULE misqp_interface
