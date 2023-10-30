#define TIMELOCAL 1

MODULE dhc_mod
    USE iso_c_binding
    USE scattersearchtypes
    USE common_functions
    USE funcevalinterface
CONTAINS


SUBROUTINE insert_in_bounds(xv)
    IMPLICIT NONE
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: xv
    CHARACTER(len = 2) :: typesearch
    INTEGER, DIMENSION(:), ALLOCATABLE :: indexfind
     
    typesearch = "LT"
    CALL find_in_vector(xv, 0d0, indexfind, typesearch)
    if (ALLOCATED(indexfind)) then
        CALL reajust_index(indexfind)
        xv(indexfind) = 0
        DEALLOCATE(indexfind)
    end if

    typesearch = "GT"
    CALL find_in_vector(xv, 1d0, indexfind, typesearch)
    if (ALLOCATED(indexfind)) then
        CALL reajust_index(indexfind)
        xv(indexfind) = 1
        DEALLOCATE(indexfind)
    end if

END SUBROUTINE insert_in_bounds



SUBROUTINE check_bounds(xvreal, XL, XU, pen2)
    IMPLICIT NONE
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: xvreal, XL, XU  
    INTEGER :: auxsize
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: pen2
    INTEGER, DIMENSION(:), ALLOCATABLE :: index_result2
    CHARACTER(len = 2) :: typesearch    
    INTEGER, DIMENSION(:), ALLOCATABLE :: aa
        

    ALLOCATE(index_result2( size(xvreal) ))
    typesearch = "LT"
    index_result2 = compare_multiple_vector(xvreal, XL, auxsize, typesearch)
    ALLOCATE(aa(auxsize))
    CALL index_of_ones(index_result2, aa)
    CALL reajust_index(aa)
    DEALLOCATE(index_result2)
    if (ALLOCATED(aa)) then
        pen2 = pen2 + SUM(XL(aa) - xvreal(aa))
        DEALLOCATE(aa)
    end if

    ALLOCATE(index_result2(size(xvreal)))
    typesearch = "GT"
    index_result2 = compare_multiple_vector(xvreal, XU, auxsize, typesearch)
    ALLOCATE(aa(auxsize))
    CALL index_of_ones(index_result2, aa)
    CALL reajust_index(aa)
    DEALLOCATE(index_result2)
    if (ALLOCATED(aa)) then
        pen2 = pen2 + SUM(xvreal(aa) - XU(aa))
        DEALLOCATE(aa)
    end if
END SUBROUTINE check_bounds              

SUBROUTINE  dhc(problem1,exp1,opts1,fitnessfunction,X0,initsize,thres,budget,numeval,fx)
        IMPLICIT NONE
        TYPE(problem), INTENT(INOUT) :: problem1
        TYPE(opts), INTENT(INOUT) :: opts1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction 
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: initsize, fx
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: X0        
        INTEGER, INTENT(INOUT) :: budget
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: thres
        ! OTHER PARAMETERS
        INTEGER :: n_out, NDIM
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: pen2
        INTEGER(C_LONG), INTENT(INOUT) :: numeval
        INTEGER(C_LONG) :: maxiter, iter
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE:: v,u,xreal,xv, xvreal, cxv, cx   
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: vvec, vr, fxv, penalty
        REAL(C_DOUBLE) :: callfitnessfunctionfortranopenmp2
        INTEGER ::  exitwhile, exitwhile2, vi
        INTEGER(C_INT) :: mig, dest, cooperativempitestess
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: NAN
        INTEGER(C_INT) :: error, setNAN
        INTEGER(C_INT) :: checkcooperativemigrationcriteriacessinner
        INTEGER :: ii, NPROC
        INTEGER :: dist_criteria

        CALL getdistcriteria(exp1,dist_criteria)
 
        numeval = 0
        error = setNAN(NAN)
        
        if ( (ALLOCATED(problem1%CL) ) .OR. (ALLOCATED(problem1%CU)) ) then
            n_out = 2
            ALLOCATE(cxv(problem1%ineq))
            ALLOCATE(cx(problem1%ineq))             
        else           
            n_out = 1
        end if
        
        NDIM = size(x0)

        ALLOCATE(v(NDIM))
        ALLOCATE(u(NDIM))
        ALLOCATE(xv(NDIM))
        ALLOCATE(xreal (NDIM) )
        ALLOCATE(xvreal(NDIM))
        v = 0d0
        u = v
        vi = -1
        vvec = 1
        vr = - initsize
         
        xreal = X0 * (problem1%XU - problem1%XL) + problem1%XL

        
        if (n_out .GT. 1) then
            fx = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xreal, cx, 0)
            penalty = ssm_penalty_function(cx,problem1%CL,problem1%CU, opts1%useroptions%tolc)
            fx = fx + opts1%useroptions%weight * penalty
        else
            fx = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xreal, cx, 0)
        end if
        
        numeval = numeval + 1
        fxv = 1d30
        exitwhile2 = 0
        mig=0
       ! !print *, "ABS(vr)", ABS(vr), ".GE.", "thres", thres,"exitwhile2",exitwhile2
        do while (( ABS(vr) .GE. thres  ) .AND. (exitwhile2 .NE. 1) .AND. (mig .NE. 1) ) 
            
          if (mig .NE. 1) THEN  
            if (ABS(vr) .LT. 2d0*thres ) then 
                maxiter = 2 * NDIM
            else
                maxiter = 2
            end if
            iter = 0
            exitwhile = 0
                    
            do while ( (fxv .GE. fx ) .AND. ( iter .LT. maxiter ) .AND.(exitwhile .NE. 1) .AND. (mig .NE. 1))
#ifdef MPI2        
             CALL setNPROC(exp1, NPROC)
             if (dist_criteria .EQ. 1) then
                do ii=1,NPROC
                        dest = ii-1
                        mig = cooperativempitestess(exp1, dest )
                        if (mig .eq. 1) then
                                EXIT
                        end if
                end do 
!                if (mig .NE. 1) mig =  checkcooperativemigrationcriteriacessinner(exp1);
             end if
#endif 
            if (mig .NE. 1) then
                if (iter .EQ. 0) then
                    xv = x0
                else
                    xv(vi+1)=xv(vi+1)-vr
                end if
                if (vvec .EQ. 1) vvec = 0
                vr = -vr
                if (vr .GT. 0) vi = MOD(vi+1, NDIM )
                
                xv(vi+1) = xv(vi+1)+vr
                
                CALL insert_in_bounds(xv)   

                xvreal = xv * (problem1%XU - problem1%XL) + problem1%XL

                if (n_out .GT. 1) then
                    fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal, cxv, 0)
                    penalty = ssm_penalty_function(cxv,problem1%CL,problem1%CU, opts1%useroptions%tolc)
                    fxv = fxv + opts1%useroptions%weight * penalty
                else
                    fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal, cxv, 0)
                end if 

                !print *, " 2 - fxv", fxv,"fx", fx, "vr", vr
                pen2 = 0d0
                
                CALL check_bounds(xvreal, problem1%XL, problem1%XU, pen2)               
                
                fxv = fxv + opts1%useroptions%weight * pen2
                
                numeval = numeval + 1
                iter=iter+1
                
                if (numeval .GE. budget) then
                    vr=thres/10d0
                    exitwhile = 1
                end if
            end if    
            end do
            
            if ( (fxv .GE. fx) .OR. (fxv .EQ. NAN) ) then
                fxv = 1d30
                vr=vr/2d0
            else
                fx=fxv
                x0=xv
              
                if (opts1%localoptions%iterprint .EQ. 1) print *, "****Nevals:",numeval,"Fobj:",fx
                if (iter==0) then
                    if (vvec .EQ. 1) then
                        u = u+v
                        v=v*2d0
                        xv=xv+v
                        vr=vr*2d0
                    else
                        u(vi+1)=u(vi+1)+vr     
                        vr=vr*2d0
                        xv(vi+1)=xv(vi+1)+vr                      
                    end if
                    
                    CALL insert_in_bounds(xv)
                       
                    xvreal = xv * (problem1%XU - problem1%XL) + problem1%XL
                    
                    if (n_out .GT. 1) then
                        fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal,  cxv, 0)
                        penalty = ssm_penalty_function(cxv,problem1%CL,problem1%CU, opts1%useroptions%tolc)
                        fxv = fxv + opts1%useroptions%weight * penalty
                    else
                        fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal,  cxv, 0)
                        !!print *, "      fxv 3", fxv                        
                    end if
                    !print *, " 3 - fxv", fxv,"fx", fx, "vr", vr
                    pen2=0d0
                    CALL check_bounds(xvreal, problem1%XL, problem1%XU, pen2)    
                        
                    fxv = fxv + opts1%useroptions%weight * pen2 
                    numeval=numeval+1
                    
                    if (numeval .GE. budget) then
                        vr=thres/10d0
                        exitwhile2 = 1
                    end if
                else    
                    xv=xv+u
                    xv(vi+1)=xv(vi+1)+vr
                    CALL insert_in_bounds(xv)
                    
                    xvreal = xv * (problem1%XU - problem1%XL) + problem1%XL
                    
                    if (n_out .GT. 1) then
                        fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal,  cxv, 0)
                        penalty = ssm_penalty_function(cxv,problem1%CL,problem1%CU, opts1%useroptions%tolc)
                        fxv = fxv + opts1%useroptions%weight * penalty
                    else
                        fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal, cxv, 0)
                    end if                    
                    
                    pen2=0d0
                    CALL check_bounds(xvreal, problem1%XL, problem1%XU, pen2)                      
                    !print *, " 4 - fxv", fxv,"fx", fx, "vr", vr
                    fxv = fxv + opts1%useroptions%weight * pen2 
                    numeval=numeval+1
                    
                    if (numeval .GE. budget) then
                        vr=thres/10d0
                        exitwhile2 = 1
                    end if
                    !!print *, "fxv4.2", fxv,"fx", fx, "vr", vr
                    if ( ((fxv .GE. fx).AND.(exitwhile2 .NE. 1)) .OR. ((fxv .EQ. NAN).AND.(exitwhile2 .NE. 1))) then
                        u=0d0
                        xv=x0
                        u(vi+1)=vr              
                        vr=vr*2d0
                        xv(vi+1)=xv(vi+1)+vr   
                        
                        CALL insert_in_bounds(xv)   
                        xvreal = xv * (problem1%XU - problem1%XL) + problem1%XL
                        
                        if (n_out .GT. 1) then
                            fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal, cxv, 0)
                            penalty = ssm_penalty_function(cxv,problem1%CL,problem1%CU, opts1%useroptions%tolc)
                            fxv = fxv + opts1%useroptions%weight * penalty
                        else
                            fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal,  cxv, 0)
                        end if  
                        !print *, " 5 - fxv", fxv,"fx", fx, "vr", vr
                        pen2=0d0
                        CALL check_bounds(xvreal, problem1%XL, problem1%XU, pen2)                      
                        fxv = fxv + opts1%useroptions%weight * pen2 
                        numeval=numeval+1    
                        
                        if (numeval .GE. budget) then
                            vr=thres/10d0
                            exitwhile2 = 1
                        end if                        
                            
                    else 
                        if (exitwhile2 .NE. 1) then
                            X0=xv
                            fx=fxv
                            u(vi+1)=u(vi+1)+vr
                            v=2d0*u
                            vvec=1
                            xv=xv+v                       
                        
                            CALL insert_in_bounds(xv)   
                            xvreal = xv * (problem1%XU - problem1%XL) + problem1%XL                        
                        
                            if (n_out .GT. 1) then
                                fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal,  cxv, 0)
                                penalty = ssm_penalty_function(cxv,problem1%CL,problem1%CU, opts1%useroptions%tolc)
                                fxv = fxv + opts1%useroptions%weight * penalty
                            else
                                fxv = callfitnessfunctionfortranopenmp2(fitnessfunction, exp1, xvreal, cxv, 0)
                            end if                          
                            pen2=0d0
                            CALL check_bounds(xvreal, problem1%XL, problem1%XU, pen2)  
                            fxv = fxv + opts1%useroptions%weight * pen2 
                            numeval=numeval+1
                        
                            if (numeval .GE. budget) then
                                vr=thres/10d0
                                exitwhile2 = 1
                            end if    
                        
                            vr=0d0
                            vr=SUM(v**2d0)
                            vr=SQRT(vr)
                        end if
                    end if
                end if
            end if
          end if  
        end do
        
        !print *, "END fx", fx
        
        if ((fxv .LT. fx) .AND. (fxv .NE. NAN) .AND.( mig .NE. 1 ) ) then
            fx = fxv
            x0 = xv
        end if
        
        x0=x0*(problem1%XU - problem1%XL)+problem1%XL
        
        if (ALLOCATED(cxv)) DEALLOCATE(cxv)
        if (ALLOCATED(cx)) DEALLOCATE(cx)
        DEALLOCATE(v)
        DEALLOCATE(u)
        DEALLOCATE(xv)
        DEALLOCATE(xreal)
        DEALLOCATE(xvreal)

    END SUBROUTINE dhc
    
    
    
END MODULE dhc_mod
