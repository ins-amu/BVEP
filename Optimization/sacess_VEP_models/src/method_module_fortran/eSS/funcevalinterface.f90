!#define DETERMINISTIC 0

MODULE funcevalinterface
    USE iso_c_binding
    USE scattersearchtypes
    USE common_functions
#ifdef DETERMINISTIC
    USE HDF5   
#endif       
CONTAINS
! ----------------------------------------------------------------------
! SUBROUTINES SSM_EVALFC
! ----------------------------------------------------------------------
    FUNCTION ssm_evalfc(exp1, fitnessfunction,x,problem1, opts1, nconst, idopenmp) result(outt)
        IMPLICIT NONE
        TYPE(opts), INTENT(IN) :: opts1
        TYPE(problem), INTENT(IN) :: problem1
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction ! fobj
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: x
        INTEGER, INTENT(IN) :: nconst, idopenmp
        REAL(C_DOUBLE) :: value1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: value_penalty, pena
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: nlc
        TYPE(outfuction) :: outt
        REAL(C_DOUBLE) :: callfitnessfunctionfortranopenmp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: xcp
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: INF, NAN, DBL_MAX
        INTEGER(C_INT) :: error1, error2, error3, setinfinity, setnan, setdblmax
        error1 = setinfinity(INF)
        error2 = setnan(NAN)
        error3 = setdblmax(DBL_MAX)        
        
        !CALL setpoint(x)
        outt%include = 1
   
 
        if (( problem1%int_var .GT. 0) .OR. (problem1%bin_var .GT. 0) ) then
            ALLOCATE(xcp(size(x),1))
            xcp(:,1) = x
            call ssm_round_int(xcp, problem1%int_var + problem1%bin_var, problem1%XL, problem1%XU)
            x = xcp(:,1)
            DEALLOCATE(xcp)
        end if
        
        ! Adjust to bounds
        ! A correccion dos bounds está dentro da función de avaliación
        if (nconst .GT. 0) then
            ALLOCATE(nlc(problem1%ineq))
            value1 = callfitnessfunctionfortranopenmp(fitnessfunction, exp1, x, size(x), problem1%XL, problem1%XU, nlc, idopenmp)

            pena = ssm_penalty_function(nlc,problem1%CL,problem1%CU, opts1%useroptions%tolc)
            value_penalty=REAL(value1,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))+ &
                   REAL(opts1%useroptions%weight,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) *pena
        else      
            value1 = callfitnessfunctionfortranopenmp(fitnessfunction, exp1, x, size(x), problem1%XL, problem1%XU, nlc, idopenmp)
            pena = REAL(0.0d0, SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
            value_penalty=value1
        end if

!        if ((value1 .eq. INF ) .or. ( value1 .eq. NAN ) .or. (value1 .eq. DBL_MAX))  then
        if ((value1 .eq. INF ) .or. ( value1 .eq. NAN ) )  then
            outt%include = 0
            CALL sumfailevals(exp1)
        end if
        
!        print *, outt%value
        
        outt%value = REAL(value1,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        outt%value_penalty = value_penalty
        outt%pena = pena
        ALLOCATE(outt%x(size(x)))
        outt%x = REAL(x,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        if (nconst .GT. 0) then 
            ALLOCATE(outt%nlc(nconst))
            outt%nlc = REAL(nlc,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        end if
        if (ALLOCATED(nlc) ) DEALLOCATE(nlc)
   
    END FUNCTION ssm_evalfc
    
    
! ----------------------------------------------------------------------
! SUBROUTINES SSM_PENALTY_FUNCTION
! ----------------------------------------------------------------------
    FUNCTION ssm_penalty_function(nlc,c_L,c_U,tolc) result (res)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: nlc
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: c_U, c_L
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: tolc
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: res
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: P
        INTEGER, DIMENSION(:), ALLOCATABLE :: index_result, a_index, b_index
        INTEGER :: auxsize, counter, i
        CHARACTER(len = 2) :: typesearch
        
        res = 0
        if (ALLOCATED(nlc)) then
            if (size(nlc) .GT. 0 ) then
                if (size(nlc) .eq. size(c_U)) then
                    ALLOCATE(index_result(size(nlc)))
                    typesearch = "LT"
                    index_result = compare_multiple_vector(nlc,c_L, auxsize, typesearch)
                    ALLOCATE(a_index(auxsize))
                    call index_of_ones(index_result, a_index)
                    CALL reajust_index(a_index)
                    DEALLOCATE(index_result)
                   
                    ALLOCATE(index_result(size(nlc)))
                    typesearch = "GT"
                    index_result = compare_multiple_vector(nlc,c_U,auxsize,typesearch)
                    ALLOCATE(b_index(auxsize))
                    call index_of_ones(index_result, b_index)
                    CALL reajust_index(b_index)
                    DEALLOCATE(index_result)       
                    ALLOCATE(P(size(a_index)+size(b_index)))
                    P=REAL(0.0d0 ,KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
                    counter = 1
                    
                    do i=1, size(a_index)
                       P(counter)=(c_L(a_index(i))-nlc(a_index(i)))
                       counter=counter+1
                    end do
                    
                    do i=1, size(b_index)
                       P(counter)=(nlc(b_index(i))-c_U(b_index(i)))
                       counter=counter+1
                    end do
                
                end if
            end if
            
            if ( ALLOCATED(P) .and. ( maxval(P) .GT. tolc ) ) then
                res = maxval(P)
            else 
                res = 0
            end if
           
   
            if (ALLOCATED(P) ) DEALLOCATE(P)
            if (ALLOCATED(a_index) ) DEALLOCATE(a_index)
            if (ALLOCATED(b_index) ) DEALLOCATE(b_index)
        end if
            
    END FUNCTION ssm_penalty_function    
    
    
! ----------------------------------------------------------------------
! SUBROUTINES SSM_ROUND_INT
! ----------------------------------------------------------------------
    SUBROUTINE ssm_round_int(x, index1, x_L, x_U)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: x
        INTEGER, INTENT(IN) :: index1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: x_L, x_U
        INTEGER, DIMENSION(:), ALLOCATABLE :: index_vector
        INTEGER :: index2, i, DIMEN(2)
        
        DIMEN = shape(x)
        index2 = DIMEN(1)-index1+1
        ALLOCATE(index_vector(index1))
        index_vector = (/ (i, i = index2, DIMEN(1) ) /)  
        
        do i=1,DIMEN(2)
            x(index_vector,i) = x_L(index_vector)+floor(0.5+x(index_vector,i)-x_L(index_vector))
        end do
    END SUBROUTINE ssm_round_int
    
    
! ----------------------------------------------------------------------
! SUBROUTINES destroy_outfuction
! ----------------------------------------------------------------------          
    SUBROUTINE destroy_outfuction(outfunct)
        TYPE(outfuction), INTENT(INOUT) :: outfunct
        
        if ( ALLOCATED(outfunct%x) ) DEALLOCATE(outfunct%x)
        if ( ALLOCATED(outfunct%nlc) ) DEALLOCATE(outfunct%nlc)
    END SUBROUTINE destroy_outfuction
    
END MODULE funcevalinterface
