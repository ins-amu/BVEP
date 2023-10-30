MODULE common_functions
    USE scattersearchtypes
    CONTAINS
    
    
    
! ----------------------------------------------------------------------
! inicializa_Vector_Value
! ----------------------------------------------------------------------
    SUBROUTINE inicializa_Vector_Value(MM, D, Value)
	! declaración de argumentos
        IMPLICIT NONE
        INTEGER, INTENT(IN):: D
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN):: Value
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT), POINTER:: MM                 
        ! declaración de variables locales
        INTEGER:: error
        ! sentencias Fortran
        ALLOCATE(MM(D), STAT=error)
        IF (error/=0) THEN
                PRINT *,"ERROR DE INICIALIZACION"
        END IF
        MM = Value
        PRINT *,"INIT"
        END SUBROUTINE inicializa_Vector_Value
        
           
     
        SUBROUTINE inicializa_Vector(MM, D)
        ! declaración de argumentos
                IMPLICIT NONE;
                INTEGER, INTENT(IN):: D
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT), ALLOCATABLE:: MM
        ! declaración de variables locales
                INTEGER:: error
        ! sentencias Fortran
                ALLOCATE(MM(D), STAT=error)
                IF (error/=0) THEN
                        PRINT *,"ERROR DE INICIALIZACION"
                END IF
        END SUBROUTINE inicializa_Vector        


! ----------------------------------------------------------------------
! destroy_Vector
! ----------------------------------------------------------------------
        SUBROUTINE destroy_Vector(MM)
        ! declaración de argumentos
                IMPLICIT NONE;
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT), ALLOCATABLE:: MM
        ! declaración de variables locales
                INTEGER:: error
        ! sentencias Fortran
                DEALLOCATE(MM, STAT=error)
                IF (error/=0) THEN
                        PRINT *,"ERROR EN LA DESTRUCCION"
                END IF
                !PRINT *,"DESTROY"
        END SUBROUTINE destroy_Vector


        SUBROUTINE slave_information(idp, dim_refset, ndiverse, MODE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: idp, dim_refset, ndiverse, MODE
            
#ifdef MPI2
            CALL mpibarrieress()            
#endif
            if (idp .NE. 0) then
                WRITE(*, '(A10, I3, A20, I4, A15, I5)') "ID SLAVE: ", idp, &
                         "Refset size:", dim_refset , "ndiverse:", ndiverse
            end if
            
            if (MODE .EQ. 1) then
                if (idp .EQ. 0) then
                WRITE(*, '(A10, I3, A20, I4, A15, I5)') "ID SLAVE: ", idp, &
                         "Refset size:", dim_refset , "ndiverse:", ndiverse
                end if                
            end if
            
            
#ifdef MPI2
            CALL mpibarrieress()      
            ! To order the ouput messages
            ! CALL sleepmpi(1)
#endif              
        END SUBROUTINE slave_information
        
        
        
        SUBROUTINE optimization_begins(idp) 
            INTEGER, INTENT(INOUT) :: idp
            if (idp .EQ. 0 ) then
                print *, ""
                WRITE(*, '(A51)') "************ OPTIMIZATION BEGINS ******************"
                print *, ""
            end if
        END SUBROUTINE optimization_begins
        
        
        
        SUBROUTINE ending_solver_message(idp,fin) 
            INTEGER, INTENT(IN) :: idp
            INTEGER, INTENT(IN) :: fin

            if (idp .EQ. 0) then
                print *, ""
                WRITE(*, '(A51)') "************ RESULTS ******************************"
                print *, ""
                
                   if (fin .eq. 1) then
                       WRITE(*, '(A67)') "[PROGRAM FINISHED] Maximum number of function evaluations achieved."
                   else if ( fin .eq. 2) then
                       WRITE(*, '(A55)') "[PROGRAM FINISHED] Maximum computation time achieved."                    
                   else if ( fin .eq. 3) then     
                       WRITE(*, '(A51)') "[PROGRAM FINISHED] Desired function value achieved."
                   else if ( fin .eq. 4) then
                       WRITE(*, '(A52)') "[PROGRAM FINISHED] Optimization stopped by the user."
                   else if ( fin .eq. 5) then
                       WRITE(*, '(A80)') "[PROGRAM FINISHED] Maximum number of iterations without significant improvement."
                   end if
                print *, ""
            endif
            
            
        END SUBROUTINE ending_solver_message        
        
        
        
        SUBROUTINE seed_recount(exp1,id)
            INTEGER(C_INT) :: returnseedcounter, counterseed
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: seed
            INTEGER, INTENT(IN) :: id
            TYPE(C_PTR), INTENT(INOUT) :: exp1
                    
            counterseed = returnseedcounter(exp1)
            ALLOCATE( seed(counterseed))
            CALL returnseed(exp1, seed)
            WRITE(*, '(A2, I4, A25, F25.1)') 'ID' , id, ' -- SEED RANDOM NUMBERS:', seed(1)
            

            DEALLOCATE( seed)
        
        
        END SUBROUTINE seed_recount
! ----------------------------------------------------------------------
! print_Vector
! ----------------------------------------------------------------------
        SUBROUTINE print_Vector(MM)
        ! declaración de argumentos
                IMPLICIT NONE;
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE:: MM
        ! declaración de variables locales
                INTEGER:: i
        ! sentencias Fortran
                DO i=1,size(MM,1)
                        PRINT *, MM(i)
                END DO
        END SUBROUTINE print_Vector
        

! ----------------------------------------------------------------------
! random_Vector_5
! ----------------------------------------------------------------------
    SUBROUTINE random_Vector_5(exp1, X, D, Xl, Xu, cont)
    ! declaración de argumentos
        IMPLICIT NONE;
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN):: Xl, Xu
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: X
        INTEGER, INTENT(IN):: D, cont
        INTEGER :: i
    ! declaración de variables locales
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)):: random
        

        DO i=1,D
            random = RETURN_RANDOM_NUMBER(exp1)
            random = (random + REAL(cont,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) ) / 5d0
            X(i) = (Xu(i)-Xl(i))*random+Xl(i)
        END DO

        
    END SUBROUTINE random_Vector_5
        
               

! ----------------------------------------------------------------------
! all_positive
! ----------------------------------------------------------------------
    SUBROUTINE all_positive(X)
    ! declaración de argumentos
        IMPLICIT NONE;
                INTEGER :: i
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT), ALLOCATABLE:: X


        DO i=1,size(X)
             if (X(i) < 0) then
                X(i) = REAL(-1.0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) * X(i)
             end if
        END DO

    END SUBROUTINE all_positive

! ----------------------------------------------------------------------
! random_Vector_SEED
! ----------------------------------------------------------------------
        SUBROUTINE random_Vector_SEED(exp1, X, D, Xl, Xu)
        ! declaración de argumentos
                IMPLICIT NONE;
                INTEGER :: i
                TYPE(C_PTR), INTENT(INOUT) :: exp1
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN)::  Xl, Xu
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: X
                INTEGER, INTENT(IN):: D                 
        ! declaración de variables locales
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)):: random
                
                
                DO i=1,D                
                        random = RETURN_RANDOM_NUMBER(exp1)
                        X(i) = (Xu(i)-Xl(i))*random+Xl(i)
                END DO
        END SUBROUTINE random_Vector_SEED    
        
! ----------------------------------------------------------------------
! random_Vector_SEED_10
! ----------------------------------------------------------------------
        SUBROUTINE random_Vector_SEED_10(exp1,X,D)
        ! declaración de argumentos
                IMPLICIT NONE;
                INTEGER :: i
                TYPE(C_PTR), INTENT(INOUT) :: exp1
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT), ALLOCATABLE:: X        
                INTEGER, INTENT(IN) :: D
        ! declaración de variables locales
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)):: random
                
                
                DO i=1,D                
                        random = RETURN_RANDOM_NUMBER(exp1)
                        X(i) = REAL(random,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
                         
                END DO

        END SUBROUTINE random_Vector_SEED_10           
      
! ----------------------------------------------------------------------
! random_Matrix_SEED
! ----------------------------------------------------------------------
        SUBROUTINE random_Matrix_SEED_10(exp1,X,sizei,sizej)
        ! declaración de argumentos
                IMPLICIT NONE
                INTEGER :: i,j
                TYPE(C_PTR), INTENT(INOUT) :: exp1
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), INTENT(INOUT), ALLOCATABLE:: X
                INTEGER, INTENT(IN) :: sizei, sizej
        ! declaración de variables locales
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)):: random
                
                
                DO j=1,sizej        
                    DO i=1,sizei
                        random = RETURN_RANDOM_NUMBER(exp1)
                        X(i,j) = REAL(random,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
                    END DO
                END DO

        END SUBROUTINE random_Matrix_SEED_10   
               
        
! ----------------------------------------------------------------------
! RETURN_RANDOM_NUMBER
! ----------------------------------------------------------------------
        FUNCTION RETURN_RANDOM_NUMBER(exp1) RESULT (res)
            IMPLICIT NONE
            TYPE(C_PTR), INTENT(INOUT) :: exp1
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: getrngrandomserial,res
            
            res = REAL(getrngrandomserial(exp1),KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
        
        END FUNCTION
        
! ----------------------------------------------------------------------
! griewank
! ----------------------------------------------------------------------
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) FUNCTION ex1 (X)
        ! declaración de argumentos
                IMPLICIT NONE;
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), POINTER:: X
                
                ex1 = 4*x(1)*x(1)-2.1*x(1)**4+1/3*x(1)**6+x(1)*x(2)-4*x(2)*x(2)+4*x(2)**4;
                
        END FUNCTION ex1
        
        
! ----------------------------------------------------------------------
! INSERT_DINAMIC_VECTOR_DOUBLE_PRECISSION
! ----------------------------------------------------------------------
    SUBROUTINE INSERT_DINAMIC_VECTOR_DOUBLE_PRECISSION (X,currentsize)
            IMPLICIT NONE
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: X
            INTEGER, INTENT(IN) :: currentsize
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: X_temp
            
            ALLOCATE(X_temp(currentsize+1))
            X_temp(1:currentsize) = X
            !DEALLOCATE(X)
            CALL MOVE_ALLOC(X_temp,X)
            
    END SUBROUTINE INSERT_DINAMIC_VECTOR_DOUBLE_PRECISSION

! ----------------------------------------------------------------------
! FUSION_MATRIX
! ----------------------------------------------------------------------
    SUBROUTINE FUSION_MATRIX (matrix1,matrix2)
            IMPLICIT NONE
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: matrix1
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: matrix2
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: X_temp
            INTEGER :: DIM1(2), DIM2(2)

            DIM1 = shape(matrix1)
            DIM2 = shape(matrix2)

            ALLOCATE(X_temp(DIM1(1), DIM1(2)+DIM2(2)))
            X_temp(:,1:DIM1(2)) = matrix1
            X_temp(:,DIM1(2)+1:DIM1(2)+DIM2(2)) = matrix2

            CALL MOVE_ALLOC(X_temp,matrix2)

    END SUBROUTINE FUSION_MATRIX

! ----------------------------------------------------------------------
! FUSION_VECTOR
! ----------------------------------------------------------------------
    SUBROUTINE FUSION_VECTOR (vector1,vector2)
            IMPLICIT NONE
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vector1
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: vector2
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE :: X_temp
            INTEGER :: DIM1, DIM2

            DIM1 = size(vector1)
            DIM2 = size(vector2)

            ALLOCATE(X_temp(DIM1+DIM2))
            X_temp(1:DIM1) = vector1
            X_temp(DIM1+1:DIM1+DIM2) = vector2

            CALL MOVE_ALLOC(X_temp,vector2)

    END SUBROUTINE FUSION_VECTOR
    
! ----------------------------------------------------------------------
! FUSION_VECTOR_INT
! ----------------------------------------------------------------------
    SUBROUTINE FUSION_VECTOR_INT (vector1,vector2)
            IMPLICIT NONE
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vector1
            INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: vector2
            INTEGER, DIMENSION(:), ALLOCATABLE :: X_temp
            INTEGER :: DIM1, DIM2

            DIM1 = size(vector1)
            DIM2 = size(vector2)

            ALLOCATE(X_temp(DIM1+DIM2))
            X_temp(1:DIM1) = vector1
            X_temp(DIM1+1:DIM1+DIM2) = vector2

            CALL MOVE_ALLOC(X_temp,vector2)

    END SUBROUTINE FUSION_VECTOR_INT    

! ----------------------------------------------------------------------
! PRINT_MATRIX
! ----------------------------------------------------------------------
    SUBROUTINE PRINT_MATRIX (matrix)
            IMPLICIT NONE
            REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: matrix
            INTEGER :: DIMEN(2), i,j

            DIMEN = shape(matrix)
            do, j=1,DIMEN(2)
                write(*, *) ( matrix(i,j), i=1,DIMEN(1))
            enddo
    END SUBROUTINE PRINT_MATRIX
    
! ----------------------------------------------------------------------
! PRINT_MATRIX_INTEGER
! ----------------------------------------------------------------------
    SUBROUTINE PRINT_MATRIX_INTEGER (matrix)
            IMPLICIT NONE
            INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: matrix
            INTEGER :: DIMEN(2), i,j

            DIMEN = shape(matrix)
            do, j=1,DIMEN(2)
                write(*, *) ( matrix(i,j), i=1,DIMEN(1))
            enddo
    END SUBROUTINE PRINT_MATRIX_INTEGER

! ----------------------------------------------------------------------
! SQUARE_ROOTS_POLINOM
! ----------------------------------------------------------------------
    SUBROUTINE SQUARE_ROOTS_POLINOM (matrix1, result1)
        REAL(SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(3), INTENT(IN) :: matrix1
        REAL(SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(2), INTENT(INOUT) :: result1

        result1(1) = (-matrix1(2) + SQRT((matrix1(2)** REAL(2.0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) )&
        -4.0*matrix1(1)*matrix1(3))) / &
                REAL(2.0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) *matrix1(1)
        result1(2) = (-matrix1(2) - SQRT((matrix1(2)** REAL(2.0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) )&
        -4.0*matrix1(1)*matrix1(3))) / &
                REAL(2.0,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) *matrix1(1)
    END SUBROUTINE SQUARE_ROOTS_POLINOM


! ----------------------------------------------------------------------
! choose
! ----------------------------------------------------------------------    
    function choose(vector) result (res)

        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vector
        INTEGER, DIMENSION(:), ALLOCATABLE :: vectint
        INTEGER :: res, i, j, contador

        ALLOCATE(vectint(size(vector)))

        contador = 0
        vectint = vectint * 0
        do j = 1, size(vector)
            vectint(j) = j
            do i = 1, size(vector)
                if (findinteger(vectint, i) .eq. 0) then
                    contador = contador + 1
                end if
            end do
        end do
        
        res = contador

        DEALLOCATE(vectint)

    end function choose

! ----------------------------------------------------------------------
! findinteger
! ----------------------------------------------------------------------
    function findinteger(vector,value1) result (res)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(IN) :: vector
        INTEGER, INTENT(IN) :: value1
        INTEGER :: i, res, status

        if (ALLOCATED(vector)) then
            res = 0
            do i = 1, size(vector)
                if (vector(i) .eq. value1) then
                    res = 1
                end if
            end do
        else
            print *, "ERROR : array not allocated in findinteger"
            CALL EXIT(status)
        end if
    end function

! ----------------------------------------------------------------------
! find_in_vector
! ----------------------------------------------------------------------
    SUBROUTINE find_in_vector(vector,value1,vectorout,opt)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE,  INTENT(IN) :: vector
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: value1
        INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(INOUT) :: vectorout
        INTEGER, DIMENSION(:), ALLOCATABLE :: temp
        character(len=2), intent(in)  :: opt
        INTEGER :: i, contador, empty, status

        if (ALLOCATED(vector)) then
            if (.not.ALLOCATED(vectorout)) ALLOCATE(vectorout(size(vector)))
            vectorout = vectorout * 0
            contador = 1
            empty = 0
            do i = 1, size(vector)
                if (opt .eq. "EQ") then
                    if (vector(i) .eq. value1) then
                        vectorout(contador) = i
                        contador = contador + 1
                        empty = 1
                    end if
                else if (opt .eq. "LT") then
                    if (vector(i) .lt. value1) then
                    vectorout(contador) = i
                    contador = contador + 1
                    empty = 1
                    end if
                else if (opt .eq. "GT") then
                    if (vector(i) .gt. value1) then
                    vectorout(contador) = i
                    contador = contador + 1
                    empty = 1
                    end if
                else if (opt .eq. "GE") then
                    if (vector(i) .ge. value1) then
                    vectorout(contador) = i
                    contador = contador + 1
                    empty = 1
                    end if                    
                end if
            end do
            if (empty .EQ. 1) then
                ALLOCATE(temp(contador-1))
                temp = vectorout(1:contador-1)
                CALL MOVE_ALLOC(temp, vectorout)
            else
                DEALLOCATE(vectorout)
            end if
        else
            print *, "ERROR : array not allocated in find_in_vector"
            CALL EXIT(status)
        end if

    END SUBROUTINE find_in_vector
    
    
! ----------------------------------------------------------------------
! find_in_vector_int
! ----------------------------------------------------------------------
    SUBROUTINE find_in_vector_int(vector,value1,vectorout,opt)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(IN) :: vector
        INTEGER, INTENT(IN) :: value1
        INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(INOUT) :: vectorout
        INTEGER, DIMENSION(:), ALLOCATABLE :: temp
        character(len=2), intent(in)  :: opt
        INTEGER :: i, contador, empty, status

        if (ALLOCATED(vector)) then
            if (.not.ALLOCATED(vectorout)) ALLOCATE(vectorout(size(vector)))
            vectorout = vectorout * 0
            contador = 1
            empty = 0
            do i = 1, size(vector)
                if (opt .eq. "EQ") then
                    if (vector(i) .eq. value1) then
                        vectorout(contador) = i
                        contador = contador + 1
                        empty = 1
                    end if
                else if (opt .eq. "LT") then
                    if (vector(i) .LT. value1) then
                    vectorout(contador) = i
                    contador = contador + 1
                    empty = 1
                    end if
                else if (opt .eq. "GT") then
                    if (vector(i) .GT. value1) then
                    vectorout(contador) = i
                    contador = contador + 1
                    empty = 1
                    end if
                else if (opt .eq. "GE") then
                    if (vector(i) .GE. value1) then
                    vectorout(contador) = i
                    contador = contador + 1
                    empty = 1
                    end if                    
                end if
            end do
            if (empty .EQ. 1) then
                ALLOCATE(temp(contador-1))
                temp = vectorout(1:contador-1)
                CALL MOVE_ALLOC(temp, vectorout)
            else
                DEALLOCATE(vectorout)
            end if
        else
            print *, "ERROR : array not allocated in find_in_vector_int"
            CALL EXIT(status)
        end if

    END SUBROUTINE find_in_vector_int
        
    
! ----------------------------------------------------------------------
! replace_zeros
! ----------------------------------------------------------------------
    SUBROUTINE replace_matrix(matrix,value1)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE,  INTENT(INOUT) :: matrix
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN)  :: value1
        INTEGER :: i,j, DIMEN(2), status

        if (ALLOCATED(matrix)) then
            DIMEN = shape(matrix)

            do j = 1, DIMEN(2)
                do i = 1, DIMEN(1)
                    if (matrix(i, j) .eq. REAL(0.0, KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))) then
                        matrix(i, j) = REAL(value1, KIND = SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
                    end if
                end do
            end do

        else
            print *, "ERROR : array not allocated in replace_matrix"
            CALL EXIT(status)
        end if
        

    END SUBROUTINE replace_matrix    

! ----------------------------------------------------------------------
! nchoosek_fpar
! ----------------------------------------------------------------------
    SUBROUTINE nchoosek_fpar (vector, pars, vectorout)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(IN) :: vector
        INTEGER, INTENT(IN) :: pars
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: vectorout
        INTEGER, DIMENSION(:), ALLOCATABLE :: vectint
        INTEGER :: sizer, i, j, contador, status

        if (ALLOCATED(vector)) then
            
            sizer = choose(vector)
            ALLOCATE(vectorout(sizer,pars))
            ALLOCATE(vectint(size(vector)))

            contador = 1
            vectint = vectint * 0
            do j = 1, size(vector)
                vectint(j) = j
                do i = 1, size(vector)
                    if (findinteger(vectint, i) .eq. 0) then
                        vectorout(contador,1) = vector(j)
                        vectorout(contador,2) = vector(i)
                        contador = contador + 1
                    end if
                end do
            end do

            DEALLOCATE(vectint)

        else
            print *, "ERROR : array not allocated in nchoosek_fpar"
            CALL EXIT(status)
        end if
        
    END SUBROUTINE nchoosek_fpar

! ----------------------------------------------------------------------
! DELETE_VALUES_VECTOR_INT
! ----------------------------------------------------------------------
    SUBROUTINE delete_values_vector_INT (vector1,num)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: vector1
        INTEGER, DIMENSION(:), ALLOCATABLE :: temp
        INTEGER, INTENT(IN) :: num
        INTEGER ::  status

        if (ALLOCATED(vector1)) then
            ALLOCATE(temp(size(vector1) - num))
            temp = vector1(num + 1:size(vector1))
            CALL MOVE_ALLOC(temp, vector1)
        else
            print *, "ERROR : array not allocated in delete_values_vector_INT"
            CALL EXIT(status)
        end if

    END SUBROUTINE delete_values_vector_INT

! ----------------------------------------------------------------------
! RANDPERM_INT
! ----------------------------------------------------------------------
    SUBROUTINE randperm_int (exp1, vector1)
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: vector1
        INTEGER :: i,valor, exit1, status
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: random

        if (ALLOCATED(vector1))then
            
            vector1 = 0 * vector1
            
            
            do i = 1, size(vector1)
                exit1 = 0
                do while (exit1 .eq. 0)
                    !CALL RANDOM_NUMBER(random)
                    random = RETURN_RANDOM_NUMBER(exp1)
                    valor = MOD(CEILING(random*(size(vector1)+1)), size(vector1)+1)
                    if (findinteger(vector1, valor) .eq. 0) then
                        vector1(i) = valor
                        exit1 = 1
                    end if
                end do
            end do
            
        else
            print *, "ERROR : array not allocated in randperm_int"
            CALL EXIT(status)
        end if

    END SUBROUTINE randperm_int



! ----------------------------------------------------------------------
! VARIANCE_VECTOR
! ----------------------------------------------------------------------
    FUNCTION variance_vector (vector1) result(variance)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vector1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: variance, avg
        INTEGER :: status

        if (ALLOCATED(vector1)) then
            avg = REAL(SUM(vector1)/size(vector1), KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))
            variance = SUM((vector1 - avg)**2)/size(vector1)
        else
            print *, "ERROR : array not allocated in variance_vector"
            CALL EXIT(status)
        end if        
        
    END FUNCTION variance_vector
    
    
! ----------------------------------------------------------------------
! abs_matrix
! ----------------------------------------------------------------------
    FUNCTION abs_matrix (vector1) result(res)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: vector1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: res
        INTEGER :: i,j, DIMEN(2), status
        
        if (ALLOCATED(vector1)) then
            DIMEN = shape(vector1)
            ALLOCATE(res(DIMEN(1), DIMEN(2)))
            i = 1
            j = 1
            res = reshape((/ ((abs(vector1(i, j)), i = 1, DIMEN(1)), j = 1, DIMEN(2)) /), (/ DIMEN(1), DIMEN(2) /))
        else
            print *, "ERROR : array not allocated in abs_matrix"
            CALL EXIT(status)
        end if
        
    END FUNCTION abs_matrix    
    
! ----------------------------------------------------------------------
! MAX_MATRIX
! ----------------------------------------------------------------------
    FUNCTION max_matrix (matrix1,matrix2) result(res)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: matrix1, matrix2
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE :: res
        INTEGER :: i,j, DIMEN(2), status
        
        DIMEN = shape(matrix1)
        ALLOCATE(res(DIMEN(1),DIMEN(2)))
        
        if (ALLOCATED(matrix1) .and. ALLOCATED(matrix2)) then
            do j = 1, DIMEN(2)
                do i = 1, DIMEN(1)
                    if (matrix1(i, j) .GE. matrix2(i, j)) then
                        res(i, j) = matrix1(i, j)
                    else
                        res(i, j) = matrix2(i, j)
                    end if
                end do
            end do

        else
            print *, "ERROR : array not allocated in max_matrix"
            CALL EXIT(status)
        end if
              
        
    END FUNCTION max_matrix      
    
! ----------------------------------------------------------------------
! compare_vector
! ----------------------------------------------------------------------
    FUNCTION compare_vector (vector1,value1,opt, sizer) result(res1)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vector1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: value1
        INTEGER, INTENT(INOUT) :: sizer
        character(len=2), intent(in)  :: opt
        
        INTEGER, DIMENSION(:), ALLOCATABLE :: res1
        INTEGER :: i, status
        
        if (ALLOCATED(vector1)) then
            ALLOCATE(res1(size(vector1)))
            res1 = res1 * 0
            sizer = 0
            do i = 1, size(vector1)
                if (opt .eq. "LT") then
                    if (vector1(i) .LT. value1) then
                        res1(i) = 1
                        sizer = sizer + 1
                    end if
                else if (opt .eq. "GT") then
                    if (vector1(i) .GT. value1) then
                    res1(i) = 1
                    sizer = sizer + 1
                    end if
                else if (opt .eq. "EQ") then
                    res1(i) = 1
                    sizer = sizer + 1
                end if
            end do
            
            !if (sizer .eq. 0 ) DEALLOCATE(res1)
        else
            print *, "ERROR : array not allocated in compare_vector"
            CALL EXIT(status)
        end if
      
    END FUNCTION compare_vector   
    
! ----------------------------------------------------------------------
! compare_matrix
! ----------------------------------------------------------------------
    FUNCTION compare_matrix (matrix1,value1,opt) result(res1)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: matrix1
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(IN) :: value1
        character(len=2), intent(in)  :: opt
        
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: res1
        INTEGER :: i,j, DIMEN(2), status, size
        
        
        if (ALLOCATED(matrix1)) then
            DIMEN = shape(matrix1)
            ALLOCATE(res1(DIMEN(1), DIMEN(2)))
            res1 = res1 * 0
            size = 0
            do j = 1, DIMEN(2)
                do i = 1, DIMEN(1)
                    if (opt .eq. "LT") then
                        if (matrix1(i, j) .LT. value1) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    else if (opt .eq. "GT") then
                        if (matrix1(i, j) .GT. value1) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    else if (opt .eq. "EQ") then
                        if (matrix1(i, j) .EQ. value1) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    end if
                end do
            end do
            
            !if ( size .eq. 0 ) DEALLOCATE(res1)
            
        else
            print *, "ERROR : array not allocated in compare_matrix"
            CALL EXIT(status)
        end if
      
    END FUNCTION compare_matrix     

! ----------------------------------------------------------------------
! compare_multiple_matrix
! ----------------------------------------------------------------------
    FUNCTION compare_multiple_matrix (matrix1,matrix2, size, opt) result(res1)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: matrix1, matrix2
        INTEGER, INTENT(INOUT) :: size
        character(len=2), intent(in)  :: opt
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: res1
        INTEGER :: i,j, DIMEN(2), status
        
        
            
        if (ALLOCATED(matrix1) .and. ALLOCATED(matrix2)) then
            DIMEN = shape(matrix1)
            ALLOCATE(res1(DIMEN(1), DIMEN(2)))
            res1 = res1 * 0        
            size = 0
            do j = 1, DIMEN(2)
                do i = 1, DIMEN(1)
                    if (opt .eq. "LT") then
                        if (matrix1(i, j) .LT. matrix2(i, j)) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    else if (opt .eq. "GT") then
                        if (matrix1(i, j) .GT. matrix2(i, j)) then
                        res1(i, j) = 1
                        size = size + 1
                        end if
                    else if (opt .eq. "EQ") then
                        if (matrix1(i, j) .EQ. matrix2(i, j)) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    else if (opt .eq. "NE") then
                        if (matrix1(i, j) .NE. matrix2(i, j)) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    end if                    
                end do
                
            end do
            
            !if ( size .eq. 0 ) DEALLOCATE(res1)
        else
            print *, "ERROR : array not allocated in compare_multiple_vector"
            CALL EXIT(status)
        end if
      
    END FUNCTION compare_multiple_matrix     
    
    
    
    
! ----------------------------------------------------------------------
! compare_multiple_matrix
! ----------------------------------------------------------------------
    FUNCTION compare_multiple_matrix_int (matrix1,matrix2, size, opt) result(res1)
        IMPLICIT NONE
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: matrix1, matrix2
        INTEGER, INTENT(INOUT) :: size
        character(len=2), intent(in)  :: opt
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: res1
        INTEGER :: i,j, DIMEN(2), status
        
        
            
        if (ALLOCATED(matrix1) .and. ALLOCATED(matrix2)) then
            DIMEN = shape(matrix1)
            ALLOCATE(res1(DIMEN(1), DIMEN(2)))
            res1 = 0        
            size = 0
            do j = 1, DIMEN(2)
                do i = 1, DIMEN(1)
                    if (opt .eq. "LT") then
                        if (matrix1(i, j) .LT. matrix2(i, j)) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    else if (opt .eq. "GT") then
                        if (matrix1(i, j) .GT. matrix2(i, j)) then
                        res1(i, j) = 1
                        size = size + 1
                        end if
                    else if (opt .eq. "EQ") then
                        if (matrix1(i, j) .EQ. matrix2(i, j)) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    else if (opt .eq. "EB") then
                        if ((matrix1(i, j) .EQ. matrix2(i, j)) .AND. (matrix1(i, j) .EQ. 1) ) then
                            res1(i, j) = 1
                            size = size + 1
                        end if  
                    else if (opt .eq. "NE") then
                        if (matrix1(i, j) .NE. matrix2(i, j)) then
                            res1(i, j) = 1
                            size = size + 1
                        end if
                    end if                    
                end do
                
            end do
            
            !if ( size .eq. 0 ) DEALLOCATE(res1)
        else
            print *, "ERROR : array not allocated in compare_multiple_vector"
            CALL EXIT(status)
        end if
      
    END FUNCTION compare_multiple_matrix_int    
        
    
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! compare_multiple_vector
! ----------------------------------------------------------------------
    FUNCTION compare_multiple_vector(vector1, vector2, sizer, opt) result(res1)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vector1, vector2
        INTEGER, INTENT(INOUT) :: sizer
        character(len = 2), intent(in) :: opt

        INTEGER, DIMENSION(:), ALLOCATABLE :: res1
        INTEGER :: i, status

        ALLOCATE(res1(size(vector1)))
        res1 = res1 * 0
        if (ALLOCATED(vector1) .and. ALLOCATED(vector2)) then
            sizer = 0
            do i = 1, size(vector1)
                if (opt .eq. "LT") then
                    if (vector1(i) .LT. vector2(i)) then
                        res1(i) = 1
                        sizer = sizer + 1
                    end if
                else if (opt .eq. "GT") then
                    if (vector1(i) .GT. vector2(i)) then
                    res1(i) = 1
                    sizer = sizer + 1
                    end if
                else if (opt .eq. "EQ") then
                    if (vector1(i) .EQ. vector2(i)) then
                    res1(i) = 1
                    sizer = sizer + 1
                    end if
                end if
            end do
            
            !if (sizer .eq. 0 ) DEALLOCATE(res1)
        else
            print *, "ERROR : array not allocated in compare_multiple_vector"
            CALL EXIT(status)
        end if

    END FUNCTION compare_multiple_vector

! ----------------------------------------------------------------------
! row_one
! ----------------------------------------------------------------------
    SUBROUTINE row_one (res1,res2,sizer)
        IMPLICIT NONE
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: res1
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: res2
        INTEGER, INTENT(INOUT) :: sizer
        INTEGER :: status
        
        INTEGER :: i,j, DIMEN(2), cont

        if (ALLOCATED(res1) .and. ALLOCATED(res2)) then
            DIMEN = shape(res1)
            res2 = res2 * 0
            sizer = 0
            do j = 1, DIMEN(2)
                cont = 0
                do i = 1, DIMEN(1)
                    if (res1(i, j) .eq. 1) then
                        cont = cont + 1
                    end if
                end do

                if (cont .eq. DIMEN(1)) then
                    res2(j) = 1
                    sizer = sizer + 1
                end if
            end do
        else
            print *, "ERROR : array not allocated in row_one"
            CALL EXIT(status)
        end if
      
    END SUBROUTINE row_one    
        
! ----------------------------------------------------------------------
! index_of_ones
! ----------------------------------------------------------------------
    SUBROUTINE index_of_ones (res1,BBB)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: res1
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: BBB
        INTEGER :: i, cont, status
        
        
        if (ALLOCATED(res1)) then
            BBB = BBB * 0
            cont = 1
            do i = 1, size(res1)
                if (res1(i) .eq. 1) then
                    BBB(cont) = i
                    cont = cont + 1
                end if
            end do
        else
            print *, "ERROR : array not allocated in index_of_ones"
            CALL EXIT(status)
        end if
      
    END SUBROUTINE index_of_ones
    
! ----------------------------------------------------------------------
! reajust_index
! ----------------------------------------------------------------------
    SUBROUTINE reajust_index (res1)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: res1
        INTEGER, DIMENSION(:), ALLOCATABLE :: outv
        INTEGER :: i , sizer, status
        
        if (ALLOCATED(res1)) then
            sizer = 0
            do i = 1, size(res1)
                if (res1(i) > 0) then
                    sizer = sizer + 1
                end if
            end do

            if (sizer > 0) then
                ALLOCATE(outv(sizer))
                outv = res1(1:sizer)
                DEALLOCATE(res1)
                ALLOCATE(res1(sizer))
                CALL MOVE_ALLOC(outv,res1)
            else 
                DEALLOCATE(res1)
            end if
        else
            print *, "ERROR : array not allocated in reajust_index"
            CALL EXIT(status)
        end if
      
    END SUBROUTINE reajust_index    
    
     
    
! ----------------------------------------------------------------------
! index_of_ones_matrix
! ----------------------------------------------------------------------
    SUBROUTINE index_of_ones_matrix (res1,indexvector)
        IMPLICIT NONE
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: res1
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: indexvector
        INTEGER :: status
        INTEGER :: i,j, cont, cont2, DIMEN(2)
        
        if (ALLOCATED(res1) .and. ALLOCATED(indexvector)) then
            cont = 1
            cont2 = 1
            DIMEN = shape(res1)
            indexvector = indexvector * 0
            do i = 1, DIMEN(1)
                do j = 1, DIMEN(2)
                    if (res1(i, j) .eq. 1) then
                        indexvector(cont) = cont2
                        cont = cont + 1
                    end if
                    cont2 = cont2 + 1
                end do
            end do
        else
            print *, "ERROR : array not allocated in index_of_ones_matrix"
            CALL EXIT(status)
        end if
        
    END SUBROUTINE index_of_ones_matrix    
    


! ----------------------------------------------------------------------
! return_values_matrix_index
! ----------------------------------------------------------------------

    SUBROUTINE return_values_matrix_index (vx,hyperx,aaa)
        IMPLICIT NONE
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: hyperx
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: vx
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: aaa
        INTEGER :: DIMEN(2), i, j, status
        
        DIMEN  = shape(hyperx)
        
        if ( ALLOCATED(vx) .and. ALLOCATED(hyperx) .and. ALLOCATED(aaa) ) then
        do i=1,DIMEN(1)
            do j=1,DIMEN(2)
                if (aaa(i,j) .EQ. 1 ) then
                    vx(i,j) = hyperx(i,j)
                end if
            end do
        end do
        else 
            print  *,"ERROR : array not allocated in return_values_matrix_index"
            CALL EXIT(status)
        end if
        
        
        
    END SUBROUTINE return_values_matrix_index
    
    
! ----------------------------------------------------------------------
! set_one_matrix
! ----------------------------------------------------------------------    
    SUBROUTINE set_one_matrix(matrix,vector, onesm) 
        ! declaración de argumentos
                IMPLICIT NONE;
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), INTENT(INOUT), ALLOCATABLE:: matrix
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT), ALLOCATABLE:: vector
                INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: onesm
                INTEGER :: DIMEN(2), contvect, i, j
                
                DIMEN = SHAPE(matrix)
                contvect = 1
                matrix = 0
                DO i=1,DIMEN(1)
                    DO j=1,DIMEN(2)
                        IF (onesm(i,j) .EQ. 1) THEN
                            matrix(i,j) = vector(contvect)
                            contvect = contvect  +  1
                        END IF
                    END DO
                END DO
        
    END SUBROUTINE set_one_matrix
    
    
! ----------------------------------------------------------------------
! random_det_Vector_SEED_10
! ----------------------------------------------------------------------
        SUBROUTINE random_det_Vector_SEED_10(X,D,contador,M)
        ! declaración de argumentos
                IMPLICIT NONE;
                INTEGER :: i
                INTEGER, INTENT(INOUT) :: contador
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), INTENT(INOUT), ALLOCATABLE:: X        
                INTEGER, INTENT(IN) :: D
        ! declaración de variables locales
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: M
                

        
                DO i=1,D          
                    if (contador .GE. 10000) then
                        contador = 1
                    end if
                    X(i) =  M(contador)
                        contador = contador + 1
                END DO

        END SUBROUTINE random_det_Vector_SEED_10           
      
        
! ----------------------------------------------------------------------
! random_Vector_SEED
! ----------------------------------------------------------------------
        SUBROUTINE random_det_Vector_SEED(X, D, Xl, Xu, contador,M)
        ! declaración de argumentos
                IMPLICIT NONE;
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN):: Xl, Xu
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: X
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: M
                
                INTEGER, INTENT(IN):: D                 
                INTEGER, INTENT(INOUT) :: contador
                ! declaración de variables locales
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)):: random
                
                INTEGER :: i
                

        
                DO i=1,D                
                     if (contador .GE. 10000) then
                        contador = 1
                     end if
                    random = M(contador)
                    contador = contador +1
                        X(i) = (Xu(i)-Xl(i))*random+Xl(i) 
                END DO

        END SUBROUTINE random_det_Vector_SEED   
        
! ----------------------------------------------------------------------
! random_Vector_5
! ----------------------------------------------------------------------
    SUBROUTINE random_det_Vector_5(X, D, Xl, Xu, cont, contador, M)
    ! declaración de argumentos
        IMPLICIT NONE;
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN):: Xl, Xu
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: X
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: M
        INTEGER, INTENT(INOUT) :: contador          
        INTEGER, INTENT(IN):: D, cont
        INTEGER :: i
    ! declaración de variables locales
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)):: random      
        

        
        DO i=1,D
            if (contador .GE. 10000) then
                contador = 1
            end if
            !CALL RANDOM_NUMBER(random)
            random = (M(contador) + REAL(cont,SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D))-1 ) / 5
            contador = contador + 1
            X(i) = (Xu(i)-Xl(i))*random+Xl(i)
        END DO

    END SUBROUTINE random_det_Vector_5
    
! ----------------------------------------------------------------------
! random_Matrix_SEED
! ----------------------------------------------------------------------
        SUBROUTINE random_det_Matrix_SEED_10(X,sizei,sizej, contador,M)
        ! declaración de argumentos
                IMPLICIT NONE;
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:,:), INTENT(INOUT), ALLOCATABLE:: X
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: M
                INTEGER, INTENT(INOUT) :: contador  
                INTEGER, INTENT(IN) :: sizei, sizej
        ! declaración de variables locales
                INTEGER :: i,j
                

                
                DO j=1,sizej        
                    DO i=1,sizei
                        if (contador .GE. 10000) then
                            contador = 1
                        end if
                        X(i,j) = M(contador)
                        contador = contador + 1
                    END DO
                END DO

        END SUBROUTINE random_det_Matrix_SEED_10       
        
! ----------------------------------------------------------------------
! random_Matrix_SEED
! ----------------------------------------------------------------------
        SUBROUTINE random_det_NUMBER(num,contador,M)
        ! declaración de argumentos
                IMPLICIT NONE;
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT) :: num
        ! declaración de variables locales
                REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: M
                INTEGER, INTENT(INOUT) :: contador  
                
                if (contador .GE. 10000) then
                    contador = 1
                end if
                        
                
                num = M(contador)
                contador = contador + 1

        END SUBROUTINE random_det_NUMBER        
    
    
END MODULE common_functions    
