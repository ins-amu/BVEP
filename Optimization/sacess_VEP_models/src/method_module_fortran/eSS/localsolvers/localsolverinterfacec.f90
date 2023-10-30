MODULE localsolverinterfacec
    USE iso_c_binding
    USE scattersearchtypes
    USE common_functions
    USE misqp_interface
    USE qsort_module
    USE funcevalinterface
    USE dhc_mod
CONTAINS
    SUBROUTINE  calldhcinterface( exp1,D,tolerance,fitnessfunction, point, fval, numeval )
        USE scattersearchfunctions
        IMPLICIT NONE
        INTEGER(C_INT), INTENT(INOUT) :: D
        TYPE(C_PTR), INTENT(INOUT) :: exp1
        INTEGER(C_INT), INTENT(INOUT) :: tolerance
        TYPE(C_FUNPTR), INTENT(INOUT) :: fitnessfunction
        TYPE(C_PTR), INTENT(INOUT) :: point
        REAL(C_DOUBLE), INTENT(INOUT) ::fval
        INTEGER(C_LONG), INTENT(INOUT) :: numeval
        REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: x0
        TYPE(problem):: problem1
        TYPE(opts) ::  default1
        TYPE(opts)  :: opts1
        INTEGER(C_LONG) :: eval
        INTEGER :: budget
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: initsize
        REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: thres

        CALL create_new_problem(problem1)
        CALL create_new_opts(opts1)
        problem1%empty = 0
        if (.not. ALLOCATED(problem1%XU))  ALLOCATE(problem1%XU(D))
        if (.not. ALLOCATED(problem1%XL))  ALLOCATE(problem1%XL(D))
        ALLOCATE(x0(D))
        CALL returnboundscfortran(exp1,problem1%XU,problem1%XL,D)
        CALL returnvectorcfortran(point,x0,D)
        opts1%localoptions%tol = tolerance
        default1 = ssm_default()
        opts1 = ssm_optset(default1, opts1)
        opts1%localoptions%iterprint = 0
        X0 = (X0 - problem1%XL)/(problem1%XU-problem1%XL)
        initsize = 1d-1
        budget = size(x0) * 10

        if (opts1%localoptions%tol .EQ. 1) then
           thres = 1d-6
        else if (opts1%localoptions%tol .EQ. 2) then
           thres = 1d-8
        else if (opts1%localoptions%tol .EQ. 3) then
           thres = 1d-10
        end if
        CALL dhc(problem1,exp1,opts1,fitnessfunction,x0,initsize,thres,budget,eval, fval)
        numeval = numeval + eval
        CALL setvectorcfortran(point,x0,D)
        DEALLOCATE(x0)
        CALL destroy_opts(opts1)
        CALL destroy_opts(default1)
        CALL destroy_problem(problem1)
    END SUBROUTINE calldhcinterface


END MODULE localsolverinterfacec
