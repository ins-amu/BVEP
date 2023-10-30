#ifdef MATLAB
#ifdef GNU

#include "fintrf.h"

MODULE matlabproblem
   USE iso_c_binding
   USE scattersearchtypes

  CONTAINS

  SUBROUTINE openmatlab(ep) 
    IMPLICIT NONE
    integer*8, INTENT(INOUT) :: ep
    integer*8 :: engOpen
    ep = engOpen('matlab -R -nodisplay -nojvm -nosplash')
    if (ep .eq. 0) then
         write(6,*) 'Can''t start MATLAB engine'
         STOP
    endif
    

  END SUBROUTINE

  SUBROUTINE closematlab(ep)
   IMPLICIT NONE
   integer*8 :: engClose
   integer*8 :: status
   integer*8, INTENT(INOUT) :: ep

   status=engClose(ep)
   IF (status .ne. 0) THEN
          WRITE(*,*) 'Failed to close MATLAB engine'
          STOP
   ELSE
          WRITE(*,*) 'Closed MATLAB engine'
   ENDIF
  
   RETURN

  END SUBROUTINE

  SUBROUTINE matlabobjfunc(x,m,res,omt,ep)
  IMPLICIT NONE                   
  TYPE(C_PTR) ,INTENT(INOUT) ::  omt
  REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)),  INTENT(INOUT) :: res 
  INTEGER, INTENT(INOUT)  :: m
  integer*8, INTENT(INOUT) :: ep
  INTEGER*8 :: n, m1, size_r, size_g, size_j
  REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), INTENT(INOUT)  :: x(m)
  REAL*8, DIMENSION(:) , ALLOCATABLE  :: res2
  REAL*8, DIMENSION(1,m) :: vars
  integer*8  :: mxCreateDoubleMatrix, mxGetN
  integer*8 :: mxGetPr
  integer*8 :: complexity = 0
  integer*8 :: engGetVariable, status
  integer*8 :: engPutVariable, engEvalString
  integer*8 :: temp
  integer*8 :: V 
  integer*8 ::  V_ptr
 
  integer*8 :: F, F_ptr
  integer*8 :: R, R_ptr
  integer*8 :: G, G_ptr
  real*8, DIMENSION(:), ALLOCATABLE :: resi, gfort
  real*8 :: dimm(2)
  integer*8 :: dimm_ptr

  n=1
  m1=m
  ALLOCATE(res2(1))
  vars(1,:)=REAL(x,8)

  V = mxCreateDoubleMatrix(n, m1, complexity)
  V_ptr = mxGetPr(V)
  CALL mxCopyReal8ToPtr(vars, V_ptr, m1)
  status = engPutVariable(ep, 'V', V)
  IF (status .ne. 0) THEN
      WRITE(6,*) 'engPutVariable failed'
      STOP
  ENDIF
  status = engEvalString(ep,"addpath(genpath('../bitbucket/benchmarks/matlab'))");
  IF (engEvalString(ep, '[y g R] = matlabCost(V);') .ne. 0) THEN
         write(6,*) 'engEvalString failed'
         STOP
  ENDIF

   F = engGetVariable(ep, 'y')
   F_ptr = mxGetPr(F)
   CALL mxCopyPtrToReal8(F_ptr, res2, n)

   R = engGetVariable(ep, 'R')
   R_ptr = mxGetPr(R)
   m1=mxGetN(R)
   size_r = m1
   if (m1 .GT. 0 ) then
         ALLOCATE(resi(m1))
         CALL mxCopyPtrToReal8(R_ptr,resi,m1)
   end if

   G = engGetVariable(ep, 'g')
   G_ptr = mxGetPr(G)
   m1=mxGetN(G)
   size_g = m1
   if (m1 .GT. 0 ) then
         ALLOCATE(gfort(m1))
         CALL mxCopyPtrToReal8(G_ptr,gfort,m1)
   end if

   CALL allocatevectorrg(omt,size_r,resi,size_g,gfort)
   if (size_r .GT. 0 ) DEALLOCATE(resi)
   if (size_g .GT. 0 ) DEALLOCATE(gfort)

   res=res2(1)
   DEALLOCATE(res2)
   CALL mxDestroyArray(V)
   CALL mxDestroyArray(R)
   CALL mxDestroyArray(G)      
   CALL mxDestroyArray(F)

   RETURN


   END SUBROUTINE matlabobjfunc

END MODULE matlabproblem

#endif
#endif
