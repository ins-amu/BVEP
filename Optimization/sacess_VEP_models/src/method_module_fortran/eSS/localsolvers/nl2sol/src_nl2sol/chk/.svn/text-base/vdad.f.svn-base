C$TEST VDAD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE VDAD
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM DVDSS1
C
C***********************************************************************
C  DSIN(X) + X**2 + 3
C
       DOUBLE PRECISION X, A(10), F, FPRIME, PT, H, SOL
       INTEGER I1MACH, IWRITE, N, I
       DOUBLE PRECISION U, L, DIST, DSOL
       DOUBLE PRECISION DSIN, DCOS
C
       IWRITE = I1MACH(2)
       N = 10
       X = 1.5D0
       H = 1.D0/ DBLE(FLOAT(N-1))
       L = 1.0D0
       U = 5.0D0
       DIST = U - L
C SET UP THE SPLINE COEFFICIENTS
       DO 100 I=1,N
              PT = L + DIST*DBLE(FLOAT(I-1))*H
              A(I) = DSIN(PT) + PT**2 + 3.0D0
 100   CONTINUE
       CALL DVDSS1 (X,N,U,L,A,F,FPRIME)
C CHECK SOLUTION
       SOL = DSIN(X) + X**2 + 3.0D0
       DSOL = DCOS(X) + 2.0D0*X
       WRITE (IWRITE,101)
 101    FORMAT (34H                ACTUAL    COMPUTED)
       WRITE (IWRITE,102) SOL,F
 102    FORMAT (15H          F(X)=,2D15.8)
       WRITE (IWRITE,103) DSOL,FPRIME
 103    FORMAT (15H    DERIVATIVE=,2D15.8)
       STOP
       END
