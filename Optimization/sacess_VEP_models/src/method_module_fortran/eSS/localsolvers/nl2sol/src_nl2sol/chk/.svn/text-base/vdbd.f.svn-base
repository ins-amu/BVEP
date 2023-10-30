C$TEST VDBD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE VDBD
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM DVDSS2
C
C***********************************************************************
C  DSIN(X(1)) + X(2)**2 + 3
C
       DOUBLE PRECISION X(2), A(10,10), F, FPRIME(2), PT1, PT2, H1, H2
       DOUBLE PRECISION SOL, DSOL1, DSOL2
       INTEGER I1MACH, IWRITE, N1, N2, NA1, I, J
       DOUBLE PRECISION U(2),L(2),DIST(2)
C
       DOUBLE PRECISION  DSIN, DCOS
C
C
       IWRITE = I1MACH(2)
       N1 = 10
       N2 = 10
       X(1) = 1.5D0
       X(2) = 1.5D0
       H1 = 1.D0/ DBLE(FLOAT(N1-1))
       H2 = 1.D0/ DBLE(FLOAT(N2-1))
       L(1) = 1.0D0
       U(1) = 5.0D0
       L(2) = 1.0D0
       U(2) = 5.0D0
       DIST(1) = U(1) - L(1)
       DIST(2) = U(2) - L(2)
C SET UP THE SPLINE COEFFICIENTS
       DO 100 I=1,N1
              PT1 = L(1) + DIST(1)*DBLE(FLOAT(I-1))*H1
              DO 100 J=1,N2
                     PT2 = L(2) + DIST(2)*DBLE(FLOAT(J-1))*H2
                     A(I,J) = DSIN(PT1) + PT2**2 + 3.0D0
 100   CONTINUE
       NA1 = 10
       CALL DVDSS2 (X,N1,N2,U,L,A,NA1,F,FPRIME)
C CHECK SOLUTION
       SOL = DSIN(X(1)) + X(2)**2 + 3.0D0
       DSOL1 = DCOS(X(1))
       DSOL2 = 2.0D0*X(2)
       WRITE (IWRITE,101)
 101    FORMAT (34H                ACTUAL    COMPUTED)
       WRITE (IWRITE,102) SOL,F
 102    FORMAT (15H          F(X)=,2D15.8)
       WRITE (IWRITE,103) DSOL1,FPRIME(1)
 103    FORMAT (15H     PARTIAL X=,2D15.8)
       WRITE (IWRITE,104) DSOL2,FPRIME(2)
 104    FORMAT (15H     PARTIAL Y=,2D15.8)
       STOP
       END
