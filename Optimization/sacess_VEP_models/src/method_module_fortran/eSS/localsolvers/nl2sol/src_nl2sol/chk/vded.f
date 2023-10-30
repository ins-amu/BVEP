C$TEST VDED
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE VDED
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM DVDSS3
C
C***********************************************************************
C  DSIN(X(1)) + X(2)**2 + 2*X(3) + 3
C
       DOUBLE PRECISION X(3), A(10,10,10), F, FPRIME(3), PT1, PT2, PT3
       DOUBLE PRECISION H1, H2, H3, SOL
       DOUBLE PRECISION DSOL1, DSOL2, DSOL3
       INTEGER I1MACH, IWRITE, N1, N2, N3, NA1, NA2
       INTEGER I,J,K
       DOUBLE PRECISION U(3),L(3),DIST(3)
C
       IWRITE = I1MACH(2)
       N1 = 10
       N2 = 10
       N3 = 10
       X(1) = 2.5D0
       X(2) = 2.5D0
       X(3) = 2.5D0
       H1 = 1.D0/ DBLE(FLOAT(N1-1))
       H2 = 1.D0/ DBLE(FLOAT(N2-1))
       H3 = 1.D0/ DBLE(FLOAT(N3-1))
       L(1) = 1.0D0
       U(1) = 3.0D0
       L(2) = 2.0D0
       U(2) = 4.0D0
       L(3) = 1.5D0
       U(3) = 3.5D0
       DIST(1) = U(1) - L(1)
       DIST(2) = U(2) - L(2)
       DIST(3) = U(3) - L(3)
C SET UP THE SPLINE COEFFICIENTS
       DO 100 I=1,N1
              PT1 = L(1) + DIST(1)*DBLE(FLOAT(I-1))*H1
              DO 100 J=1,N2
                     PT2 = L(2) + DIST(2)*DBLE(FLOAT(J-1))*H2
                     DO 100 K=1,N3
                            PT3 = L(3) + DIST(3)*DBLE(FLOAT(K-1))*H3
                            A(I,J,K) = PT1*PT2 + PT3**2
 100   CONTINUE
       NA1 = 10
       NA2 = 10
       CALL DVDSS3 (X,N1,N2,N3,U,L,A,NA1,NA2,F,FPRIME)
C CHECK SOLUTION
       SOL = X(1)*X(2) + X(3)**2
       DSOL1 = X(2)
       DSOL2 = X(1)
       DSOL3 = 2.0D0*X(3)
       WRITE (IWRITE,101)
 101    FORMAT (34H                ACTUAL    COMPUTED)
       WRITE (IWRITE,102) SOL,F
 102    FORMAT (15H          F(X)=,2D15.8)
       WRITE (IWRITE,103) DSOL1,FPRIME(1)
 103    FORMAT (15H     PARTIAL X=,2D15.8)
       WRITE (IWRITE,104) DSOL2,FPRIME(2)
 104    FORMAT (15H     PARTIAL Y=,2D15.8)
       WRITE (IWRITE,105) DSOL3,FPRIME(3)
 105    FORMAT (15H     PARTIAL Z=,2D15.8)
       STOP
       END
