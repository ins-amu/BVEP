C$TEST SVBD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE SVBD
C***********************************************************************
C
C  SECOND TEST OF THE PORT PROGRAM DSVD
C
C***********************************************************************
C TDSVD1.F--SECOND TEST PROGRAM FOR DSVD
       INTEGER NAU, NV, I, J, K
       INTEGER M, N, IWRITE, I1MACH
       DOUBLE PRECISION A(15,5), V(5,5), W(5)
       DOUBLE PRECISION U(15,5)
       DOUBLE PRECISION TMP(15,5), EPS
       DOUBLE PRECISION TEMP, D1MACH, MAX
       LOGICAL MATU, MATV
       DOUBLE PRECISION DSTAK(500)
       COMMON /CSTAK/DSTAK
C
       IWRITE = I1MACH(2)
       CALL ISTKIN(500,4)
       M = 4
       N = 3
       NAU = 15
       NV = 5
       MATU = .TRUE.
       MATV = .TRUE.
C SET UP THE RANDOM MATRIX TO BE DECOMPOSED
       A(1,1) = 1.
       A(1,2) = 2.
       A(1,3) = 3.
       A(2,1) = 3.
       A(2,2) = 4.
       A(2,3) = 7.
       A(3,1) = 5.
       A(3,2) = 6.
       A(3,3) = 11.
       A(4,1) = 7.
       A(4,2) = 8.
       A(4,3) = 15.
C CALL THE SINGULAR VALUE DECOMPOSITION PACKAGE
       CALL DSVD(M,N,A,NAU,U,MATU,W,V,NV,MATV)
       EPS = D1MACH(4) * 1.E03
C CHECK THAT THE DECOMPOSITION RETURNED IS CORRECT
C COMPUTE U*W
       MAX = 0.
       DO 10 J=1,N
          DO 10 I=1,M
             TMP(I,J) = U(I,J) * W(J)
 10    CONTINUE
C COMPUTE U*W*TRANSPOSE(V)
       DO 40 I=1,N
          DO 40 J=1,M
              DO 20 K=1,N
                  TEMP = TEMP + TMP(J,K)*V(I,K)
 20           CONTINUE
C CALCULATE INFINITY NORM
              TEMP = TEMP - A(J,I)
              IF (DABS(TEMP).GT.MAX) MAX = DABS(TEMP)
 30           TEMP = 0
 40    CONTINUE
C PRINT OUT RESULTS
       WRITE (IWRITE,50) MAX
 50    FORMAT (33H NORM OF (A - U*W*TRANSPOSE(V))- ,E15.8)
       WRITE (IWRITE,60) (W(I), I=1,N)
 60       FORMAT (3H W-,5E15.8)
       STOP
       END
