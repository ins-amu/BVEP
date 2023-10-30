C$TEST SVDA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE SVDA
C***********************************************************************
C
C  FIRST TEST OF THE PORT PROGRAM SVD
C
C***********************************************************************
C TSVD.F--FIRST TEST PROGRAM FOR SVD
       INTEGER NAU, NV, I, J, K
       INTEGER M, N, IWRITE, I1MACH
       REAL A(15,5), V(5,5), W(5)
       REAL U(15,5)
       REAL TMP, MAX, UNI
       LOGICAL MATU, MATV
       DOUBLE PRECISION DSTAK(500)
       REAL RSTAK(1000)
       COMMON /CSTAK/DSTAK
       EQUIVALENCE (DSTAK(1),RSTAK(1))
C
       CALL ISTKIN(500,4)
       IWRITE = I1MACH(2)
       M = 15
       N = 5
       NAU = 15
       NV = 5
       MATU = .TRUE.
       MATV = .TRUE.
C SET UP THE RANDOM MATRIX TO BE DECOMPOSED
       DO 10 I=1,M
          DO 10 J=1,N
             A(I,J) = UNI(0) * 10.0
 10    CONTINUE
C CALL THE SINGULAR VALUE DECOMPOSITION PACKAGE
       CALL SVD(M,N,A,NAU,U,MATU,W,V,NV,MATV)
C CHECK THAT THE DECOMPOSITION RETURNED IS CORRECT
C COMPUTE U*W*TRANSPOSE(V)
       MAX = 0.
       DO 30 I=1,N
          DO 30 J=1,M
              TMP = 0
              DO 20 K=1,N
                  TMP = TMP + U(J,K)*V(I,K)*W(K)
 20           CONTINUE
C CALCULATE INFINITY NORM
              TMP = TMP - A(J,I)
              IF (ABS(TMP).GT.MAX) MAX = ABS(TMP)
 30    CONTINUE
C PRINT OUT RESULTS
       WRITE (IWRITE,40) MAX
 40       FORMAT (33H NORM OF (A - U*W*TRANSPOSE(V)) -,E15.8)
       WRITE (IWRITE,50) (W(I), I=1,N)
 50       FORMAT (3H W-,5E15.8)
       STOP
       END
