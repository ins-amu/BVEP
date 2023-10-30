C$TEST SVDC
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE SVDC
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM SVDLS
C
C***********************************************************************
C TSVDLS.F--TEST PROGRAM FOR SVDLS
       INTEGER I,J,K
       INTEGER M, N, IWRITE, I1MACH
       REAL A(10,2), V(2,2), W(2)
       REAL B(10), U(10,2)
       REAL TMP, MAX
C
       IWRITE = I1MACH(2)
       M = 6
       N = 2
C SET THE FIRST COLUMN OF THE A ARRAYS TO THE ACTUAL X,
C AND THE SECOND COLUMN TO 1.0
       DO 10 I=1,M
          A(I,1) = FLOAT(I)
          A(I,2) = 1.0
 10    CONTINUE
C SET THE VALUES OF THE RIGHT HAND SIDE, B
       B(1) = .3
       B(2) = .95
       B(3) = 2.6
       B(4) = 2.5
       B(5) = 2.3
       B(6) = 3.95
C CALL THE SINGULAR VALUE DECOMPOSITION PACKAGE
       CALL SVDLS(M,N,A,10,U,.TRUE.,W,V,2,.TRUE.,B,1)
       WRITE (IWRITE,20)
 20    FORMAT (4H SVD)
       WRITE (IWRITE,30) B(1), B(2)
 30    FORMAT (7H C(1)- , E15.8, 7H C(2): , E15.8)
       WRITE (IWRITE,40) W(1), W(2)
 40    FORMAT (7H W(1)- , E15.8, 7H W(2): , E15.8)
C CHECK THAT THE DECOMPOSITION RETURNED IS CORRECT
C COMPUTE U*W*TRANSPOSE(V)
       MAX = 0.
       DO 60 I=1,N
          DO 60 J=1,M
              TMP = 0
              DO 50 K=1,N
                  TMP = TMP + U(J,K)*V(I,K)*W(K)
 50           CONTINUE
C CALCULATE INFINITY NORM
              TMP = TMP - A(J,I)
              IF (ABS(TMP).GT.MAX) MAX = ABS(TMP)
 60    CONTINUE
C PRINT OUT RESULTS
       WRITE (IWRITE,70) MAX
 70       FORMAT (33H NORM OF (A - U*W*TRANSPOSE(V)) - ,E15.8)
       STOP
       END
