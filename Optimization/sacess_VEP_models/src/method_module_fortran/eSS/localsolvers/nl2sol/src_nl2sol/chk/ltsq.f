C$TEST LTSQ
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LTSQC
C***********************************************************************
C
C  TEST OF THE PORT PROGRAMS LTSQ AND LSTSQ
C
C***********************************************************************
C TEST PROGRAM TO COMPARE LTSQ AND LSTSQ
       INTEGER I
       INTEGER M, N, IWRITE, I1MACH
       REAL A(10,2)
       REAL B(10), W(10)
       REAL X(10,2), Y(10), C(2)
C
       IWRITE = I1MACH(2)
       M = 6
       N = 2
C SET THE FIRST COLUMN OF THE X AND A ARRAYS TO THE ACTUAL X,
C AND THE SECOND COLUMN TO 1.0
       DO 151 I=1,M
          A(I,1) = FLOAT(I)
          X(I,1) = A(I,1)
          A(I,2) = 1.0
          X(I,2) = A(I,2)
 151   CONTINUE
C SET THE VALUES OF THE RIGHT HAND SIDES, Y AND B
       Y(1) = .3
       Y(2) = .95
       Y(3) = 2.6
       Y(4) = 2.5
       Y(5) = 2.3
       Y(6) = 3.95
       DO 145 I=1,M
          B(I) = Y(I)
 145   CONTINUE
C CALL THE LEAST SQUARES PACKAGE
       CALL LSTSQ(10, 2, M, N, X, Y, 1, C)
       WRITE (IWRITE,103)
 103    FORMAT (6H LSTSQ)
       WRITE (IWRITE,100) C(1), C(2)
 100    FORMAT (7H C(1)- ,E15.8, 7H C(2)- , E15.8)
C CALL THE SINGULAR VALUE DECOMPOSITION PACKAGE
       CALL LTSQ(M,N,A,10,B,1,W)
       WRITE (IWRITE,104)
 104    FORMAT (4H SVD)
       WRITE (IWRITE,100) B(1), B(2)
       STOP
       END
