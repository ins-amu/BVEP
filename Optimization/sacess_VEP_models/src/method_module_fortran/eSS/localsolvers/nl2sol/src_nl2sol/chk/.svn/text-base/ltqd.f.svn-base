C$TEST LTQD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LTQD
C***********************************************************************
C
C  TEST OF THE PORT PROGRAMS DLTSQ AND DLSTSQ
C
C***********************************************************************
C TDLTSQ2.F--TEST PROGRAM TO COMPARE DLTSQ AND DLSTSQ
       INTEGER I
       INTEGER M, N, IWUNIT, I1MACH
       DOUBLE PRECISION A(10,2)
       DOUBLE PRECISION B(10), W(10)
       DOUBLE PRECISION X(10,2), Y(10), C(2)
       DOUBLE PRECISION DFLOAT
C
       IWUNIT = I1MACH(2)
       M = 6
       N = 2
C SET THE FIRST COLUMN OF THE X AND A ARRAYS TO THE ACTUAL X,
C AND THE SECOND COLUMN TO 1.0
       DO 10 I=1,M
          A(I,1) = DFLOAT(I)
          X(I,1) = A(I,1)
          A(I,2) = 1.0
          X(I,2) = A(I,2)
 10    CONTINUE
C SET THE VALUES OF THE RIGHT HAND SIDES, Y AND B
       Y(1) = .3
       Y(2) = .95
       Y(3) = 2.6
       Y(4) = 2.5
       Y(5) = 2.3
       Y(6) = 3.95
       DO 20 I=1,M
          B(I) = Y(I)
 20    CONTINUE
C CALL THE LEAST SQUARES PACKAGE
       CALL DLSTSQ(10, 2, M, N, X, Y, 1, C)
       WRITE (IWUNIT,30)
 30     FORMAT (7H DLSTSQ)
       WRITE (IWUNIT,40) C(1), C(2)
 40     FORMAT (7H C(1)- ,E15.8, 7H C(2)- , E15.8)
C CALL THE SINGULAR VALUE DECOMPOSITION PACKAGE
       CALL DLTSQ(M,N,A,10,B,1,W)
       WRITE (IWUNIT,50)
 50     FORMAT (5H DSVD)
       WRITE (IWUNIT,40) B(1), B(2)
       STOP
       END
