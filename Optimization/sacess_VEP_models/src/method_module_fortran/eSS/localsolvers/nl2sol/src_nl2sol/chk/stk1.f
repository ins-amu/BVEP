C$TEST STK1
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE STK1
C **********************************************************************
C
C  FIRST STORAGE ALLOCATOR TEST
C TESTS THE STORAGE ALLOCATOR WITH DEFAULT INITIALIZATION LENGTH.
C
C **********************************************************************
C NUMBER OF OUTSTANDING ALLOCATIONS.
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IS(1000), ISTKMD, ISTKGT, ISTKQU, ISTKST, I
      INTEGER J, K, NALOCS, I1MACH
      REAL RS(1000), R1MACH
      LOGICAL LS(1000)
C/R
C     REAL CS(2,500), R2(2)
C/C
      COMPLEX CS(500), CMPLX
C/
      DOUBLE PRECISION D1MACH
      INTEGER TEMP
C/R
C     EQUIVALENCE (DS(1), CS(1,1), RS(1), IS(1), LS(1))
C/C
      EQUIVALENCE (DS(1), CS(1), RS(1), IS(1), LS(1))
C/
      NALOCS = 0
      TEMP = I1MACH(2)
      WRITE (TEMP,  1)
   1  FORMAT (
     1  61H AN ERROR BELOW INDICATES TROUBLE WITH THE STORAGE ALLOCATOR,
     2   )
      TEMP = I1MACH(2)
      WRITE (TEMP,  2)
   2  FORMAT (42H WHEN USING THE STACK WITH DEFAULT LENGTH.//)
      I = 1
         GOTO  4
   3     I = I+1
   4     K = 5
C GET THE ENTIRE STACK, I ITEMS AT A TIME.
C   DO THE ALLOCATIONS IN ORDER OF
C   COMPLEX, LONG REAL, REAL, INTEGER AND LOGICAL.
            GOTO  6
   5        K = K-1
   6        IF (K .LT. 1) GOTO  16
            IF (ISTKQU(K) .LT. I) GOTO  17
C GET ALL THE REMAINING STACK.
            J = ISTKGT(ISTKQU(K), K)
C TRUNCATE TO I ITEMS.
            J = ISTKMD(I)
            GOTO  12
C FILL THE SPACE UP ACCORDINGLY.
C/R
C  7           R2(1) = R1MACH(2)
C              R2(2) = R2(1)
C              CALL SETC(I, R2, CS(1,J))
C/C
   7           CALL SETC(I, CMPLX(R1MACH(2), R1MACH(2)), CS(J))
C/
               GOTO  13
   8           CALL SETD(I, D1MACH(2), DS(J))
               GOTO  13
   9           CALL SETR(I, R1MACH(2), RS(J))
               GOTO  13
  10           CALL SETI(I, I1MACH(9), IS(J))
               GOTO  13
  11           CALL SETL(I, .TRUE., LS(J))
               GOTO  13
  12           IF (K .EQ. 1) GOTO  11
               IF (K .EQ. 2) GOTO  10
               IF (K .EQ. 3) GOTO  9
               IF (K .EQ. 4) GOTO  8
               IF (K .EQ. 5) GOTO  7
  13        NALOCS = NALOCS+1
            IF (ISTKST(1) .EQ. NALOCS) GOTO 15
               TEMP = I1MACH(2)
               WRITE (TEMP,  14)
  14           FORMAT (24H ISTKST(1) IS INCORRECT.)
  15        CONTINUE
            GOTO  5
  16     CONTINUE
         GOTO  3
  17  IF (NALOCS .GE. 6) GOTO 19
         TEMP = I1MACH(2)
         WRITE (TEMP,  18) NALOCS
  18     FORMAT (30H THE DEFAULT STACK ONLY HOLDS , I1, 7H ITEMS,,
     1      32H IT SHOULD HOLD SEVERAL HUNDRED.)
C RELEASE THE ALLOCATIONS, ONE BY ONE.
  19  I = 1
         GOTO  21
  20     I = I+1
  21     IF (I .GT. NALOCS) GOTO  25
         IF (ISTKST(1) .LE. 0) GOTO 22
            CALL ISTKRL(1)
            GOTO  24
  22        TEMP = I1MACH(2)
            WRITE (TEMP,  23)
  23        FORMAT (
     1         50H ALLOCATOR OUT OF ALLOCATIONS BEFORE IT SHOULD BE.)
  24     CONTINUE
         GOTO  20
  25  IF (ISTKST(1) .EQ. 0) GOTO 27
         TEMP = I1MACH(2)
         WRITE (TEMP,  26)
  26     FORMAT (49H AFTER DE-ALLOCATING ALL THE ITEMS, ITEMS REMAIN.)
  27  TEMP = I1MACH(2)
      WRITE (TEMP,  28)
  28  FORMAT (
     1   39H FIRST STORAGE ALLOCATOR TEST COMPLETE.//)
      TEMP = I1MACH(2)
      WRITE (TEMP,  29)
  29  FORMAT (49H NOW FORCE AN ERROR BY REQUESTING TOO MUCH SPACE.)
      I = ISTKGT(2*ISTKQU(5)+10, 5)
      STOP
      END
      SUBROUTINE SETC(N,V,B)
C
C     SETC SETS THE N COMPLEX ITEMS IN B TO V
C
C/R
C     REAL B(2,N), V(2)
C/C
      COMPLEX B(N),V
C/
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
C/R
C       B(1,I) = V(1)
C10     B(2,I) = V(2)
C/C
 10     B(I) = V
C/
C
      RETURN
C
      END
      SUBROUTINE SETD(N,V,B)
C
C     SETD SETS THE N DOUBLE PRECISION ITEMS IN B TO V
C
      DOUBLE PRECISION B(N),V
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
 10     B(I) = V
C
      RETURN
C
      END
      SUBROUTINE SETI(N,V,B)
C
C     SETI SETS THE N INTEGER ITEMS IN B TO V
C
      INTEGER B(N),V
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
 10     B(I) = V
C
      RETURN
C
      END
      SUBROUTINE SETL(N,V,B)
C
C     SETL SETS THE N LOGICAL ITEMS IN B TO V
C
      LOGICAL B(N),V
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
 10     B(I) = V
C
      RETURN
C
      END
      SUBROUTINE SETR(N,V,B)
C
C     SETR SETS THE N REAL ITEMS IN B TO V
C
      REAL B(N),V
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
 10     B(I) = V
C
      RETURN
C
      END
