C$TEST RNRM
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE RNRM
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM RNORM
C
C***********************************************************************
      INTEGER IWRITE,J, K
      INTEGER COUNT(80)
      REAL X, RNORM
C
C  ZERO THE COUNTERS
C
      DO 5 K=1,80
  5    COUNT(K) = 0
C
C  SET UP THE OUTPUT WRITE UNIT
C
      IWRITE = I1MACH(2)
C
      DO 10 K=1,20000
      X = RNORM(0)
      IF (X .LT. (-3.E0)) COUNT(79)=COUNT(79)+1
      IF (X .GT. 3.E0) COUNT(80)=COUNT(80)+1
      IF (ABS(X) .GT. 3.E0) GO TO 10
C
      J = (X + 4.E0)*10.E0 + 1.E0
      COUNT(J) = COUNT(J)+1
 10     CONTINUE
C
      WRITE(IWRITE,99) COUNT
 99      FORMAT(1H0,10I5)
      STOP
      END
