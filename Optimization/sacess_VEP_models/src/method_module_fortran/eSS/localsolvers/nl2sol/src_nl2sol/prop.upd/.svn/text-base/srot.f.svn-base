      SUBROUTINE SROT(N, SX, INCX, SY, INCY, SC, SS)
      INTEGER N, INCX, INCY
      REAL SX(INCX, 1), SY(INCY, 1), SC, SS
      INTEGER J
      REAL SW
C     APPLY  ( SC SS)  TO THE 2 BY N MATRIX (SX(1) ... SX(N))
C            (-SS SC)                       (SY(1) ... SY(N))
      IF (N .EQ. 0) RETURN
C/6S
C     IF (N .LT. 0) CALL SETERR(17H  SROT - N .LT. 0, 17, 1, 2)
C     IF (INCX .LE. 0) CALL SETERR(20H  SROT - INCX .LT. 0, 20, 1, 2)
C     IF (INCY .LE. 0) CALL SETERR(20H  SROT - INCY .LT. 0, 20, 1, 2)
C/7S
      IF (N .LT. 0) CALL SETERR('  SROT - N .LT. 0', 17, 1, 2)
      IF (INCX .LE. 0) CALL SETERR('  SROT - INCX .LT. 0', 20, 1, 2)
      IF (INCY .LE. 0) CALL SETERR('  SROT - INCY .LT. 0', 20, 1, 2)
C/
      DO  1 J = 1, N
         SW = SC*SX(1, J)+SS*SY(1, J)
         SY(1, J) = (-SS)*SX(1, J)+SC*SY(1, J)
         SX(1, J) = SW
   1     CONTINUE
      RETURN
      END
