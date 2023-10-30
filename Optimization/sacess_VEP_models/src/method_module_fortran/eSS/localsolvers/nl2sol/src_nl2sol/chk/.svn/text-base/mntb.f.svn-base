C$TEST MNTB
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE MNTB
C***********************************************************************
C
C  TEST OF THE PORT MONOTONICITY ROUTINES
C
C***********************************************************************
      INTEGER IX(12), I, INC, IWRITE, I1MACH
      REAL RX(12)
      LOGICAL SMONOD, SMONOI, SMONOR, MONOD, MONOI, MONOR
      DOUBLE PRECISION DX(12)
      DATA INC/1/
      DATA IX(1)/-4/
      DATA IX(2)/-4/
      DATA IX(3)/-2/
      DATA IX(4)/0/
      DATA IX(5)/2/
      DATA IX(6)/4/
      DATA IX(7)/4/
      DATA IX(8)/2/
      DATA IX(9)/0/
      DATA IX(10)/-2/
      DATA IX(11)/-4/
      DATA IX(12)/-4/
C     THIS PROCEDURE TESTS THE MONOTONICITY ROUTINES
      IWRITE = I1MACH(2)
      WRITE (IWRITE,  1)
   1  FORMAT (47H0THIS PROCEDURE TESTS THE MONOTONICITY ROUTINES)
      IF ((.NOT. MONOI(IX, 0, INC)) .OR. (.NOT. MONOI(IX, 1, INC))
     1    .OR. (.NOT. MONOI(IX, 2, INC)) .OR. (.NOT. MONOI(IX, 3, INC))
     2    .OR. (.NOT. MONOI(IX, 4, INC)) .OR. (.NOT. MONOI(IX, 5, INC))
     3    .OR. (.NOT. MONOI(IX(2), 5, INC)) .OR. (.NOT. MONOI(IX(3), 5
     4   , INC)) .OR. MONOI(IX(4), 5, INC) .OR. MONOI(IX(5), 5, INC)
     5    .OR. (.NOT. MONOI(IX(6), 5, INC)) .OR. (.NOT. MONOI(IX(7), 5
     6   , INC)) .OR. (.NOT. MONOI(IX(8), 5, INC)) .OR. (.NOT. MONOI(IX(
     7   9), 4, INC)) .OR. (.NOT. MONOI(IX(10), 3, INC))) GOTO 3
         WRITE (IWRITE,  2)
   2     FORMAT (18H MONOI PASSED TEST)
         GOTO  5
   3     WRITE (IWRITE,  4)
   4     FORMAT (18H MONOI FAILED TEST)
   5  DO  6 I = 1, 12
         RX(I) = IX(I)
   6     CONTINUE
      IF ((.NOT. MONOR(RX, 0, INC)) .OR. (.NOT. MONOR(RX, 1, INC))
     1    .OR. (.NOT. MONOR(RX, 2, INC)) .OR. (.NOT. MONOR(RX, 3, INC))
     2    .OR. (.NOT. MONOR(RX, 4, INC)) .OR. (.NOT. MONOR(RX, 5, INC))
     3    .OR. (.NOT. MONOR(RX(2), 5, INC)) .OR. (.NOT. MONOR(RX(3), 5
     4   , INC)) .OR. MONOR(RX(4), 5, INC) .OR. MONOR(RX(5), 5, INC)
     5    .OR. (.NOT. MONOR(RX(6), 5, INC)) .OR. (.NOT. MONOR(RX(7), 5
     6   , INC)) .OR. (.NOT. MONOR(RX(8), 5, INC)) .OR. (.NOT. MONOR(RX(
     7   9), 4, INC)) .OR. (.NOT. MONOR(RX(10), 3, INC))) GOTO 8
         WRITE (IWRITE,  7)
   7     FORMAT (18H MONOR PASSED TEST)
         GOTO  10
   8     WRITE (IWRITE,  9)
   9     FORMAT (18H MONOR FAILED TEST)
  10  DO  11 I = 1, 12
         DX(I) = IX(I)
  11     CONTINUE
      IF ((.NOT. MONOD(DX, 0, INC)) .OR. (.NOT. MONOD(DX, 1, INC))
     1    .OR. (.NOT. MONOD(DX, 2, INC)) .OR. (.NOT. MONOD(DX, 3, INC))
     2    .OR. (.NOT. MONOD(DX, 4, INC)) .OR. (.NOT. MONOD(DX, 5, INC))
     3    .OR. (.NOT. MONOD(DX(2), 5, INC)) .OR. (.NOT. MONOD(DX(3), 5
     4   , INC)) .OR. MONOD(DX(4), 5, INC) .OR. MONOD(DX(5), 5, INC)
     5    .OR. (.NOT. MONOD(DX(6), 5, INC)) .OR. (.NOT. MONOD(DX(7), 5
     6   , INC)) .OR. (.NOT. MONOD(DX(8), 5, INC)) .OR. (.NOT. MONOD(DX(
     7   9), 4, INC)) .OR. (.NOT. MONOD(DX(10), 3, INC))) GOTO 13
         WRITE (IWRITE,  12)
  12     FORMAT (18H MONOD PASSED TEST)
         GOTO  15
  13     WRITE (IWRITE,  14)
  14     FORMAT (18H MONOD FAILED TEST)
  15  IF ((.NOT. SMONOI(IX, 0, INC)) .OR. (.NOT. SMONOI(IX, 1, INC))
     1    .OR. SMONOI(IX, 2, INC) .OR. SMONOI(IX, 3, INC) .OR. SMONOI(
     2   IX, 4, INC) .OR. SMONOI(IX, 5, INC) .OR. (.NOT. SMONOI(IX(2), 5
     3   , INC)) .OR. SMONOI(IX(3), 5, INC) .OR. SMONOI(IX(4), 5, INC)
     4    .OR. SMONOI(IX(5), 5, INC) .OR. SMONOI(IX(6), 5, INC) .OR. (
     5   .NOT. SMONOI(IX(7), 5, INC)) .OR. SMONOI(IX(8), 5, INC) .OR.
     6   SMONOI(IX(9), 4, INC) .OR. SMONOI(IX(10), 3, INC)) GOTO 17
         WRITE (IWRITE,  16)
  16     FORMAT (20H  SMONOI PASSED TEST)
         GOTO  19
  17     WRITE (IWRITE,  18)
  18     FORMAT (20H  SMONOI FAILED TEST)
  19  IF ((.NOT. SMONOR(RX, 0, INC)) .OR. (.NOT. SMONOR(RX, 1, INC))
     1    .OR. SMONOR(RX, 2, INC) .OR. SMONOR(RX, 3, INC) .OR. SMONOR(
     2   RX, 4, INC) .OR. SMONOR(RX, 5, INC) .OR. (.NOT. SMONOR(RX(2), 5
     3   , INC)) .OR. SMONOR(RX(3), 5, INC) .OR. SMONOR(RX(4), 5, INC)
     4    .OR. SMONOR(RX(5), 5, INC) .OR. SMONOR(RX(6), 5, INC) .OR. (
     5   .NOT. SMONOR(RX(7), 5, INC)) .OR. SMONOR(RX(8), 5, INC) .OR.
     6   SMONOR(RX(9), 4, INC) .OR. SMONOR(RX(10), 3, INC)) GOTO 21
         WRITE (IWRITE,  20)
  20     FORMAT (20H  SMONOR PASSED TEST)
         GOTO  23
  21     WRITE (IWRITE,  22)
  22     FORMAT (20H  SMONOR FAILED TEST)
  23  IF ((.NOT. SMONOD(DX, 0, INC)) .OR. (.NOT. SMONOD(DX, 1, INC))
     1    .OR. SMONOD(DX, 2, INC) .OR. SMONOD(DX, 3, INC) .OR. SMONOD(
     2   DX, 4, INC) .OR. SMONOD(DX, 5, INC) .OR. (.NOT. SMONOD(DX(2), 5
     3   , INC)) .OR. SMONOD(DX(3), 5, INC) .OR. SMONOD(DX(4), 5, INC)
     4    .OR. SMONOD(DX(5), 5, INC) .OR. SMONOD(DX(6), 5, INC) .OR. (
     5   .NOT. SMONOD(DX(7), 5, INC)) .OR. SMONOD(DX(8), 5, INC) .OR.
     6   SMONOD(DX(9), 4, INC) .OR. SMONOD(DX(10), 3, INC)) GOTO 25
         WRITE (IWRITE,  24)
  24     FORMAT (20H  SMONOD PASSED TEST)
         GOTO  27
  25     WRITE (IWRITE,  26)
  26     FORMAT (20H  SMONOD FAILED TEST)
  27  STOP
      END
