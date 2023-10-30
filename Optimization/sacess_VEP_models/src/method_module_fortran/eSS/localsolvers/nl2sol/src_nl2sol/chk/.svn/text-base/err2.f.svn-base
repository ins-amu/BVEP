C$TEST ERR2
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE ERR2
C **********************************************************************
C
C  TO TEST THE ERROR HANDLING PACKAGE
C TEST NUMBER 2: FATAL ERRORS.
C
C **********************************************************************
      INTEGER I1MACH
      INTEGER TEMP
C SET AN ERROR.
C/6S
      CALL SETERR(
     1   39HSECOND ERROR TEST - FATAL ERROR CHECKED, 39, 1, 2)
C/7S
C     CALL SETERR(
C    1   'SECOND ERROR TEST - FATAL ERROR CHECKED', 39, 1, 2)
C/
      TEMP = I1MACH(2)
      WRITE (TEMP,  1)
   1  FORMAT (39H A FATAL ERROR FAILS TO HALT EXECUTION.)
      STOP
      END
