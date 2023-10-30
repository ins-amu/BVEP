C$TEST SDBA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE SDBA
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM DL2SF
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /PASS/ KCOM, ITCOM, NTCOM, IACOM
      INTEGER KCOM, ITCOM, NTCOM, IACOM
      INTEGER IA, NA(2), IS(1000), IT, NT, IX
      INTEGER IUMB, IY, NX, IW, ISTKGT, IFB
      INTEGER I, K, M, NDT, IPUMB, I1MACH
      REAL XB(3), WS(1000), RS(1000)
      LOGICAL LS(1000)
      INTEGER TEMP, TEMP1, TEMP2
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      CALL ISTKIN(1000, 3)
C
C SET OUTPUT UNIT
C
      IWRITE = I1MACH(2)
C
      XB(1) = -1
      XB(2) = 0
      XB(3) = 2
      DO  9 K = 2, 4
         KCOM = K
         TEMP = 4*K+3
         DO  8 NDT = 2, TEMP
            CALL ENTER(1)
            NA(1) = (NDT+1)/2
            NA(2) = NDT+1-NA(1)
            IF (NDT .NE. 2) GOTO 1
               IT = IUMB(XB(1), XB(3), NDT, K, NT)
               GOTO  2
   1           IT = IPUMB(XB, 3, NA, K, NT)
   2        DO  7 IFB = 1, 2
               CALL ENTER(1)
               IACOM = ISTKGT(NT-K, 3)
               TEMP1 = NT-K
               DO  3 I = 1, TEMP1
                  TEMP2 = IACOM+I
                  WS(TEMP2-1) = (-1)**(IFB+1)*I+(IFB-1)*(NT-K+1)
   3              CONTINUE
               IA = ISTKGT(NT-K, 3)
               NA(1) = (NT-K+1)/2
               NA(2) = NT-K+1-NA(1)
               IF (NT-K .NE. 2) GOTO 4
                  IX = IUMB(XB(1), XB(3), NT-K, K, NX)
                  GOTO  5
   4              IX = IPUMB(XB, 3, NA, K, NX)
   5           IY = ISTKGT(2*NX, 3)
               IW = ISTKGT(2*NX, 3)
               TEMP1 = K-1
               DO  6 M = 1, K, TEMP1
                  ITCOM = IT+K-M
                  NTCOM = NT-2*(K-M)
                  CALL CHECK(K, WS(ITCOM), NTCOM, WS(IA), WS(IACOM), WS(
     1               IX), NX, WS(IY), WS(IW))
   6              CONTINUE
               CALL LEAVE
   7           CONTINUE
            CALL LEAVE
   8        CONTINUE
   9     CONTINUE
      CALL WRAPUP
      WRITE (IWRITE,  10)
  10  FORMAT (24H  DL2SF TEST IS COMPLETE)
      STOP
      END
      SUBROUTINE CHECK(K, T, NT, A, ACOM, X, NX, Y, W)
      INTEGER NT, NX
      INTEGER K
      REAL T(NT), A(1), ACOM(1), X(NX), Y(NX, 2), W(NX, 2)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IS(1000), I, I1MACH
      REAL RS(1000), WS(500), ERRBND, AMAX1, ABS, EXP
      REAL FLOAT, ERROR, R1MACH
      LOGICAL LS(1000)
      INTEGER TEMP
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      IF (NT .LE. K) RETURN
      ERROR = 0
      ERRBND = 1.0E+1*EXP(2.0E0)*FLOAT(NT*K**2)*R1MACH(4)
      CALL F(X, NX, Y, W)
      CALL DL2SF(X, Y, NX, K, T, NT, A)
      TEMP = NT-K
      DO  1 I = 1, TEMP
         ERROR = AMAX1(ERROR, ABS(A(I)-ACOM(I)))
   1     CONTINUE
      IF (ERROR .LE. ERRBND) GOTO 3
         WRITE (IWRITE,  2) ERROR
   2     FORMAT (23H DDL2SF FAILS FOR F BY , 1PE10.2)
         ERROR = ERRBND
   3  CALL FD(X, NX, 2, Y, W)
      CALL DL2SH(X, Y, NX, 2, W, K, T, NT, A)
      TEMP = NT-K
      DO  4 I = 1, TEMP
         ERROR = AMAX1(ERROR, ABS(A(I)-ACOM(I)))
   4     CONTINUE
      IF (ERROR .LE. ERRBND) GOTO 6
         WRITE (IWRITE,  5) ERROR
   5     FORMAT (24H DDL2SH FAILS FOR FD BY , 1PE10.2)
         ERROR = ERRBND
   6  CALL FW(X, NX, 1, Y, W)
      CALL DL2SW(X, Y, NX, W, K, T, NT, A)
      TEMP = NT-K
      DO  7 I = 1, TEMP
         ERROR = AMAX1(ERROR, ABS(A(I)-ACOM(I)))
   7     CONTINUE
      IF (ERROR .LE. ERRBND) GOTO 9
         WRITE (IWRITE,  8) ERROR
   8     FORMAT (24H DDL2SW FAILS FOR FW BY , 1PE10.2)
         ERROR = ERRBND
   9  CALL FWD(X, NX, 2, Y, W)
      CALL DL2SH(X, Y, NX, 2, W, K, T, NT, A)
      TEMP = NT-K
      DO  10 I = 1, TEMP
         ERROR = AMAX1(ERROR, ABS(A(I)-ACOM(I)))
  10     CONTINUE
      IF (ERROR .LE. ERRBND) GOTO 12
         WRITE (IWRITE,  11) ERROR
  11     FORMAT (25H DDL2SH FAILS FOR FWD BY , 1PE10.2)
         ERROR = ERRBND
  12  RETURN
      END
      SUBROUTINE F(X, NX, FX, WX)
      INTEGER NX
      REAL X(NX), FX(NX), WX(NX)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /PASS/ K, IT, NT, IA
      INTEGER K, IT, NT, IA
      INTEGER IS(1000)
      REAL RS(1000), WS(500)
      LOGICAL LS(1000)
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      CALL SPLNE(K, WS(IT), NT, WS(IA), X, NX, FX)
      CALL SETR(NX, 1.0E0, WX)
      RETURN
      END
      SUBROUTINE FW(X, NX, MD, FX, WX)
      INTEGER MD, NX
      REAL X(NX), FX(NX, MD), WX(NX, MD)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /PASS/ K, IT, NT, IA
      INTEGER K, IT, NT, IA
      INTEGER IS(1000), I, J
      REAL RS(1000), WS(500)
      LOGICAL LS(1000)
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      CALL SPLNE(K, WS(IT), NT, WS(IA), X, NX, FX)
      DO  2 I = 1, NX
         DO  1 J = 1, MD
            WX(I, J) = X(I)**2+1.0E0
   1        CONTINUE
   2     CONTINUE
      RETURN
      END
      SUBROUTINE FD(X, NX, MD, FX, WX)
      INTEGER MD, NX
      REAL X(NX), FX(NX, MD), WX(NX, MD)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /PASS/ K, IT, NT, IA
      INTEGER K, IT, NT, IA
      INTEGER IS(1000)
      REAL RS(1000), WS(500)
      LOGICAL LS(1000)
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      CALL SPLND(K, WS(IT), NT, WS(IA), X, NX, MD, FX)
      CALL SETR(NX*MD, 1.0E0, WX)
      RETURN
      END
      SUBROUTINE FWD(X, NX, MD, FX, WX)
      INTEGER MD, NX
      REAL X(NX), FX(NX, MD), WX(NX, MD)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /PASS/ K, IT, NT, IA
      INTEGER K, IT, NT, IA
      INTEGER IS(1000), I, J
      REAL RS(1000), WS(500)
      LOGICAL LS(1000)
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      CALL SPLND(K, WS(IT), NT, WS(IA), X, NX, MD, FX)
      DO  2 I = 1, NX
         DO  1 J = 1, MD
            WX(I, J) = X(I)**2+1.0E0
   1        CONTINUE
   2     CONTINUE
      RETURN
      END
      SUBROUTINE WRAPUP
      INTEGER LUSED, I1MACH, NERROR, ISTKST, LMAX, NERR
      INTEGER IWRITE
C TO WRAP-UP A RUN BY CHECKING THE STACK AND ERROR STATES.
C/6S
C     IF (NERROR(NERR) .NE. 0) CALL SETERR(
C    1   36HWRAPUP - AN ERROR STATE IS LEFT OVER, 36, 1, 2)
C     IF (ISTKST(1) .GT. 0) CALL SETERR(
C    1   32HWRAPUP - STUFF LEFT ON THE STACK, 32, 2, 2)
C/7S
      IF (NERROR(NERR) .NE. 0) CALL SETERR(
     1   'WRAPUP - AN ERROR STATE IS LEFT OVER', 36, 1, 2)
      IF (ISTKST(1) .GT. 0) CALL SETERR(
     1   'WRAPUP - STUFF LEFT ON THE STACK', 32, 2, 2)
C/
      LUSED = ISTKST(3)
      LMAX = ISTKST(4)
      IWRITE = I1MACH(2)
      WRITE (IWRITE,  99) LUSED, LMAX
  99  FORMAT (6H USED , I8, 3H / , I8, 22H OF THE STACK ALLOWED.)
      RETURN
      END
