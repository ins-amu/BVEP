C$TEST SPLE
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE SPLE
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM SPLNE
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER NA(2), IS(1000), IT, NT, IUMB, K
      INTEGER NDT, IPUMB
      INTEGER TEMP, TEMP1, IWRITE, I1MACH
      REAL XB(3), RS(1000), WS(1000)
      LOGICAL LS(1000)
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      CALL ISTKIN(1000, 3)
C
C  SET OUTPUT UNIT
C
      IWRITE = I1MACH(2)
      XB(1) = -1
      XB(2) = 0
      XB(3) = 2
      DO  4 K = 2, 4
         TEMP = 4*K+3
         DO  3 NDT = 2, TEMP
            CALL ENTER(1)
            IT = IUMB(XB(1), XB(3), NDT, K, NT)
            CALL USPLCK(K, WS(IT), NT)
            CALL LEAVE
            CALL ENTER(1)
            IF (NDT .NE. 2) GOTO 1
               IT = IUMB(XB(1), XB(3), NDT, K, NT)
               GOTO  2
   1           NA(1) = (NDT+1)/2
               NA(2) = NDT-NA(1)+1
               IT = IPUMB(XB, 3, NA, K, NT)
   2        CALL OSPLCK(K, WS(IT), NT)
            TEMP1 = IT+K
            CALL MSPLCK(K, WS(IT), NT, WS(TEMP1-1), NT-2*(K-1))
            CALL LEAVE
   3        CONTINUE
   4     CONTINUE
      CALL WRAPUP
      TEMP = I1MACH(2)
      WRITE (IWRITE,  5)
   5  FORMAT (21H DSPLNE TEST COMPLETE)
      STOP
      END
      SUBROUTINE USPLCK(K, T, NT)
      INTEGER NT
      INTEGER K
      REAL T(NT)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IA, IS(1000), IBASIS, ISTKGT, IBAISC, IUBSIS
      REAL RS(1000), WS(500)
      LOGICAL LS(1000)
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      CALL ENTER(1)
      IA = ISTKGT(NT-K, 3)
      IBASIS = ISTKGT(K**2*(K+1), 3)
      IUBSIS = ISTKGT(K**2*(K+1), 3)
      IBAISC = ISTKGT(K**2*(K+1), 3)
      CALL U8PLCK(K, T, NT, WS(IA), WS(IBASIS), WS(IUBSIS), WS(IBAISC))
      CALL LEAVE
      RETURN
      END
      SUBROUTINE U8PLCK(K, T, NT, A, BASIS, UBASIS, BASISC)
      INTEGER NT
      INTEGER K
      REAL T(NT), A(1), BASIS(1), UBASIS(1), BASISC(1)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER ID(1), JD, JJ, IS(1000), IUMD, IXCHEK
      INTEGER IYCHEK, NXCHEK, ISTKGT, INTRVR, I, J
      INTEGER M, IID, INA, ILEFT, IAINT, IXCECK
      INTEGER ILUMB, NXCECK, IPUMD, ILUMD, MTEMP, IXFIT
      INTEGER IYFIT, ITINT, NXFIT, NTINT, IVALU
      INTEGER IUM1
      REAL AMAXDF, RS(1000), ERRBND, TEMP, WS(500), AMAX1
      REAL ABS, FLOAT, ERROR, R1MACH
      LOGICAL LS(1000)
      INTEGER TEMP1, TEMP2, TEMP3, IWRITE, I1MACH
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      ERROR = 0
      ERRBND = 1.0E+1*10.E0**(FLOAT(K)/5.0E0)*FLOAT(K**2)*R1MACH(4)*
     1   FLOAT(NT-2*K+1)**(K-1)
      CALL ENTER(1)
      M = K+1
      INA = ISTKGT(K, 2)
      DO  1 I = 1, K
         TEMP1 = INA+I
         IS(TEMP1-1) = M
   1     CONTINUE
      IID = ISTKGT(K-1, 2)
      TEMP1 = K-1
      DO  2 I = 1, TEMP1
         TEMP3 = IID+I
         IS(TEMP3-1) = I
   2     CONTINUE
      TEMP1 = NT-K
      DO  24 I = 1, TEMP1
         CALL SETR(NT-K, 0.0E0, A)
         A(I) = 1
         CALL ENTER(1)
         ITINT = ILUMB(T, NT, 2, K+1, NTINT)
         IAINT = ISTKGT(NTINT-(K+1), 3)
         IXFIT = ILUMD(T, NT, K+1, NXFIT)
         IYFIT = ISTKGT(NXFIT, 3)
         CALL SPLNI(K, T, NT, A, WS(IXFIT), NXFIT, WS(IYFIT))
         ERROR = AMAX1(ERROR, ABS(WS(IYFIT)))
         IF (ERROR .LE. ERRBND) GOTO 4
            IWRITE = I1MACH(2)
            WRITE (IWRITE,  3) ERROR
   3        FORMAT (32H DSPLNI(T(1)) DIFFERS FROM 0 BY , 1PE10.2)
            ERROR = ERRBND
   4     CALL DL2SF(WS(IXFIT), WS(IYFIT), NXFIT, K+1, WS(ITINT), NTINT
     1      , WS(IAINT))
         CALL SPLNE(K, T, NT, A, WS(IXFIT), NXFIT, WS(IYFIT))
         IVALU = ISTKGT(NXFIT, 3)
         ID(1) = 1
         CALL SPLN1(K+1, WS(ITINT), NTINT, WS(IAINT), WS(IXFIT), NXFIT
     1      , ID, 1, WS(IVALU))
         ERROR = AMAX1(ERROR, AMAXDF(WS(IYFIT), 1, WS(IVALU), 1, NXFIT))
         IF (ERROR .LE. ERRBND) GOTO 6
            WRITE (IWRITE,  5) ERROR
   5        FORMAT (
     1         49H THE DERIVATIVE OF DSPLNI DIFFERS FROM DSPLNE BY , 1P
     2         E10.2)
            ERROR = ERRBND
   6     CALL LEAVE
         CALL ENTER(1)
         TEMP3 = I+K
         ILEFT = INTRVR(NT, T, T(TEMP3-1))
         IXCHEK = ILUMD(T(ILEFT), 2, M, NXCHEK)
         IYCHEK = ISTKGT(NXCHEK, 3)
         CALL SPLNI(K, T, NT, A, WS(IXCHEK), NXCHEK, WS(IYCHEK))
         IVALU = ISTKGT(K*NXCHEK, 3)
         CALL BSPLI(K, T, NT, WS(IXCHEK), NXCHEK, ILEFT, WS(IVALU))
         TEMP3 = IVALU+(I+K-ILEFT-1)*NXCHEK
         ERROR = AMAX1(ERROR, AMAXDF(WS(IYCHEK), 1, WS(TEMP3), 1,
     1      NXCHEK))
         IF (ERROR .LE. ERRBND) GOTO 8
            WRITE (IWRITE,  7) ERROR
   7        FORMAT (31H DSPLNI AND DBSPLI DISAGREE BY , 1PE10.2)
            ERROR = ERRBND
   8     CALL LEAVE
         CALL ENTER(1)
         IXCECK = IPUMD(T(I), K+1, IS(INA), NXCECK)
         CALL SPLNE(K, T, NT, A, WS(IXCECK), NXCECK, BASISC)
         CALL SPLND(K, T, NT, A, WS(IXCECK), NXCECK, K, BASIS)
         ERROR = AMAX1(ERROR, AMAXDF(BASISC, 1, BASIS, 1, NXCECK))
         IF (ERROR .LE. ERRBND) GOTO 10
            WRITE (IWRITE,  9) ERROR
   9        FORMAT (31H DSPLNE AND DSPLND DISAGREE BY , 1PE10.2)
            ERROR = ERRBND
  10     IF (I .EQ. K) CALL MOVEFR(K*NXCECK, BASIS, UBASIS)
         CALL MOVEFR(K*NXCECK, BASIS, BASISC)
         CALL SPLN1(K, T, NT, A, WS(IXCECK), NXCECK, IS(IID), K-1,
     1      BASIS)
         ERROR = AMAX1(ERROR, AMAXDF(BASISC(NXCECK+1), 1, BASIS, 1, (K-1
     1      )*NXCECK))
         IF (ERROR .LE. ERRBND) GOTO 12
            WRITE (IWRITE,  11) ERROR
  11        FORMAT (31H DSPLND AND DSPLN1 DISAGREE BY , 1PE10.2)
            ERROR = ERRBND
  12     IF (K .LE. I .AND. I .LT. NT-2*K+1) ERROR = AMAX1(ERROR,
     1      AMAXDF(BASISC, 1, UBASIS, 1, K*NXCECK))
         IF (ERROR .LE. ERRBND) GOTO 14
            WRITE (IWRITE,  13) ERROR, K, NT, I
  13        FORMAT (48H INTERIOR B SUB I DIFFER FROM CANONICAL FORM BY ,
     1         1PE10.2, 12H   K,NT,I = , I3, I3, I3)
            ERROR = ERRBND
  14     CALL LEAVE
         TEMP3 = I+K-1
         DO  23 J = I, TEMP3
            IF (T(J) .EQ. T(J+1)) GOTO  23
            CALL ENTER(1)
            IF (T(J+1) .EQ. T(NT)) IXCECK = IUMD(T(J), T(J+1), M)
            IF (T(J+1) .LT. T(NT)) IXCECK = IUM1(T(J), T(J+1), M+1, 1, 1
     1         , 0, MTEMP)
            CALL BSPLN(K, T, NT, WS(IXCECK), M, J, BASIS)
            CALL BSPLD(K, T, NT, WS(IXCECK), M, J, K, BASISC)
            ERROR = AMAX1(ERROR, AMAXDF(BASIS, 1, BASISC, 1, K*M))
            IF (ERROR .LE. ERRBND) GOTO 16
               WRITE (IWRITE,  15) ERROR
  15           FORMAT (33H DBSPLN DISAGREES WITH DBSPLD BY , 1PE10.2)
               ERROR = ERRBND
  16        CALL MOVEFR(K**2*M, BASISC, BASIS)
            CALL BSPL1(K, T, NT, WS(IXCECK), M, J, IS(IID), K-1, BASISC)
            TEMP2 = M*K
            ERROR = AMAX1(ERROR, AMAXDF(BASIS(TEMP2+1), 1, BASISC, 1, (K
     1         -1)*M*K))
            IF (ERROR .LE. ERRBND) GOTO 18
               WRITE (IWRITE,  17) ERROR
  17           FORMAT (33H DBSPLD DISAGREES WITH DBSPL1 BY , 1PE10.2)
               ERROR = ERRBND
  18        DO  20 JJ = 1, K
               CALL SETR(NT-K, 0.0E0, A)
               TEMP2 = J-K+JJ
               A(TEMP2) = 1
               DO  19 JD = 1, K
                  ID(1) = JD-1
                  CALL SPLN1(K, T, NT, A, WS(IXCECK), M, ID, 1, BASISC)
                  TEMP2 = (JJ-1)*M+1+(JD-1)*M*K
                  ERROR = AMAX1(ERROR, AMAXDF(BASIS(TEMP2), 1, BASISC, 1
     1               , M))
  19              CONTINUE
  20           CONTINUE
            IF (ERROR .LE. ERRBND) GOTO 22
               WRITE (IWRITE,  21) ERROR, K, NT, I, J
  21           FORMAT (33H DSPLN1 DISAGREES WITH DBSPLD BY , 1PE10.2,
     1            14H   K,NT,I,J = , I3, I3, I3, I3)
               ERROR = ERRBND
  22        CALL LEAVE
  23        CONTINUE
  24     CONTINUE
      DO  31 I = 1, K
         CALL ENTER(1)
         CALL SETR(NT-K, 0.0E0, A)
         A(I) = 1
         IXCECK = ILUMD(T(I), K+1, M, NXCECK)
         CALL SPLND(K, T, NT, A, WS(IXCECK), NXCECK, K, BASIS)
         J = NT-K+1-I
         CALL SETR(NT-K, 0.0E0, A)
         A(J) = 1
         IXCECK = ILUMD(T(J), K+1, M, MTEMP)
         TEMP1 = NXCECK/2
         DO  25 J = 1, TEMP1
            TEMP3 = IXCECK+J
            TEMP = WS(TEMP3-1)
            TEMP3 = IXCECK+J
            TEMP2 = IXCECK+NXCECK-J
            WS(TEMP3-1) = WS(TEMP2)
            TEMP2 = IXCECK+NXCECK-J
            WS(TEMP2) = TEMP
  25        CONTINUE
         CALL SPLND(K, T, NT, A, WS(IXCECK), NXCECK, K, BASISC)
         DO  27 J = 1, NXCECK
            DO  26 JJ = 2, K, 2
               TEMP1 = J+(JJ-1)*NXCECK
               BASISC(TEMP1) = BASISC(TEMP1)*(-1.0E0)
  26           CONTINUE
  27        CONTINUE
         TEMP1 = NXCECK+(K-1)*NXCECK
         TEMP2 = NXCECK-1+(K-1)*NXCECK
         BASIS(TEMP1) = BASIS(TEMP2)
         J = M
            GOTO  29
  28        J = J+M-1
  29        IF (J .GE. NXCECK) GOTO  30
            TEMP2 = J+(K-1)*NXCECK
            TEMP1 = J+1+(K-1)*NXCECK
            BASISC(TEMP2) = BASISC(TEMP1)
            GOTO  28
  30     ERROR = AMAX1(ERROR, AMAXDF(BASIS, 1, BASISC, 1, K*NXCECK))
         CALL LEAVE
  31     CONTINUE
      IF (ERROR .LE. ERRBND) GOTO 33
         WRITE (IWRITE,  32) ERROR
  32     FORMAT (44H THE END BASIS SPLINES ARE NOT SYMMETRIC BY , 1P
     1      E10.2)
         ERROR = ERRBND
  33  CALL LEAVE
      RETURN
      END
      SUBROUTINE OSPLCK(K, T, NT)
      INTEGER NT
      INTEGER K
      REAL T(NT)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IA, IB, ID(1), IS(1000), IUMD, IXCHEK
      INTEGER IYCHEK, NXCHEK, ISTKGT, I, J, L
      INTEGER INA, IXCECK, NXCECK, IVALU, IPUMD, ILUMD
      INTEGER MTEMP, I1MACH, IUM1, IWRITE
      REAL RS(1000), WS(500), ERRBND, VALU, AMAX1, ABS
      REAL FLOAT, EXP, ERROR, R1MACH
      LOGICAL LS(1000)
      INTEGER TEMP, TEMP1, TEMP2, TEMP3
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      ERRBND = 1.0E+1*10.0E0**(FLOAT(K)/5.0E0)*FLOAT(K**2)*R1MACH(4)*
     1   EXP(2.0E0)*FLOAT(NT-2*K+1)**(K-1)
      ERROR = 0
      CALL ENTER(1)
      IA = ISTKGT(NT-K, 3)
      CALL SETR(NT-K, 1.0E0, WS(IA))
      CALL ENTER(1)
      IXCHEK = ILUMD(T, NT, K+1, NXCHEK)
      IYCHEK = ISTKGT(NXCHEK, 3)
      CALL SPLNI(K, T, NT, WS(IA), WS(IXCHEK), NXCHEK, WS(IYCHEK))
      DO  1 I = 1, NXCHEK
         TEMP1 = IYCHEK+I
         TEMP = IXCHEK+I
         WS(TEMP1-1) = WS(TEMP1-1)-(WS(TEMP-1)-T(1))
         TEMP = IYCHEK+I
         ERROR = AMAX1(ERROR, ABS(WS(TEMP-1)))
   1     CONTINUE
      IF (ERROR .LE. ERRBND) GOTO 3
         IWRITE = I1MACH(2)
         WRITE (IWRITE,  2) ERROR
   2     FORMAT (34H DSPLNI(1) DIFFERS FROM X-T(1) BY , 1PE10.2)
         ERROR = ERRBND
   3  CALL LEAVE
      INA = ISTKGT(NT-1, 2)
      TEMP = NT-1
      DO  4 I = 1, TEMP
         TEMP1 = INA+I
         IS(TEMP1-1) = K+1
   4     CONTINUE
      IXCECK = IPUMD(T, NT, IS(INA), NXCECK)
      IVALU = ISTKGT(NXCECK, 3)
      DO  10 J = 1, K
         ID(1) = J-1
         CALL SPLN1(K, T, NT, WS(IA), WS(IXCECK), NXCECK, ID, 1, WS(
     1      IVALU))
         IF (J .NE. 1) GOTO 5
            VALU = 1
            GOTO  6
   5        VALU = 0
   6     DO  7 I = 1, NXCECK
            TEMP = IVALU+I
            ERROR = AMAX1(ERROR, ABS(WS(TEMP-1)-VALU))
   7        CONTINUE
         IF (ERROR .LE. ERRBND) GOTO 9
            WRITE (IWRITE,  8) ERROR
   8        FORMAT (35H B(1)+...+B(N-K) DIFFERS FROM 1 BY , 1PE10.2)
            ERROR = ERRBND
   9     CONTINUE
  10     CONTINUE
      CALL LEAVE
      CALL ENTER(1)
      IVALU = ISTKGT(K*(K+1), 3)
      TEMP = NT-K
      DO  20 I = K, TEMP
         CALL ENTER(1)
         IF (I .NE. NT-K) GOTO 11
            IXCECK = IUMD(T(I), T(I+1), K+1)
            GOTO  12
  11        IXCECK = IUM1(T(I), T(I+1), K+2, 1, 1, 0, MTEMP)
  12     DO  19 J = 1, K
            ID(1) = J-1
            CALL BSPL1(K, T, NT, WS(IXCECK), K+1, I, ID, 1, WS(IVALU))
            ERROR = 0
            IF (J .NE. 1) GOTO 13
               VALU = 1
               GOTO  14
  13           VALU = 0
  14        TEMP1 = K+1
            DO  16 L = 1, TEMP1
               DO  15 IB = 2, K
                  TEMP3 = IVALU+L
                  TEMP2 = IVALU+L-1+(IB-1)*(K+1)
                  WS(TEMP3-1) = WS(TEMP3-1)+WS(TEMP2)
  15              CONTINUE
               TEMP2 = IVALU+L
               ERROR = AMAX1(ERROR, ABS(WS(TEMP2-1)-VALU))
  16           CONTINUE
            IF (ERROR .LE. ERRBND) GOTO 18
               WRITE (IWRITE,  17) ERROR
  17           FORMAT (43H BASIS(1)+...+BASIS(N-K) DIFFERS FROM 1 BY ,
     1            1PE10.2)
               ERROR = ERRBND
  18        CONTINUE
  19        CONTINUE
         CALL LEAVE
  20     CONTINUE
      CALL LEAVE
      RETURN
      END
      SUBROUTINE MSPLCK(K, T1, NT1, T2, NT2)
      INTEGER NT1, NT2
      INTEGER K
      REAL T1(NT1), T2(NT2)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IA, ID(1), IS(1000), IUMD, IBASIS, IXCHEK
      INTEGER NXCHEK, ISTKGT, IYCEK1, IYCEK2, I, J
      INTEGER IID, INA, IBAISC, IXCECK, NXCECK, IPUMD
      INTEGER ILUMD, I1MACH, MAX0, IVALU1, IVALU2
      REAL AMAXDF, RS(1000), WS(500), ERRBND, AMAX1, FLOAT
      REAL EXP, ERROR, R1MACH
      LOGICAL LS(1000)
      INTEGER TEMP, TEMP1, IWRITE
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
      ERRBND = 1.0E+1*10.0E0**(FLOAT(K)/5.0E0)*FLOAT(K**2)*R1MACH(4)*
     1   EXP(2.0E0)*FLOAT(NT1-2*K+1)**(K-1)
      ERROR = 0
      IF (NT1 .LT. 3*K-1) GOTO 8
         CALL ENTER(1)
         IA = ISTKGT(NT1-K, 3)
         CALL SETR(NT1-K, 0.0E0, WS(IA))
         TEMP = NT2-K
         DO  1 I = 1, TEMP
            TEMP1 = IA+I+K
            WS(TEMP1-2) = I
   1        CONTINUE
         CALL ENTER(1)
         IXCHEK = ILUMD(T1, NT1, K+1, NXCHEK)
         IYCEK1 = ISTKGT(NXCHEK, 3)
         IYCEK2 = ISTKGT(NXCHEK, 3)
         CALL SPLNI(K, T1, NT1, WS(IA), WS(IXCHEK), NXCHEK, WS(IYCEK1))
         TEMP = IA+K
         CALL SPLNI(K, T2, NT2, WS(TEMP-1), WS(IXCHEK), NXCHEK, WS(
     1      IYCEK2))
         ERROR = AMAX1(ERROR, AMAXDF(WS(IYCEK1), 1, WS(IYCEK2), 1,
     1      NXCHEK))
         IF (ERROR .LE. ERRBND) GOTO 3
            WRITE (IWRITE,  2) ERROR
   2        FORMAT (46H DSPLNI HANDLES MULTIPLICITIES INCORRECTLY BY ,
     1         1PE10.2)
            ERROR = ERRBND
   3     CALL LEAVE
         INA = ISTKGT(NT1-1, 2)
         TEMP = NT1-1
         DO  4 I = 1, TEMP
            TEMP1 = INA+I
            IS(TEMP1-1) = K+1
   4        CONTINUE
         IXCECK = IPUMD(T1, NT1, IS(INA), NXCECK)
         IVALU1 = ISTKGT(NXCECK, 3)
         IVALU2 = ISTKGT(NXCECK, 3)
         DO  7 J = 1, K
            ID(1) = J-1
            CALL SPLN1(K, T1, NT1, WS(IA), WS(IXCECK), NXCECK, ID, 1,
     1         WS(IVALU1))
            TEMP = IA+K
            CALL SPLN1(K, T2, NT2, WS(TEMP-1), WS(IXCECK), NXCECK, ID, 1
     1         , WS(IVALU2))
            ERROR = AMAXDF(WS(IVALU1), 1, WS(IVALU2), 1, NXCECK)
            IF (ERROR .LE. ERRBND) GOTO 6
               WRITE (IWRITE,  5) ERROR
   5           FORMAT (23H F1 DIFFERS FROM F2 BY , 1PE10.2)
               ERROR = ERRBND
   6        CONTINUE
   7        CONTINUE
         CALL LEAVE
   8  IF (NT1 .GT. 2*K) RETURN
      CALL ENTER(1)
      IBASIS = ISTKGT(K**2*(K+1), 3)
      IBAISC = ISTKGT(K**2*(K+1), 3)
      IXCECK = IUMD(T1(K), T1(K+1), K+2)
      NXCECK = K+1
      IID = ISTKGT(K, 2)
      DO  9 I = 1, K
         TEMP = IID+I
         IS(TEMP-1) = I-1
   9     CONTINUE
      DO  13 I = 1, K
         IF (I .GT. 1) T1(I-1) = R1MACH(2)
         IF (I+MAX0(2*(K+1-I), K+1) .GT. NT1) GOTO 10
            TEMP = I+MAX0(2*(K+1-I), K+1)
            T1(TEMP) = -R1MACH(2)
  10     CALL BSPL1(K, T1(I), MAX0(2*(K+1-I), K+1), WS(IXCECK), NXCECK
     1      , K+1-I, IS(IID), K, WS(IBASIS))
         IF (I .EQ. 1) CALL MOVEFR(K**2*NXCECK, WS(IBASIS), WS(IBAISC))
         ERROR = AMAX1(ERROR, AMAXDF(WS(IBASIS), 1, WS(IBAISC), 1, K**2*
     1      NXCECK))
         IF (ERROR .LE. ERRBND) GOTO 12
            WRITE (IWRITE,  11) ERROR
  11        FORMAT (
     1         50H DBSPL1 HANDLES END MULTIPLICITIES INCORRECTLY BY ,
     2         1PE10.2)
            ERROR = ERRBND
  12     CONTINUE
  13     CONTINUE
      CALL LEAVE
      RETURN
      END
      REAL FUNCTION AMAXDF(A, INCA, B, INCB, NAB)
      INTEGER INCA, INCB, NAB
      REAL A(INCA, NAB), B(INCB, NAB)
      INTEGER I
      REAL AMAX1, ABS, MAXIFF
      MAXIFF = 0
      I = 1
         GOTO  2
   1     I = I+1
   2     IF (I .GT. NAB) GOTO  3
         MAXIFF = AMAX1(MAXIFF, ABS(A(1, I)-B(1, I)))
         GOTO  1
   3  AMAXDF = MAXIFF
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
