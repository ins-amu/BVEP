C$TEST QPRA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE QPRA
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM IQP
C
C***********************************************************************
C TEST PROGRAM FOR QP/FEASA  2/12/81
        REAL X(50), Q(50,50), A(200,50), BL(50), BU(50)
        REAL C(50), B(200), TEMPQ(50,50), INIT(50)
        REAL SAVEQ(50,50), AX(200), SUM(50)
        REAL SOL(50), EPSI, SUM1, R1MACH
        INTEGER N, I, J, IPRINT, MAXITR, IQ, M, IA, IE, TEST(5)
        INTEGER IN, ICOUNT, NUM, IROLD
        INTEGER ITMP1, ITMP, ICT, IFLAG, IPROB, IWRITE
C
        DOUBLE PRECISION DSTAK(6000)
        COMMON /CSTAK/DSTAK
C
        IWRITE = I1MACH(2)
        CALL ISTKIN(6000,4)
        CALL ENTSRC(IROLD,1)
        EPSI = SQRT(R1MACH(4))
C
        DO 350 IN = 5,20,15
           WRITE (IWRITE,10) IN
 10        FORMAT (/3HN  ,I5)
        DO 340 ICT=1,2
           WRITE (IWRITE,20) ICT
 20        FORMAT(/19HRANDOM GENERATION -,I2)
           IE = ICT - 1
        DO 330 IPROB = 1,2
           ICOUNT = 0
           DO 30 I=1,5
              TEST(I) = 0
 30        CONTINUE
           N = IN
           M = N
C PROBLEM -1
           IQ = 50
           IA = 200
           IPRINT = 1
           IF (IN.GT.5.OR.ICT.GT.1)IPRINT=0
           MAXITR = 150
           IF (IPROB.NE.1) GO TO 50
           WRITE (IWRITE,40)
 40        FORMAT(48HTEST -1  Q IS POSITIVE DEFINITE  A,X ARE RANDOM.)
           CALL SETUP1(X, Q, A, BL, BU, C, B, TEMPQ,
     *          SAVEQ, IQ, IA, N, M, INIT)
           IFLAG = 0
           GO TO 70
 50        WRITE (IWRITE,60)
 60        FORMAT(41HTEST -1  Q IS INDEFINITE  A,X ARE RANDOM.)
           CALL SETUP2(X, Q, A, BL, BU, C, B, TEMPQ,
     *          SAVEQ, IQ, IA, N, M)
           IFLAG = 1
C
 70        CONTINUE
           NUM = 1
           WRITE (IWRITE,80) M, IE
 80        FORMAT (13H CONSTRAINTS ,I5,10H EQUALITY ,I5)
           CALL IQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
C ERROR IN FEAS
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
C ERROR IN IQP
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
           DO 90 I=1,N
               SOL(I) = X(I)
 90        CONTINUE
           CALL FUNCT(N, IQ, X, C, SUM, Q)
           CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM), EPSI)
           CALL CHBND(N, BL, BU, X, TEST(NUM), EPSI)
C
C          IF (N.GE.10) GO TO 340
           WRITE (IWRITE,100)
 100       FORMAT (15H FINAL SOLUTION)
           DO 120 I=1,N
              WRITE (IWRITE,110) X(I)
 110          FORMAT (F12.3)
 120       CONTINUE
C PROBLEM -2 --PUT IN THE SOLUTION FROM PROBLEM -1
 130       WRITE (IWRITE,140)
 140       FORMAT(48HTEST -2  SOLVING WITH SOLUTION AS INITIAL GUESS.)
           NUM = 2
           DO 150 I=1,N
               X(I) = SOL(I)
 150       CONTINUE
           DO 160 I=1,N
              DO 160 J=1,N
                 Q(I,J) = SAVEQ(I,J)
 160       CONTINUE
C
           CALL IQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
           CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
           CALL FUNCT(N, IQ, X, C, SUM, Q)
           CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM), EPSI)
           CALL CHBND(N, BL, BU, X, TEST(NUM), EPSI)
C
C PROBLEM -3 --VIOLATE A SIMPLE CONSTRAINT
           WRITE (IWRITE,170)
 170       FORMAT(41HTEST -3  VIOLATE FIRST SIMPLE CONSTRAINT.)
           NUM = 3
           DO 180 I=1,N
               X(I) = SOL(I)
 180       CONTINUE
           X(1) = BL(1) - 10.E0
           DO 190 I=1,N
              DO 190 J=1,N
                 Q(I,J) = SAVEQ(I,J)
 190       CONTINUE
           CALL IQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
           CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
           CALL FUNCT(N, IQ, X, C, SUM, Q)
           CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM), EPSI)
           CALL CHBND(N, BL, BU, X, TEST(NUM), EPSI)
C
C PROBLEM -4 --VIOLATE A GENERAL CONSTRAINT
           WRITE (IWRITE,200)
 200       FORMAT(42HTEST -4  VIOLATE FIRST GENERAL CONSTRAINT.)
           NUM = 4
           DO 210 I=1,N
               X(I) = SOL(I)
 210       CONTINUE
           SUM1 = 0.
           DO 220 I=2,N
              SUM1 = SUM1 + A(1,I)*X(I)
 220      CONTINUE
          X(1) = 10.E0 + (SUM1 + B(1))/A(1,1)
          DO 230 I=1,N
             DO 230 J=1,N
                Q(I,J) = SAVEQ(I,J)
 230      CONTINUE
          CALL IQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
          CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
          CALL FUNCT(N, IQ, X, C, SUM, Q)
          CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM), EPSI)
          CALL CHBND(N, BL, BU, X, TEST(NUM), EPSI)
C
C PROBLEM -5 --MAKE SIMPLE CONSTRAINTS INTO GENERAL
          WRITE (IWRITE,240)
 240      FORMAT(46HTEST -5  MAKE SIMPLE CONSTRAINTS INTO GENERAL.)
          NUM = 5
          DO 250 I=1,N
              X(I) = INIT(I)
 250      CONTINUE
          M = M + 2*N
          DO 260 I=1,N
             ITMP = N+I
             ITMP1 = 2*N+I
             A(ITMP,I) = 1.E0
             A(ITMP1,I) = -1.E0
             B(ITMP) = BL(I)
             B(ITMP1) = -BU(I)
             BL(I) = BL(I) - 10.E0
             BU(I) = BU(I) + 10.E0
 260      CONTINUE
          DO 270 I=1,N
             DO 270 J=1,N
                Q(I,J) = SAVEQ(I,J)
 270      CONTINUE
          CALL IQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
          CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
          CALL FUNCT(N, IQ, X, C, SUM, Q)
          CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM), EPSI)
          CALL CHBND(N, BL, BU, X, TEST(NUM), EPSI)
C
          DO 280 I=1,N
             BL(I) = BL(I) + 10.
             BU(I) = BU(I) - 10.
 280      CONTINUE
C
          DO 300 I=1,5
              IF (TEST(I).EQ.0)  GOTO 300
              WRITE (IWRITE,290) I
 290          FORMAT (25HTHE PROGRAM FAILED TEST -,I3)
              ICOUNT = ICOUNT + 1
 300      CONTINUE
          IF (ICOUNT.EQ.0) WRITE (IWRITE,310)
 310      FORMAT (25HALL TESTS WERE SUCCESSFUL/)
          IF (ICOUNT.NE.0) WRITE (IWRITE,320) ICOUNT
 320      FORMAT (19HTHE PROGRAM FAILED ,I3,6H TESTS/)
 330   CONTINUE
 340   CONTINUE
 350   CONTINUE
       STOP
       END
C PROGRAM TO SET UP PROBLEMS FOR TESTCOMP.F 6/7/82
       SUBROUTINE SETUP1(X, Q, A, BL, BU, C, B, TEMPQ, SAVEQ,
     *              IQ, IA, N, M, INIT)
        INTEGER N, I, J, IQ, M, IA
        REAL X(N), Q(IQ,N), A(IA,N), BL(N), BU(N)
        REAL C(N), B(M), TEMPQ(IQ,N), INIT(N)
        REAL SAVEQ(IQ,N), SUM
C
        DO 20 I=1,N
            X(I) = UNI(0) * 10.
            INIT(I) = X(I)
            C(I) = 0.E0
            DO 10 J=1,I
                TEMPQ(J,I) = 0.E0
                TEMPQ(I,J) = UNI(0) * 10. + 1.
 10         CONTINUE
 20     CONTINUE
        DO 40 I=1,M
            B(I) = 0.E0
            DO 30 J=1,N
                A(I,J) = UNI(0) * 10. + 1.
                B(I) = B(I) + A(I,J)
 30         CONTINUE
        B(I) = B(I)/2.E0
 40     CONTINUE
C
        DO 60 I=1,N
            ITMP = M+I
            B(ITMP) = 0.E0
            ITMP1 = M+N+I
            B(ITMP1) = 0.E0
            DO 60 K=1,N
                A(ITMP,K) = 0.E0
                A(ITMP1,K) = 0.E0
                SUM = 0.E0
                DO 50 J=1,N
                    SUM = SUM + TEMPQ(K,J)*TEMPQ(I,J)
 50             CONTINUE
                Q(K,I) = SUM
                SAVEQ(K,I) = SUM
 60     CONTINUE
C
        DO 70 I=1,N
            BL(I) = 0.E0
            BU(I) = 10.E0
 70     CONTINUE
        RETURN
        END
C
C
C
        SUBROUTINE SETUP2(X, Q, A, BL, BU, C, B, TEMPQ, SAVEQ,
     *              IQ, IA, N, M)
        INTEGER N, I, J, IQ, M, IA
        REAL X(N), Q(IQ,N), A(IA,N), BL(N), BU(N)
        REAL C(N), B(M), TEMPQ(IQ,N)
        REAL SAVEQ(IQ,N)
C
        DO 10 I=1,N
           DO 10 J=1,N
                Q(I,J) = UNI(0) * 10. + 1.
                SAVEQ(I,J) = Q(I,J)
 10     CONTINUE
        RETURN
        END
         SUBROUTINE CHANGE(N, SOL, X, EPSI, ICOUNT, IFLAG)
         INTEGER N, ICOUNT, I, IFLAG, IWRITE
         REAL SOL(N), X(N), RNRM, RMAX, EPSI
C
         IWRITE = I1MACH(2)
         RMAX = 0.
         DO 10 I=1,N
            RNRM = ABS(SOL(I) - X(I))
            IF (RNRM.GT.RMAX) RMAX = RNRM
 10      CONTINUE
         WRITE (IWRITE,20) RMAX
 20      FORMAT (21H CHANGE IN SOLUTION  ,E15.5)
         IF (IFLAG.EQ.1) RETURN
         IF (RMAX.GT.EPSI) ICOUNT = 1
         RETURN
         END
        SUBROUTINE FUNCT(N, IQ, X, C, SUM, Q)
        INTEGER N, IQ, ITMP, IWRITE
        REAL X(N), C(N), Q(IQ,N), SUM(N)
        REAL F, CTX
        IWRITE = I1MACH(2)
        CTX = 0.
        DO 10 I=1,N
            CTX = X(I) * C(I) + CTX
 10     CONTINUE
        DO 20 J=1,N
            SUM(J) = 0.
            SUM(1) = SUM(1) + X(J)*Q(1,J)
 20     CONTINUE
        DO 50 I=2,N
            DO 30 J=I,N
                SUM(I) = SUM(I) + X(J)*Q(I,J)
 30         CONTINUE
            ITMP = I-1
            DO 40 J=1,ITMP
                SUM(I) = SUM(I) + X(J)*Q(J,I)
 40         CONTINUE
 50     CONTINUE
        F = 0.
        DO 60 I=1,N
            F = SUM(I) * X(I) + F
 60     CONTINUE
        F = F/2. + CTX
        WRITE (IWRITE,70) F
 70     FORMAT (22H FINAL FUNCTION VALUE , E14.5)
        RETURN
        END
C
        SUBROUTINE MULT(A, X, AX, B, N, M, IA, TEST, EPSI)
        INTEGER IA, IWRITE
        INTEGER N, M, I, J, TEST, MTEST
        REAL A(IA,N), X(N), AX(N), B(N)
        REAL  RMAX, RES, EPSI
C
        IWRITE = I1MACH(2)
        MTEST = 0
        RMAX = 0.
        DO 30 I=1,M
           AX(I) = 0.
           DO 10 J=1,N
                AX(I) = AX(I) + A(I,J) * X(J)
 10        CONTINUE
           RES = B(I) - AX(I)
           IF (RES.LE.EPSI) GO TO 30
           WRITE (IWRITE,20) I, RES
 20        FORMAT (11H CONSTRAINT,I5,24H IS VIOLATED.  RESIDUAL ,E15.5)
           IF (RMAX.LT.RES) RMAX = RES
           TEST = 1
           MTEST = 1
 30     CONTINUE
        IF (MTEST.EQ.0) WRITE (IWRITE,40)
 40     FORMAT (30H NO CONSTRAINTS WERE VIOLATED.)
        IF (MTEST.GT.0) WRITE (IWRITE,50) RMAX
 50     FORMAT (38H RESIDUAL OF THE VIOLATED CONSTRAINTS ,E15.5)
        RETURN
        END
       SUBROUTINE CHBND(N, BL, BU, X, TEST, EPSI)
       INTEGER N, ILOW, IUP, I, TEST, IWRITE
       REAL BL(N), BU(N), X(N), EPSI
C
       IWRITE = I1MACH(2)
       ILOW = 0
       IUP = 0
       DO 60 I=1,N
           IF ((BL(I)-X(I)).GT.EPSI) GO TO 20
 10        IF ((X(I)-BU(I)).GT.EPSI) GO TO 40
           GO TO 60
C
 20        WRITE (IWRITE,30) I
 30        FORMAT (13H LOWER BOUND ,I5,9H VIOLATED)
           ILOW = ILOW + 1
           GO TO 60
 40        WRITE (IWRITE,50) I
 50        FORMAT (13H UPPER BOUND ,I5,9H VIOLATED)
           IUP = IUP + 1
 60     CONTINUE
        IF ((ILOW.NE.0).OR.(IUP.NE.0)) GO TO 80
        WRITE (IWRITE,70)
 70     FORMAT (40H NO UPPER OR LOWER BOUNDS WERE VIOLATED.)
        RETURN
 80     WRITE (IWRITE,90) ILOW
 90     FORMAT (1H ,I5,27H LOWER BOUNDS WERE VIOLATED)
        WRITE (IWRITE,100) IUP
 100    FORMAT (1H ,I5,27H UPPER BOUNDS WERE VIOLATED)
        TEST = 1
        RETURN
       END
