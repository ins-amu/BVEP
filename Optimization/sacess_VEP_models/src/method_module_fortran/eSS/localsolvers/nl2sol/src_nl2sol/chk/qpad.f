C$TEST QPAD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE QPAD
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM DIQP
C
C***********************************************************************
C TEST PROGRAM FOR QP/FEASA  2/12/81
        DOUBLE PRECISION X(50), Q(50,50), A(200,50), BL(50), BU(50)
        DOUBLE PRECISION C(50), B(200), TEMPQ(50,50), INIT(50)
        DOUBLE PRECISION SAVEQ(50,50), AX(200), SUM(50)
        DOUBLE PRECISION SOL(50), EPSI, SUM1, D1MACH
        INTEGER N, I, J, IPRINT, MAXITR, IQ, M, IA, IE, TEST(5)
        INTEGER IN, ICOUNT, NUM, IROLD
        INTEGER ITMP1, ITMP, ICT, IFLAG, IPROB, IWRITE
C
        DOUBLE PRECISION DSTAK
        COMMON /CSTAK/DSTAK(6000)
C
        IWRITE = I1MACH(2)
        CALL ISTKIN(6000,4)
        CALL ENTSRC(IROLD,1)
        EPSI = DSQRT(D1MACH(4))
C
        DO 100 IN = 5,20,15
           WRITE (IWRITE,110) IN
 110       FORMAT (/3HN  ,I5)
        DO 200 ICT=1,2
           WRITE (IWRITE,210) ICT
 210       FORMAT(/19HRANDOM GENERATION -,I2)
           IE = ICT - 1
        DO 300 IPROB = 1,2
           ICOUNT = 0
           DO 301 I=1,5
              TEST(I) = 0
 301       CONTINUE
           N = IN
           M = N
C PROBLEM -1
           IQ = 50
           IA = 200
           IPRINT = 0
           IF (N.LE.5.AND.ICT.EQ.1)IPRINT=1
           MAXITR = 150
           IF (IPROB.NE.1) GO TO 305
           WRITE (IWRITE,310)
 310       FORMAT(48HTEST -1  Q IS POSITIVE DEFINITE  A,X ARE RANDOM.)
           CALL SETUP1(X, Q, A, BL, BU, C, B, TEMPQ,
     *          SAVEQ, IQ, IA, N, M, INIT)
           IFLAG = 0
           GO TO 315
 305       WRITE (IWRITE,320)
 320       FORMAT(41HTEST -1  Q IS INDEFINITE  A,X ARE RANDOM.)
           CALL SETUP2(X, Q, A, BL, BU, C, B, TEMPQ,
     *          SAVEQ, IQ, IA, N, M)
           IFLAG = 1
C
 315       CONTINUE
           NUM = 1
           WRITE (IWRITE,321) M, IE
 321       FORMAT (13H CONSTRAINTS ,I5,10H EQUALITY ,I5)
           CALL DIQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
C ERROR IN DFEAS
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
C ERROR IN DIQP
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
           DO 330 I=1,N
               SOL(I) = X(I)
 330       CONTINUE
           CALL FUNCT(N, IQ, X, C, SUM, Q)
           CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM))
           CALL CHBND(N, BL, BU, X, TEST(NUM))
C
C          IF (N.GE.10) GO TO 340
           WRITE (IWRITE,350)
 350       FORMAT (15H FINAL SOLUTION)
           DO 360 I=1,N
              WRITE (IWRITE,370) X(I)
 370          FORMAT (F12.3)
 360       CONTINUE
C PROBLEM -2 --PUT IN THE SOLUTION FROM PROBLEM -1
 340       WRITE (IWRITE,380)
 380       FORMAT(48HTEST -2  SOLVING WITH SOLUTION AS INITIAL GUESS.)
           NUM = 2
           DO 385 I=1,N
               X(I) = SOL(I)
 385       CONTINUE
           DO 390 I=1,N
              DO 390 J=1,N
                 Q(I,J) = SAVEQ(I,J)
 390       CONTINUE
C
           CALL DIQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
           CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
           CALL FUNCT(N, IQ, X, C, SUM, Q)
           CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM))
           CALL CHBND(N, BL, BU, X, TEST(NUM))
C
C PROBLEM -3 --VIOLATE A SIMPLE CONSTRAINT
           WRITE (IWRITE,400)
 400       FORMAT(41HTEST -3  VIOLATE FIRST SIMPLE CONSTRAINT.)
           NUM = 3
           DO 405 I=1,N
               X(I) = SOL(I)
 405       CONTINUE
           X(1) = BL(1) - 10.D0
           DO 410 I=1,N
              DO 410 J=1,N
                 Q(I,J) = SAVEQ(I,J)
 410       CONTINUE
           CALL DIQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
           CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
           CALL FUNCT(N, IQ, X, C, SUM, Q)
           CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM))
           CALL CHBND(N, BL, BU, X, TEST(NUM))
C
C PROBLEM -4 --VIOLATE A GENERAL CONSTRAINT
           WRITE (IWRITE,420)
 420       FORMAT(42HTEST -4  VIOLATE FIRST GENERAL CONSTRAINT.)
           NUM = 4
           DO 425 I=1,N
               X(I) = SOL(I)
 425       CONTINUE
           SUM1 = 0.
           DO 430 I=2,N
              SUM1 = SUM1 + A(1,I)*X(I)
 430      CONTINUE
          X(1) = 10.D0 + (SUM1 + B(1))/A(1,1)
          DO 440 I=1,N
             DO 440 J=1,N
                Q(I,J) = SAVEQ(I,J)
 440      CONTINUE
          CALL DIQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
          CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
          CALL FUNCT(N, IQ, X, C, SUM, Q)
          CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM))
          CALL CHBND(N, BL, BU, X, TEST(NUM))
C
C PROBLEM -5 --MAKE SIMPLE CONSTRAINTS INTO GENERAL
          WRITE (IWRITE,450)
 450      FORMAT(46HTEST -5  MAKE SIMPLE CONSTRAINTS INTO GENERAL.)
          NUM = 5
          DO 455 I=1,N
              X(I) = INIT(I)
 455      CONTINUE
          M = M + 2*N
          DO 460 I=1,N
             ITMP = N+I
             ITMP1 = 2*N+I
             A(ITMP,I) = 1.D0
             A(ITMP1,I) = -1.D0
             B(ITMP) = BL(I)
             B(ITMP1) = -BU(I)
             BL(I) = BL(I) - 10.D0
             BU(I) = BU(I) + 10.D0
 460      CONTINUE
          DO 470 I=1,N
             DO 470 J=1,N
                Q(I,J) = SAVEQ(I,J)
 470      CONTINUE
          CALL DIQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IE)
           IF ((NERROR(NERR).EQ.10).OR.(NERROR(NERR).EQ.8)
     1        .OR.(NERROR(NERR).EQ.9)) CALL ERROFF
           IF ((NERROR(NERR).EQ.6).OR.(NERROR(NERR).EQ.7)) CALL ERROFF
          CALL CHANGE(N, SOL, X, EPSI, TEST(NUM), IFLAG)
          CALL FUNCT(N, IQ, X, C, SUM, Q)
          CALL MULT(A, X, AX, B, N, M, IA, TEST(NUM))
          CALL CHBND(N, BL, BU, X, TEST(NUM))
C
          DO 480 I=1,N
             BL(I) = BL(I) + 10.
             BU(I) = BU(I) - 10.
 480      CONTINUE
C
          DO 481 I=1,5
              IF (TEST(I).EQ.0)  GOTO 481
              WRITE (IWRITE,482) I
 482          FORMAT (25HTHE PROGRAM FAILED TEST -,I3)
              ICOUNT = ICOUNT + 1
 481      CONTINUE
          IF (ICOUNT.EQ.0) WRITE (IWRITE,490)
 490      FORMAT (25HALL TESTS WERE SUCCESSFUL/)
          IF (ICOUNT.NE.0) WRITE (IWRITE,500) ICOUNT
 500      FORMAT (19HTHE PROGRAM FAILED ,I3,6H TESTS/)
 300   CONTINUE
 200    CONTINUE
 100   CONTINUE
       STOP
       END
         SUBROUTINE CHANGE(N, SOL, X, EPSI, ICOUNT, IFLAG)
         INTEGER N, ICOUNT, I, IFLAG, IWRITE
         DOUBLE PRECISION SOL(N), X(N), RNRM, RMAX, EPSI
C
         IWRITE = I1MACH(2)
         RMAX = 0.
         DO 1 I=1,N
            RNRM = DABS(SOL(I) - X(I))
            IF (RNRM.GT.RMAX) RMAX = RNRM
 1       CONTINUE
         WRITE (IWRITE,2) RMAX
 2       FORMAT (21H CHANGE IN SOLUTION  ,E15.5)
         IF (IFLAG.EQ.1) RETURN
         IF (RMAX.GT.EPSI) ICOUNT = 1
         RETURN
         END
        SUBROUTINE FUNCT(N, IQ, X, C, SUM, Q)
        INTEGER N, IQ, ITMP, IWRITE
        DOUBLE PRECISION X(N), C(N), Q(IQ,N), SUM(N)
        DOUBLE PRECISION F, CTX
        IWRITE = I1MACH(2)
        CTX = 0.
        DO 5 I=1,N
            CTX = X(I) * C(I) + CTX
 5      CONTINUE
        DO 6 J=1,N
            SUM(J) = 0.
            SUM(1) = SUM(1) + X(J)*Q(1,J)
 6      CONTINUE
        DO 7 I=2,N
            DO 8 J=I,N
                SUM(I) = SUM(I) + X(J)*Q(I,J)
 8          CONTINUE
            ITMP = I-1
            DO 9 J=1,ITMP
                SUM(I) = SUM(I) + X(J)*Q(J,I)
 9          CONTINUE
 7      CONTINUE
        F = 0.
        DO 10 I=1,N
            F = SUM(I) * X(I) + F
 10     CONTINUE
        F = F/2. + CTX
        WRITE (IWRITE,1002) F
 1002   FORMAT (22H FINAL FUNCTION VALUE , D14.5)
        RETURN
        END
C
        SUBROUTINE MULT(A, X, AX, B, N, M, IA, TEST)
        INTEGER IA, IWRITE
        INTEGER N, M, I, J, TEST, MTEST
        DOUBLE PRECISION A(IA,N), X(N), AX(N), B(N)
        DOUBLE PRECISION  RMAX, RES, EPSI, D1MACH
C
        IWRITE = I1MACH(2)
        MTEST = 0
        RMAX = 0.
        EPSI = DSQRT(D1MACH(4))
        DO 1 I=1,M
           AX(I) = 0.
           DO 3 J=1,N
                AX(I) = AX(I) + A(I,J) * X(J)
 3         CONTINUE
           RES = B(I) - AX(I)
           IF (RES.LE.EPSI) GO TO 1
           WRITE (IWRITE,10) I, RES
 10        FORMAT (11H CONSTRAINT,I5,24H IS VIOLATED.  RESIDUAL ,E15.5)
           IF (RMAX.LT.RES) RMAX = RES
           TEST = 1
           MTEST = 1
 1      CONTINUE
        IF (MTEST.EQ.0) WRITE (IWRITE,20)
 20     FORMAT (30H NO CONSTRAINTS WERE VIOLATED.)
        IF (MTEST.GT.0) WRITE (IWRITE,30) RMAX
 30      FORMAT (38H RESIDUAL OF THE VIOLATED CONSTRAINTS ,E15.5)
        RETURN
        END
       SUBROUTINE CHBND(N, BL, BU, X, TEST)
       INTEGER N, ILOW, IUP, I, TEST, IWRITE
       DOUBLE PRECISION BL(N), BU(N), X(N)
C
       IWRITE = I1MACH(2)
       ILOW = 0
       IUP = 0
       DO 1 I=1,N
           IF (BL(I).GT.X(I)) GO TO 100
 120       IF (BU(I).LT.X(I)) GO TO 110
           GO TO 1
C
 100       WRITE (IWRITE,150) I
 150       FORMAT (13H LOWER BOUND ,I5,9H VIOLATED)
           ILOW = ILOW + 1
           GO TO 120
 110       WRITE (IWRITE,160) I
 160       FORMAT (13H UPPER BOUND ,I5,9H VIOLATED)
           IUP = IUP + 1
 1      CONTINUE
        IF ((ILOW.NE.0).OR.(IUP.NE.0)) GO TO 210
        WRITE (IWRITE,200)
 200    FORMAT (40H NO UPPER OR LOWER BOUNDS WERE VIOLATED.)
        RETURN
 210    WRITE (IWRITE,220) ILOW
 220    FORMAT (1H ,I5,27H LOWER BOUNDS WERE VIOLATED)
        WRITE (IWRITE,230) IUP
 230    FORMAT (1H ,I5,27H UPPER BOUNDS WERE VIOLATED)
        TEST = 1
        RETURN
       END
C PROGRAM TO SET UP PROBLEMS FOR TESTCOMP.F 6/7/82
       SUBROUTINE SETUP1(X, Q, A, BL, BU, C, B, TEMPQ, SAVEQ,
     *              IQ, IA, N, M, INIT)
        INTEGER N, I, J, IQ, M, IA
        DOUBLE PRECISION X(N), Q(IQ,N), A(IA,N), BL(N), BU(N)
        DOUBLE PRECISION C(N), B(M), TEMPQ(IQ,N), INIT(N)
        DOUBLE PRECISION SAVEQ(IQ,N), SUM, DBLE
C
        DO 1 I=1,N
            X(I) = DBLE(UNI(0) * 10.)
            INIT(I) = X(I)
            C(I) = 0.D0
            DO 7 J=1,I
                TEMPQ(J,I) = 0.D0
                TEMPQ(I,J) = DBLE(UNI(0) * 10. + 1.)
 7          CONTINUE
 1      CONTINUE
        DO 9 I=1,M
            B(I) = 0.D0
            DO 11 J=1,N
                A(I,J) = DBLE(UNI(0) * 10. + 1.)
                B(I) = B(I) + A(I,J)
 11         CONTINUE
        B(I) = B(I)/2.D0
 9      CONTINUE
C
        DO 12 I=1,N
            ITMP = M+I
            B(ITMP) = 0.D0
            ITMP1 = M+N+I
            B(ITMP1) = 0.D0
            DO 12 K=1,N
                A(ITMP,K) = 0.D0
                A(ITMP1,K) = 0.D0
                SUM = 0.D0
                DO 13 J=1,N
                    SUM = SUM + TEMPQ(K,J)*TEMPQ(I,J)
 13             CONTINUE
                Q(K,I) = SUM
                SAVEQ(K,I) = SUM
 12     CONTINUE
C
        DO 3 I=1,N
            BL(I) = 0.D0
            BU(I) = 10.D0
 3      CONTINUE
        RETURN
        END
C
C
C
        SUBROUTINE SETUP2(X, Q, A, BL, BU, C, B, TEMPQ, SAVEQ,
     *              IQ, IA, N, M)
        INTEGER N, I, J, IQ, M, IA
        DOUBLE PRECISION X(N), Q(IQ,N), A(IA,N), BL(N), BU(N)
        DOUBLE PRECISION C(N), B(M), TEMPQ(IQ,N)
        DOUBLE PRECISION SAVEQ(IQ,N), DBLE
C
        DO 1973 I=1,N
           DO 1973 J=1,N
                Q(I,J) = DBLE(UNI(0) * 10. + 1.)
                SAVEQ(I,J) = Q(I,J)
 1973   CONTINUE
        RETURN
        END
