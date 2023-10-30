C$TEST PRAC
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PRAC
C***********************************************************************
C
C  TEST OF THE PORT COMPLEX SPARSE MATRIX PACKAGE
C
C***********************************************************************
C  THIS IS THE TESTER FOR LINDA KAUFMAN'S SPARSE
C  MATRIX PACKAGE - COMPLEX
C
C  MARCH 20, 1981
C
C     MAIN PROGRAM
      INTEGER IWRITE,I1MACH
C     ALLOW 5000 UNDERFLOWS.
C
C     OUTPUT UNIT NUMBER
C
      IWRITE = I1MACH(2)
C
      CALL SGETS(IWRITE)
      STOP
      END
      SUBROUTINE SGETS(IWRITE)
C     IWRITE IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        GECE,GEFS,GEBS,GEML,SPMOR,SPFOR,SPMCE,SPFCE,SPMLE,SPFLE,SPMSL
C        SPFSL,SPMNF,SPFNF,SPFSF,SPMSF,SPMIN
C         COMPLEX VERSIONS
C
C
C     SUBROUTINES AND FUNCTIONS
C
C     PORT GECE,GEFS,GEBS,GEML,SPMOR,SPFOR,SPMCE,SPFCE,SPMLE,SPFLE
C     PORT SPMSL,SPFSL,SPMNF,SPFNF,SPMS,SPFSF.SPMIN
C     PORT UTILITIES ERROFF,ENTER,LEAVE,R1MACH
C     EXTERNAL SGEXX
C     BLAS CAXPY,SDOT,SSCAL,SCASUM
C     FORTRAN ABS,AMAX1,FLOAT,MAX0,MIN0
C
C     INTERNAL VARIABLES
C
      INTEGER IIA(101), JJA(2500)
      INTEGER JSPSAV(2500),IWORK(2800),IL(101),MR(100),MC(100),FR(100)
      INTEGER FC(100),MIC(100),JA(2500),IA(101)
      INTEGER I,IPVT(100),IPVTB(100),IQ(8),I1,I2,J
      INTEGER K,KASE,KB,KBFAIL,KOUNT,KP1,KSING,KSUSP(8)
      INTEGER L,LDA,LDAB,IWRITE,M,ML,MU,N,NM1,NPRINT
      REAL AINORM,ANORM,SMACH,COND,COND1,EN,ENORM,EPS, EPS1
      REAL ETNORM,FNI,FNORM,ONEPX,RCOND,RCONDB,RNORM
      REAL RTNORM,Q(10),QS(10),SCASUM,XNORM
      REAL AIN
      REAL EB(2),XN(2),EBM(2),EBF(2),EBM2(2),EBF2(2)
      REAL AFNORM,AMNORM,ANRM,CONDF1,CONDF,CONDM1,CONDM,GROWTH
      COMPLEX A(100,100),AS(2500),AINV(100,100),ASAVE(100,100),UL(2500)
      COMPLEX XEXACT(100,2),B(100,2),BM(100,2),BM2(100,2),BF(100,2)
      COMPLEX BF2(100,2)
      COMPLEX SPSAV(2500),Z(200),AA(2500)
      LOGICAL KBF
      EXTERNAL IIROW,AROW
      COMMON /SPA/ IIA,JJA,AA
      COMMON /CSTAK/ D(6000)
      CALL ISTKIN(6000,3)
C
      LDA = 100
      IAMAX = 2500
      MAXUL = 2500
      IWMAX = 2500
      CALL ENTER(1)
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT =3
C
      WRITE (IWRITE,480)
      WRITE (IWRITE,770)
C
      DO 10 I = 1, 8
         KSUSP(I) = 0
 10   CONTINUE
      KSING = 0
      KBFAIL = 0
C
C     SET EPS TO ROUNDING UNIT
C
      EPS =R1MACH(4)
      WRITE (IWRITE,490) EPS
      WRITE (IWRITE,470)
C
C     START MAIN LOOP
C
      KASE = 1
 20   CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL CGEXX(A,LDA,N,KASE,IWRITE)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 460
         ANORM = 0.0E0
         DO 30 J = 1, N
            ANORM = AMAX1(ANORM,SCASUM(N,A(1,J),1))
 30      CONTINUE
         WRITE (IWRITE,650) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (IWRITE,470)
            DO 40 I = 1, N
               WRITE (IWRITE,690) (A(I,J), J = 1, N)
 40         CONTINUE
            WRITE (IWRITE,470)
 50      CONTINUE
C
C        GENERATE EXACT SOLUTION
C
         XEXACT(1,1) = (1.0,1.0)
         XEXACT(1,2) = (1.0,1.0)
         XEXACT(2,2) = (1.0,1.0)
         IF (N .GE. 2) XEXACT(2,1) = (0.0,0.0)
         IF (N .LE. 2) GO TO 70
            DO 60 I = 3, N
               XEXACT(I,1) = -XEXACT(I-2,1)
                XEXACT(I,2) = (1.0,1.0)
 60         CONTINUE
 70      CONTINUE
C
C        SAVE MATRIX AND GENERATE R.H.S.
C
         DO 90 I = 1, N
            DO 80 J = 1, N
               ASAVE(I,J) = A(I,J)
 80         CONTINUE
 90      CONTINUE
         CALL CGEML(N,A,LDA,XEXACT,B)
         CALL CGEML(N,A,LDA,XEXACT(1,2),B(1,2))
         NN=LDA*2
         CALL MOVEFC(NN,B,BM)
         CALL MOVEFC(NN,B,BF)
         CALL MOVEFC(NN,B,BM2)
         CALL MOVEFC(NN,B,BF2)
      CALL CGETSP(N,A, LDA, IA,JA,AS,MAXA)
      CALL CGETSP(N,A,LDA,IIA,JJA,AA,MAXA)
         DO 100 I=1,MAXA
            SPSAV(I)=AS(I)
            JSPSAV(I)=JA(I)
 100     CONTINUE
         NP1=N+1
C
C CHECK ORDERING SUBROUTINES
C
         CALL SPMOR(N,IA,JA,MC,MIC)
         CALL SPFOR(N,IIROW,FC)
         IBAD=0
         DO 120 I=1,N
            MR(I)=MC(I)
            FR(I)=FC(I)
            IF (MC(I).EQ.FC(I)) GO TO 120
            IBAD=IBAD+1
            WRITE(IWRITE,110)I,MC(I),FC(I)
 110         FORMAT(19H ORDER DISAGREEMENT,3I5)
 120    CONTINUE
        IF (IBAD.EQ.0)WRITE(IWRITE,130)
 130     FORMAT(28H NO DISAGREEMENT IN ORDERING)
C
C        FACTOR AND ESTIMATE CONDITION
C
         CALL CGECE(N,A,LDA,IPVT,RCOND)
         IF (NERROR(INE).NE.0) CALL  ERROFF
C
C
C        FACTOR SPARSE FORM AND COMPARE
C
        CALL CSPMCE(N,MR,MC,AS,IA,JA,IAMAX,IL,ISIZE,CONDM,Z)
        IF(NERROR(INM).NE.0) CALL ERROFF
         CALL CSPFCE(N,FR,FC,AROW,IWORK,UL,MAXUL,ISIZE,CONDF,Z)
         IF (NERROR(INF).NE.0) CALL ERROFF
        WRITE(IWRITE,140)
 140     FORMAT(21H CONDITION ESTIMATION)
        WRITE(IWRITE,150)
 150     FORMAT(42H    GENERAL        SPARSE M       SPARSE F)
        WRITE(IWRITE,160)RCOND,CONDM,CONDF
 160    FORMAT(1H ,3E15.5)
C
C
C           TEST FOR SINGULARITY
C
          IF (INE+INM+INF.EQ.0)GO TO 190
            WRITE(IWRITE,170)INE,INM,INF
 170         FORMAT(13H INE, INM,INF,3I5)
 180          FORMAT(25H SPARSE ROUTINES DISAGREE,2I4)
               WRITE (IWRITE,500)
               KSING = KSING + 1
            GO TO 440
 190        CONTINUE
C
C              COMPUTE INVERSE AND COND1 = TRUE CONDITION
C
               DO 210 J = 1, N
                  DO 200 I = 1, N
                     AINV(I,J) = (0.0,0.0)
 200              CONTINUE
               AINV(J,J)=(1.0,0.0)
 210           CONTINUE
               CALL CGEFS(N,A,LDA,AINV,LDA,N,IPVT)
               CALL CGEBS(N,A,LDA,AINV,LDA,N)
               AINORM=ANRM(N,ASAVE,LDA,AINV,LDA, ANORM, COND1)
               WRITE (IWRITE,520) COND1
               DO 230 I=1,N
                  DO 220 J=1,N
                     AINV(I,J)=(0.0,0.0)
 220               CONTINUE
                  AINV(I,I)=(1.0,0.0)
 230           CONTINUE
               CALL CSPMSL(N,MR,MC,IA,JA,AS,IL,AINV,LDA,N)
               AMNORM=ANRM(N,ASAVE,LDA,AINV,LDA,ANORM,CONDM1)
               WRITE(IWRITE,260)CONDM1
               DO 250 I=1,N
                   DO 240 J=1,N
                      AINV(I,J)=(0.0,0.0)
 240               CONTINUE
                   AINV(I,I)=(1.0,0.0)
 250           CONTINUE
               CALL CSPFSL(N,FR,FC,IWORK,UL,AINV,LDA,N)
               AFNORM=ANRM(N,ASAVE,LDA,AINV,LDA,ANORM,CONDF1)
               WRITE(IWRITE,270)CONDF1
 260           FORMAT(19H COND FROM AINV*A M,E15.5)
 270           FORMAT(19H COND FROM AINV*A F,E15.5)
C
C              SOLVE  A*X = B
C
               CALL CGEFS(N,A,LDA,B,LDA,2,IPVT)
               CALL CGEBS(N,A,LDA,B,LDA,2)
C
C              MORE SPAR COMPARE
C
               DO 280 I=1,N
                  MC(I)=MR(I)
                  FC(I)=FR(I)
 280           CONTINUE
               CALL SPMSF(N,MR,MIC,IA,JSPSAV,IWORK,IWMAX,IFILL)
               DO 290 J=1,N
                  I=MR(J)
                  NUM=IA(I+1)-IA(I)
                 IB=IA(I)
                   CALL CSPMIN(N,MIC,IWORK,J,SPSAV(IB),JSPSAV(IB),NUM,
     1          J,UL)
 290           CONTINUE
              EPS1=EPS*ANORM
               CALL CSPMNF(N,IWORK,UL,EPS1,GROWTH)
               IF (NERROR(IMNF).EQ.0) GO TO 300
                 CALL ERROFF
                 GO TO 310
 300           CONTINUE
               CALL CSPSOL(N,MR,MC,IWORK,UL,BM,LDA,2)
 310           CONTINUE
               CALL SPFSF(N,MR,MC,IIROW,IWORK,IWMAX,IFILL)
               CALL CSPFNF(N,MR,MC,AROW,IWORK,UL,GROWTH,EPS1)
               IF (NERROR(IFNF).EQ.0) GO TO 320
                  CALL ERROFF
                  GO TO 330
 320           CONTINUE
               CALL CSPSOL(N,MR,MC,IWORK,UL,BF,LDA,2)
 330           CONTINUE
               CALL CSPMLE(N,.TRUE.,IA,JSPSAV,SPSAV,ISIZE,BM2,LDA,2)
               IF(NERROR(IMLE).NE.0)CALL ERROFF
               CALL CSPFLE(N,.TRUE.,AROW,ISIZE,BF2,LDA,2)
                IF (NERROR(IFLE).NE.0) CALL ERROFF
C COMPUTE ERRORS IN SOLUTION WITH 2 RIGHT HAND SIDES
               IF(IFLE+IMLE+IMNF+IFNF.EQ.0)GO TO 350
                 WRITE(IWRITE,340)IMNF,IFNF,IMLE,IFLE
 340    FORMAT(42H NUMERICALLY SINGULAR- IMNF,IFNF,IMLE,IFLE,4I5)
                GO TO 440
 350           CONTINUE
               DO 360 I=1,2
                  XN(I)=SCASUM(N,B(1,I),1)
                  CALL CAXPY(N,(-1.0,0.0),XEXACT(1,I),1,B(1,I),1)
                  CALL CAXPY(N,(-1.0,0.0),XEXACT(1,I),1,BM(1,I),1)
                  CALL CAXPY(N,(-1.0,0.0),XEXACT(1,I),1,BF(1,I),1)
                  CALL CAXPY(N,(-1.0,0.0),XEXACT(1,I),1,BM2(1,I),1)
                  CALL CAXPY(N,(-1.0,0.0),XEXACT(1,I),1,BF2(1,I),1)
                  EB(I)=SCASUM(N,B(1,I),1)
                  EBM(I)=SCASUM(N,BM(1,I),1)
                  EBF(I)=SCASUM(N,BF(1,I),1)
                  EBM2(I)=SCASUM(N,BM2(1,I),1)
                  EBF2(I)=SCASUM(N,BF2(1,I),1)
 360              CONTINUE
                  WRITE(IWRITE,370)
 370              FORMAT(29H ERRORS IN SOLUTION -ABSOLUTE)
                  WRITE(IWRITE,380)
 380           FORMAT(15X,42H GENERAL SPMN     SPFN       SPMLE   SPFLE)
                  DO 400 I=1,2
                  WRITE(IWRITE,390)I,EB(I),EBM(I),EBF(I),EBM2(I),EBF2(I)
 390      FORMAT(8H PROBLEM,I3,5E15.5)
 400              CONTINUE
C
C
               Q(1) = RCOND/COND1
               Q(2) = CONDM/COND1
               Q(3) = COND1/RCOND
               Q(4)=COND1/CONDM
               Q(5) = EB(1)/(EPS*COND1*XN(1))
               Q(6) = EBM(1)/(EPS*COND1*XN(1))
               Q(7) = EBM2(1)/(EPS*COND1*XN(1))
               Q(8) = AINORM/(EPS*COND1)
               Q(9)=AMNORM/(EPS*COND1)
               WRITE (IWRITE,470)
               WRITE (IWRITE,560)
               WRITE (IWRITE,470)
               WRITE (IWRITE,620)
               WRITE (IWRITE,630)
               WRITE (IWRITE,640)
               WRITE (IWRITE,470)
               WRITE (IWRITE,680) (Q(I), I = 1, 9)
               WRITE (IWRITE,470)
C
C              LOOK FOR SUSPICIOUS RATIOS
C
               QS(1) = 1.0E0 + 4.0E0*EPS
               QS(2)=QS(1)
               QS(4) = 10.0E0
               QS(3)=10.0E0
               EN = FLOAT(N)
               IF (N .EQ. 1) EN = 2.0E0
               DO 410 I = 3, 10
                  QS(I) = EN
 410           CONTINUE
               KOUNT = 0
               DO 430 I = 1, 9
                  IQ(I) = 0
                  IF (Q(I) .LE. QS(I)) GO TO 420
                     IQ(I) = 1
                     KSUSP(I) = KSUSP(I) + 1
                     KOUNT = KOUNT + 1
 420              CONTINUE
 430           CONTINUE
               IF (KOUNT .EQ. 0) WRITE (IWRITE,750)
               IF (KOUNT .NE. 0) WRITE (IWRITE,760) (IQ(I), I = 1, 9)
               WRITE (IWRITE,470)
 440        CONTINUE
 450     CONTINUE
C
         WRITE (IWRITE,570)
         KASE = KASE + 1
      GO TO 20
 460  CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (IWRITE,580)
      KASE = KASE - 1
      WRITE (IWRITE,590) KASE
      WRITE (IWRITE,600) KSING
      WRITE (IWRITE,610) KSUSP
      WRITE (IWRITE,730)
      CALL LEAVE
      RETURN
C
C     MOST FORMATS, ALSO SOME IN SGEXX
C
 470  FORMAT (1H )
 480  FORMAT (29H1  PORT  TESTER, CGE**, CSP**)
 490  FORMAT ( / 14H EPSILON     =, 1PE13.5)
 500  FORMAT ( / 19H EXACT SINGULARITY. /)
 510  FORMAT ( / 16H MAYBE SINGULAR. /)
 520  FORMAT (14H ACTUAL COND =, 1PE13.5)
 530  FORMAT(/14H X AND XBAND =)
 540  FORMAT (14H ERROR NORMS =, 2(1PE13.5))
 550  FORMAT (14H RESID NORMS =, 2(1PE13.5))
 560  FORMAT (26H TEST RATIOS.. E = EPSILON)
 570  FORMAT ( / 14H ************* /)
 580  FORMAT (8H1SUMMARY)
 590  FORMAT (18H NUMBER OF TESTS =, I4)
 600  FORMAT (30H NUMBER OF SINGULAR MATRICES =, I4)
 610  FORMAT (30H NUMBER OF SUSPICIOUS RATIOS =, 10I4)
 620  FORMAT(40H    COND  COND(M)  ACTUAL  ACTUAL  ERROR,
     1       36H ERROR(M) ERROR(M2) A*AI-I A*AI-I(M))
 630  FORMAT (9(8H   -----))
 640  FORMAT(42H    ACTUAL ACTUAL   COND  COND(M) E*COND*X,
     1       31H E*CONDX E*COND*X E*COND E*COND)
 650  FORMAT (14H NORM(A)     =, 1PE13.5)
 660  FORMAT (14H NORM(A - LU)=, 1PE13.5)
 670  FORMAT (14H NORM(A*AI-I)=, 2(1PE13.5))
 680  FORMAT (10(1X, F7.2))
 690  FORMAT (1H , 6G11.4)
 700  FORMAT (14H 1/COND      =, 1PE13.5)
 710  FORMAT (2G14.6)
 720  FORMAT (5H ML =, I2, 6H  MU =, I2)
 730  FORMAT ( / 12H END OF TEST)
 740  FORMAT(7H COND =,1PE13.5,13H COND(BAND) =,1PE13.5 /)
 750  FORMAT (21H NO SUSPICIOUS RATIOS)
 760  FORMAT (I8, 9I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
 770  FORMAT (29H THIS VERSION DATED 03/11/78.)
      END
      SUBROUTINE CGEXX(A,LDA,N,KASE,IWRITE)
      INTEGER LDA,N,KASE,IWRITE
      COMPLEX A(LDA,1)
C
C     GENERATES COMPLEX GENERAL TEST MATRICES
C
C     EXTERNAL CMACH
C     FORTRAN CMPLX,FLOAT,MAX0
      COMPLEX T1,T2
      REAL CMACH,HUGE,TINY
      INTEGER I,J,IRAND
C
      GO TO (10, 10, 10, 60, 60, 80, 80, 80, 120, 160, 200, 240, 280,
     *       320, 360, 410, 460, 520), KASE
C
C     KASE 1, 2 AND 3
C
 10   CONTINUE
         N = 3*KASE
         WRITE (IWRITE,20) KASE,N
 20      FORMAT (5H KASE, I3, 3X, 16HHILBERT SLICE    / 4H N =, I4)
         DO 50 J = 1, N
            DO 40 I = 1, N
               A(I,J) = (0.0E0,0.0E0)
               IF (I .GT. J + 2) GO TO 30
               IF (I .LT. J - 3) GO TO 30
                  A(I,J) = (1.0E0,0.0E0)/CMPLX(FLOAT(I+J-1),1.0E0)
 30            CONTINUE
 40         CONTINUE
 50      CONTINUE
      GO TO 530
C
C     KASE 4 AND 5
C
 60   CONTINUE
         N = 1
         WRITE (IWRITE,70) KASE,N
 70      FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = (3.0E0,1.0E0)
         IF (KASE .EQ. 5) A(1,1) = (0.0E0,0.0E0)
      GO TO 530
C
C     KASE 6, 7 AND 8
C
 80   CONTINUE
         N = 15
         WRITE (IWRITE,90) KASE,N
 90      FORMAT (5H KASE, I3, 3X, 16HTRIDIAGONAL      / 4H N =, I4)
         T1 = (1.0E0,0.0E0)
         T2 = (1.0E0,0.0E0)
         IF (KASE .EQ. 7) T1 = (100.0E0,100.0E0)
         IF (KASE .EQ. 8) T2 = (100.0E0,100.0E0)
         DO 110 I = 1, N
            DO 100 J = 1, N
               A(I,J) = (0.0E0,0.0E0)
               IF (I .EQ. J) A(I,I) = (4.0E0,0.0E0)
               IF (I .EQ. J - 1) A(I,J) = T1
               IF (I .EQ. J + 1) A(I,J) = T2
 100        CONTINUE
 110     CONTINUE
      GO TO 530
C
C     KASE 9
C
 120  CONTINUE
         N = 5
         WRITE (IWRITE,130) KASE,N
 130     FORMAT (5H KASE, I3, 3X, 16HRANK ONE         / 4H N =, I4)
         DO 150 I = 1, N
            DO 140 J = 1, N
               A(I,J) = CMPLX(10.0E0**(I-J),0.0E0)
 140        CONTINUE
 150     CONTINUE
      GO TO 530
C
C     KASE 10
C
 160  CONTINUE
         N = 4
         WRITE (IWRITE,170) KASE,N
 170     FORMAT (5H KASE, I3, 3X, 16HZERO COLUMN      / 4H N =, I4)
         DO 190 I = 1, N
            DO 180 J = 1, N
               T1 = CMPLX(FLOAT(J-3),0.0E0)
               T2 = CMPLX(FLOAT(I),0.0E0)
               A(I,J) = T1/T2
 180        CONTINUE
 190     CONTINUE
      GO TO 530
C
C     KASE 11
C
 200  CONTINUE
         N = 5
         WRITE (IWRITE,210) KASE,N
 210     FORMAT (5H KASE, I3, 3X, 16HTEST COND        / 4H N =, I4)
         DO 230 I = 1, N
            DO 220 J = 1, N
               IF (I .EQ. J) A(I,J) = CMPLX(FLOAT(I),0.0E0)
               IF (I .GT. J) A(I,J) = CMPLX(FLOAT(J-2),0.0E0)
               IF (I .LT. J) A(I,J) = CMPLX(FLOAT(I-2),0.0E0)
 220        CONTINUE
 230     CONTINUE
      GO TO 530
C
C     KASE 12
C
 240  CONTINUE
         N = 3
         WRITE (IWRITE,250) KASE,N
 250     FORMAT (5H KASE, I3, 3X, 16HIDENTITY         / 4H N =, I4)
         DO 270 I = 1, N
            DO 260 J = 1, N
               IF (I .EQ. J) A(I,I) = (1.0E0,0.0E0)
               IF (I .NE. J) A(I,J) = (0.0E0,0.0E0)
 260        CONTINUE
 270     CONTINUE
      GO TO 530
C
C     KASE 13
C
 280  CONTINUE
         N = 6
         WRITE (IWRITE,290) KASE,N
 290     FORMAT (5H KASE, I3, 3X, 16HUPPER TRIANGULAR / 4H N =, I4)
         DO 310 I = 1, N
            DO 300 J = 1, N
               IF (I .GT. J) A(I,J) = (0.0E0,0.0E0)
               IF (I .LE. J) A(I,J) = CMPLX(FLOAT(J-I+1),FLOAT(J-I))
 300        CONTINUE
 310     CONTINUE
      GO TO 530
C
C     KASE 14
C
 320  CONTINUE
         N = 6
         WRITE (IWRITE,330) KASE,N
 330     FORMAT (5H KASE, I3, 3X, 16HLOWER TRIANGULAR / 4H N =, I4)
         DO 350 I = 1, N
            DO 340 J = 1, N
               IF (I .LT. J) A(I,J) = (0.0E0,0.0E0)
               IF (I .GE. J) A(I,J) = CMPLX(FLOAT(I-J+1),-1.0E0)
 340        CONTINUE
 350     CONTINUE
      GO TO 530
C
C     KASE 15
C
 360  CONTINUE
         N = 5
         WRITE (IWRITE,370) KASE,N
 370     FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY =R1MACH(1)*FLOAT(N*N)*100.0
         TINY=SQRT(TINY)
         WRITE (IWRITE,380) TINY
 380     FORMAT (14H TINY        =, 1PE13.5)
         DO 400 I = 1, N
            DO 390 J = 1, N
               A(I,J) = CMPLX(TINY*FLOAT(J)/FLOAT(MAX0(I,J)),0.0E0)
 390        CONTINUE
 400     CONTINUE
      GO TO 530
C
C     KASE 16
C
 410  CONTINUE
         N = 5
         WRITE (IWRITE,420) KASE,N
 420     FORMAT (5H KASE, I3, 3X, 16HNEAR OVERFLOW    / 4H N =, I4)
         HUGE= SQRT(R1MACH(2))/FLOAT(200*N*N)
         WRITE (IWRITE,430) HUGE
 430     FORMAT (14H HUGE        =, 1PE13.5)
         DO 450 I = 1, N
            DO 440 J = 1, N
               A(I,J) = CMPLX(HUGE*FLOAT(J)/FLOAT(MAX0(I,J)),0.0E0)
 440        CONTINUE
 450     CONTINUE
      GO TO 530
C
 460  CONTINUE
C
C THIS IS RANDOM SPARSE ROUTINE
C
        WRITE(IWRITE,470)KASE,N
 470    FORMAT(5H KASE, I3,14H RANDOM SPARSE,/4H N =, I4)
        DO 490 I=1,N
           DO 480 J=1,N
              A(I,J)=(0.0E0,0.0E0)
 480       CONTINUE
 490     CONTINUE
        K=N/2
         DO 510 I=1,N
            J=0
            IRAND = UNI(0)
 500        JJ=J+K*IRAND+1
             IF (J.LT.I.AND.JJ.GT.I) JJ=I
             J=JJ
             IF (J.GT.N) GO TO 510
               A(I,J)=CMPLX(UNI(0),UNI(0))
               GO TO 500
 510     CONTINUE
         GO TO 530
C
 520       CONTINUE
         N = 0
 530  CONTINUE
      RETURN
C
      END
      SUBROUTINE AROW(I, ROW, JCOL, NUM)
      INTEGER JCOL(200), JPIB, IA(101), JA(2500)
      COMPLEX ROW(200),A(2500)
      COMMON /SPA/IA,JA,A
       NUM=IA(I+1)-IA(I)
      IB=IA(I)-1
      DO 10 J=1,NUM
          JPIB = J+IB
          ROW(J)=A(JPIB)
          JCOL(J)=JA(JPIB)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE IIROW(I,  JCOL, NUM)
      INTEGER JCOL(200),JPIB,IA(101),JA(2500)
      COMPLEX A(2500)
      COMMON /SPA/IA,JA,A
      NUM=IA(I+1)-IA(I)
      IB=IA(I)-1
      DO 10 J=1,NUM
          JPIB = J+IB
          JCOL(J)=JA(JPIB)
 10   CONTINUE
      RETURN
      END
      REAL FUNCTION ANRM(N,AS,LDA1,AINV,LDA,ANORM,COND1)
      COMPLEX AS(LDA,N),AINV(LDA1,N)
C
C     THIS SUBROUTINE COMPUTES THE CONDITION NUMBER FROM
C     A AND AINVERSE
C
C     IT ALSO COMPUTES THE ERROR A*AINV-I AND STORES IT IN ANRM
      COMPLEX B(200),T
      ANRM=0.0E0
      DO 30 J =1,N
         DO 10 I=1,N
            B(I)=(0.0,0.0)
 10      CONTINUE
         DO 20 K=1,N
            T=AINV(K,J)
            CALL CAXPY(N,T,AS(1,K),1,B,1)
 20       CONTINUE
         B(J)=B(J)-(1.0,0.0)
         ANRM=AMAX1(ANRM,SCASUM(N,B,1))
 30   CONTINUE
      AIN=0.0E0
      DO 40 I=1,N
         AIN=AMAX1(AIN,SCASUM(N,AINV(1,I),1))
 40   CONTINUE
      COND1=ANORM*AIN
      RETURN
      END
         SUBROUTINE CGETSP(N,A,LDA,IA,JA,AS,L)
         COMPLEX A(LDA,N),AS(1)
         INTEGER JA(1),IA(1)
         L=0
         IA(1)=1
         DO 30 I=1,N
            DO 10 J=1,N
               IF (CABS1(A(I,J)).EQ.0.0E0) GO TO 10
               L=L+1
               AS(L)=A(I,J)
               JA(L)=J
 10         CONTINUE
           IF (L.GE.IA(I))GO TO 20
               L=L+1
               AS(L)=(0.0,0.0)
               JA(L)=1
 20        CONTINUE
            IA(I+1)=L+1
 30     CONTINUE
        RETURN
        END
