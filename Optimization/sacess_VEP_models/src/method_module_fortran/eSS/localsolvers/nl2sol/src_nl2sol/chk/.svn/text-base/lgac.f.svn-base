C$TEST LGAC
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LGAC
C***********************************************************************
C
C  TEST OF THE PORT PROGRAMS CGECE AND FRIENDS
C
C***********************************************************************
C     MAIN PROGRAM
      INTEGER IWRITE,I1MACH
C     ALLOW 5000 UNDERFLOWS.
C
C     SET OUTPUT UNIT NUMBER
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
C        CGECE,CGEFS,CGEBS,CGEML,CBACE,CBABS,CBAFS,CBALE,CBAML
C
C
C     SUBROUTINES AND FUNCTIONS
C
C     PORT CGECE,CGEFS,CGEBS,CGEML,CBACE,CBAFS,CBABS,CBAML,CBALE
C     PORT UTILITIES ERROFF,ENTER,LEAVE,R1MACH
C     EXTERNAL SGEXX
C     BLAS CAXPY,CSCAL,SCASUM,CABS1
C     FORTRAN AMAX1,FLOAT,MAX0,MIN0,CABS
C
C     INTERNAL VARIABLES
C
      INTEGER I,IPVT(15),IPVTB(15),IQ(8),I1,I2,J
      INTEGER K,KASE,KB,KBFAIL,KOUNT,KP1,KSING,KSUSP(8)
      INTEGER L,LDA,LDAB,IWRITE,M,ML,MU,N,NM1,NPRINT
      REAL CMACH,HUGE,TINY
      REAL AINORM,ANORM,SMACH,COND,COND1,EN,ENORM,EPS
      REAL ETNORM,FNI,FNORM,ONEPX,RCOND,RCONDB,RNORM
      REAL RTNORM,Q(10),QS(10),SCASUM,XNORM,DENOM
      COMPLEX A(15,15),AB(43,15),AINV(15,15),ASAVE(15,15),AL(43,15)
      COMPLEX B(15),BT(15),CDOTU,DET(2),DETB(2)
      COMPLEX X(15),XB(15),XEXACT(15),XT(15),XTB(15),T
      COMPLEX ABSAVE(43,15),AINVB(15,15)
      LOGICAL KBF
C
      LDA = 15
      CALL ENTER(1)
      LDAB = 43
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT =3
C
      WRITE (IWRITE,380)
      WRITE (IWRITE,670)
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
      WRITE (IWRITE,390) EPS
      WRITE (IWRITE,370)
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
         IF(N. LE.0) GO TO 360
         ANORM = 0.0E0
         DO 30 J = 1, N
            ANORM = AMAX1(ANORM,SCASUM(N,A(1,J),1))
 30      CONTINUE
         WRITE (IWRITE,550) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (IWRITE,370)
            DO 40 I = 1, N
               WRITE (IWRITE,590) (A(I,J), J = 1, N)
 40         CONTINUE
            WRITE (IWRITE,370)
 50      CONTINUE
C
C        GENERATE EXACT SOLUTION
C
         XEXACT(1) = (1.0E0,1.0E0)
         IF (N .GE. 2) XEXACT(2) = (0.0E0,0.0E0)
         IF (N .LE. 2) GO TO 70
            DO 60 I = 3, N
               XEXACT(I) = -XEXACT(I-2)
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
         CALL MOVEFC(N,B,X)
C
C        FACTOR AND ESTIMATE CONDITION
C
         CALL CGECE(N,A,LDA,IPVT,RCOND)
         IF (NERROR(INE).NE.0) CALL  ERROFF
C
C
C        FACTOR BAND FORM AND COMPARE
C
         KBF = .FALSE.
         ML = 0
         MU = 0
         DO 120 J = 1, N
            DO 110 I = 1, N
               IF (CABS1(ASAVE(I,J)) .EQ. 0.0E0) GO TO 100
                  IF (I .LT. J) MU = MAX0(MU,J-I)
                  IF (I .GT. J) ML = MAX0(ML,I-J)
 100           CONTINUE
 110        CONTINUE
 120     CONTINUE
         MLP1=ML+1
         WRITE (IWRITE,620) ML,MU
            M = ML + MU + 1
            DO 140 J = 1, N
               I1 = MAX0(1,J-MU)
               I2 = MIN0(N,J+ML)
               DO 130 I = I1, I2
                  K=ML+1+J-I
                  AB(K,I) = ASAVE(I,J)
                  ABSAVE(K,I)=ASAVE(I,J)
 130           CONTINUE
 140        CONTINUE
C
            CALL CBACE(N,MLP1,M,AB,LDAB,AL,LDAB,IPVTB,MU,RCONDB)
          WRITE(IWRITE,640)RCOND,RCONDB
C
C
C           TEST FOR SINGULARITY
C
           IF (INE+NERROR(IRE).EQ.0) GO TO 160
             IF (IRE.NE.0) CALL ERROFF
             IF (IRE.NE.INE) WRITE(IWRITE,150)
 150         FORMAT(35H BAND AND GENERAL ROUTINES DISAGREE)
               WRITE (IWRITE,400)
               KSING = KSING + 1
            GO TO 340
 160        CONTINUE
C
C              COMPUTE INVERSE AND COND1 = TRUE CONDITION
C
               DO 180 J = 1, N
                  DO 170 I = 1, N
                     AINV(I,J) = (0.E0,0.E0)
                     AINVB(I,J)=(0.E0,0.E0)
 170              CONTINUE
               AINV(J,J)=(1.E0,0.E0)
               AINVB(J,J)=(1.E0,0.E0)
 180           CONTINUE
               CALL CGEFS(N,A,LDA,AINV,LDA,N,IPVT)
               CALL CGEBS(N,A,LDA,AINV,LDA,N)
               CALL CBAFS(N,MLP1,AL,LDAB,IPVTB,AINVB,LDA,N)
               CALL CBABS(N,AB,LDAB,AINVB,LDA,N,MU)
               AINORM = 0.0E0
               DO 190 J = 1, N
                  AINORM = AMAX1(AINORM,SCASUM(N,AINV(1,J),1))
 190           CONTINUE
               COND1 = ANORM*AINORM
               WRITE (IWRITE,420) COND1
C
C              SOLVE  A*X = B
C
               CALL CGEFS(N,A,LDA,X,N,1,IPVT)
               CALL CGEBS(N,A,LDA,X,N,1)
C
C              MORE BAND COMPARE
C
C              TEST CONSISTENCY OF BAML AND BALE
C
               CALL CBAML(N,MLP1,M,ABSAVE,LDAB,XEXACT,XB)
               CALL CBALE(N,MLP1,M,ABSAVE,LDAB,XB,N,1)
               IF (NERROR(IRE).EQ.0) GO TO 200
                    CALL ERROFF
                    WRITE(IWRITE,410)
                    GO TO 340
 200            CONTINUE
               IF (N .GT. NPRINT) GO TO 220
                  WRITE (IWRITE,430)
                  DO 210 I = 1, N
                     WRITE (IWRITE,610) X(I),XB(I)
 210              CONTINUE
                  WRITE (IWRITE,370)
 220           CONTINUE
C
C              RECONSTRUCT  A  FROM TRIANGULAR FACTORS , L AND U
C
               NM1 = N - 1
               IF (NM1 .LT. 1) GO TO 250
               DO 240 KB = 1, NM1
                  K = N - KB
                  NMK=N-K
                  KP1 = K + 1
                  L = IPVT(K)
                  DO 230 J = KP1, N
                     T = -A(K,J)
                     CALL CAXPY(NMK,T,A(K+1,K),1,A(K+1,J),1)
                     T = A(L,J)
                     A(L,J) = A(K,J)
                     A(K,J) = T
 230              CONTINUE
                  T = -A(K,K)
                  CALL CSCAL(NMK,T,A(K+1,K),1)
                  T = A(L,K)
                  A(L,K) = A(K,K)
                  A(K,K) = T
 240           CONTINUE
 250           CONTINUE
C
C              COMPUTE ERRORS AND RESIDUALS
C                 E  =  X - XEXACT
C                 EB =  XB - XEXACT
C                 R  =  B - A*X
C                 F  =  A - L*U
C                 AI =  A*INV(A) - I
C                 AIB = A(BAND)*INV(A(BAND)) - I
C
               XNORM = SCASUM(N,X,1)
               ENORM = 0.0E0
               EBNORM=0.E0
               FNORM = 0.0E0
               DO 270 J = 1, N
                  ENORM = ENORM + CABS(X(J)-XEXACT(J))
                  EBNORM = EBNORM + CABS(XB(J) - XEXACT(J))
                  T = -X(J)
                  CALL CAXPY(N,T,ASAVE(1,J),1,B,1)
                  FNI = 0.0E0
                  DO 260 I = 1, N
                     FNI = FNI + CABS(ASAVE(I,J)-A(I,J))
 260              CONTINUE
                  IF (FNI .GT. FNORM) FNORM = FNI
 270           CONTINUE
               RNORM = SCASUM(N,B,1)
C
C              A*INV(A) - I
C
               AINORM = 0.0E0
               AIBNO=0.0E0
               DO 300 J = 1, N
                  DO 280 I = 1, N
                     B(I) = (0.E0,0.E0)
                     XB(I) = (0.E0,0.E0)
 280              CONTINUE
                  DO 290 K = 1, N
                     T = AINV(K,J)
                     CALL CAXPY(N,T,ASAVE(1,K),1,B,1)
                     T=AINVB(K,J)
                     CALL CAXPY(N,T,ASAVE(1,K),1,XB,1)
 290              CONTINUE
                  B(J) = B(J) - (1.0E0,0.0E0)
                  XB(J) = XB(J) -(1.0E0,0.0E0)
                  AIBNO=AMAX1(AIBNO,SCASUM(N,XB,1))
                  AINORM = AMAX1(AINORM,SCASUM(N,B,1))
 300           CONTINUE
C
               WRITE (IWRITE,440) ENORM,EBNORM
               WRITE (IWRITE,450) RNORM
               WRITE (IWRITE,560) FNORM
               WRITE (IWRITE,570) AINORM,AIBNO
C
C              COMPUTE TEST RATIOS
C
               Q(1) = RCOND/COND1
               Q(2) = RCONDB/COND1
               Q(3) = COND1/RCOND
               Q(4)=COND1/RCONDB
               Q(5) = ENORM/(EPS*RCOND*XNORM)
               Q(6) = EBNORM/(EPS*RCOND*XNORM)
              DENOM=AMAX1(100.*R1MACH(1),EPS*ANORM*XNORM)
               Q(7)=RNORM/DENOM
               DENOM=AMAX1(100.*R1MACH(1),EPS*ANORM)
               Q(8)=FNORM/DENOM
               Q(9) = AINORM/(EPS*RCOND)
               Q(10)=AIBNO/(EPS*RCOND)
               WRITE (IWRITE,370)
               WRITE (IWRITE,460)
               WRITE (IWRITE,370)
               WRITE (IWRITE,520)
               WRITE (IWRITE,530)
               WRITE (IWRITE,540)
               WRITE (IWRITE,370)
               WRITE (IWRITE,580) (Q(I), I = 1, 10)
               WRITE (IWRITE,370)
C
C              LOOK FOR SUSPICIOUS RATIOS
C
               QS(1) = 1.0E0 + 4.0E0*EPS
               QS(2)=QS(1)
               QS(4) = 10.0E0
               QS(3)=10.0E0
               EN = FLOAT(N)
               IF (N.LT.3)EN=3.0
               DO 310 I = 3, 10
                  QS(I) = EN
 310           CONTINUE
               KOUNT = 0
            IF (KASE.EQ.15)QS(7)=QS(7)*100.0
            IF (KASE.EQ.15)QS(8)=QS(8)*100.0
               DO 330 I = 1, 10
                  IQ(I) = 0
                  IF (Q(I) .LE. QS(I)) GO TO 320
                     IQ(I) = 1
                     KSUSP(I) = KSUSP(I) + 1
                     KOUNT = KOUNT + 1
 320              CONTINUE
 330           CONTINUE
               IF (KOUNT .EQ. 0) WRITE (IWRITE,650)
               IF (KOUNT .NE. 0) WRITE (IWRITE,660) (IQ(I), I = 1, 10)
               WRITE (IWRITE,370)
 340        CONTINUE
 350     CONTINUE
C
         WRITE (IWRITE,470)
         KASE = KASE + 1
      GO TO 20
 360  CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (IWRITE,480)
      KASE = KASE - 1
      WRITE (IWRITE,490) KASE
      WRITE (IWRITE,500) KSING
      WRITE (IWRITE,510) KSUSP
      WRITE (IWRITE,630)
      CALL LEAVE
      RETURN
C
C     MOST FORMATS, ALSO SOME IN SGEXX
C
 370  FORMAT (1H )
 380  FORMAT (29H1  PORT  TESTER, CGE**, CBA**)
 390  FORMAT ( / 14H EPSILON     =, 1PE17.5)
 400  FORMAT ( / 19H EXACT SINGULARITY. /)
 410  FORMAT ( / 16H MAYBE SINGULAR. /)
 420  FORMAT (14H ACTUAL COND =, 1PE17.5)
 430  FORMAT(/14H X AND XBAND =)
 440  FORMAT (14H ERROR NORMS =, 2(1PE17.5))
 450  FORMAT (14H RESID NORMS =, 2(1PE17.5))
 460  FORMAT (26H TEST RATIOS.. E = EPSILON)
 470  FORMAT ( / 14H ************* /)
 480  FORMAT (8H1SUMMARY)
 490  FORMAT (18H NUMBER OF TESTS =, I4)
 500  FORMAT (30H NUMBER OF SINGULAR MATRICES =, I4)
 510  FORMAT (30H NUMBER OF SUSPICIOUS RATIOS =, 10I4)
 520  FORMAT(42H    COND  COND(B)  ACTUAL  ACTUAL  ERROR  ,
     1       40HERROR(B)  RESID  A-LU  A*AI-I  A*AI-I(B))
 530  FORMAT (10(8H   -----))
 540  FORMAT(42H    ACTUAL ACTUAL   COND  COND(B) E*COND*X,
     1       40H E*COND*X E*A*X    E*A   E*COND  E*COND )
 550  FORMAT (14H NORM(A)     =, 1PE17.5)
 560  FORMAT (14H NORM(A - LU)=, 1PE17.5)
 570  FORMAT (14H NORM(A*AI-I)=, 2(1PE17.5))
 580  FORMAT (10(1X, F7.2))
 590  FORMAT (1H , 6E15.5)
 600  FORMAT (14H 1/COND      =, 1PE17.5)
 610  FORMAT (4G14.6)
 620  FORMAT (5H ML =, I2, 6H  MU =, I2)
 630  FORMAT ( / 12H END OF TEST)
 640  FORMAT(7H COND =,1PE17.5,13H COND(BAND) =,1PE17.5 /)
 650  FORMAT (21H NO SUSPICIOUS RATIOS)
 660  FORMAT (I8, 9I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
 670  FORMAT (29H THIS VERSION DATED 03/11/78.)
      END
      SUBROUTINE CGEXX(A,LDA,N,KASE,IWRITE)
C
C     GENERATES COMPLEX GENERAL TEST MATRICES
C
C     EXTERNAL CMACH
C     FORTRAN CMPLX,FLOAT,MAX0
      INTEGER LDA,N,KASE,IWRITE
      COMPLEX A(LDA,1)
      COMPLEX T1,T2
      GO TO (10, 10, 10, 60, 60, 80, 80, 80, 120, 160, 200, 240, 280,
     *       320, 360, 410, 460), KASE
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
      GO TO 470
C
C     KASE 4 AND 5
C
 60   CONTINUE
         N = 1
         WRITE (IWRITE,70) KASE,N
 70      FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = (3.0E0,1.0E0)
         IF (KASE .EQ. 5) A(1,1) = (0.0E0,0.0E0)
      GO TO 470
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
      GO TO 470
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
      GO TO 470
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
      GO TO 470
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
      GO TO 470
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
      GO TO 470
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
      GO TO 470
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
      GO TO 470
C
C     KASE 15
C
 360  CONTINUE
         N = 5
         WRITE (IWRITE,370) KASE,N
 370     FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY =R1MACH(1)*FLOAT(N*N)*100.0
         WRITE (IWRITE,380) TINY
 380     FORMAT (14H TINY        =, 1PE13.5)
         DO 400 I = 1, N
            DO 390 J = 1, N
               A(I,J) = CMPLX(TINY*FLOAT(J)/FLOAT(MAX0(I,J)),0.0E0)
 390        CONTINUE
 400     CONTINUE
      GO TO 470
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
      GO TO 470
C
 460  CONTINUE
         N = 0
 470  CONTINUE
      RETURN
C
      END
