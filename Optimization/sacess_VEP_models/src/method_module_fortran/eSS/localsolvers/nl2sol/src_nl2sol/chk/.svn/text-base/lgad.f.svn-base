C$TEST LGAD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LGAD
C***********************************************************************
C
C  TEST OF THE PORT PROGRAMS DGECE AND FRIENDS
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
C        DGECE,DGEFS,DGEBS,DGEML,,DBACE,,DBABS,,DBAFS,,DBALE,,DBAML
C
C
C     SUBROUTINES AND FUNCTIONS
C
C     PORT DGECE,DGEFS,DGEBS,DGEML,DBACE,DBAFS,DBABS,DBAML,DBALE
C     PORT UTILITIES ERROFF,ENTER,LEAVE,D1MACH
C     EXTERNAL SGEXX
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN ABS,DMAX1,FLOAT,MAX0,MIN0
C
C     INTERNAL VARIABLES
C
      INTEGER I,IPVT(15),IPVTB(15),IQ(8),I1,I2,J
      INTEGER K,KASE,KB,KBFAIL,KOUNT,KP1,KSING,KSUSP(8)
      INTEGER L,LDA,LDAB,IWRITE,M,ML,MU,N,NM1,NPRINT
      REAL    Q(10),QS(10)
      DOUBLE PRECISION A(15,15),AB(43,15),AINV(15,15),ASAVE(15,15)
      DOUBLE PRECISION D1MACH,AIBNO,EBNORM
      DOUBLE PRECISION B(15),BT(15),DDOT,AL(43,15)
      DOUBLE PRECISION X(15),XB(15),XEXACT(15),XT(15),XTB(15),T
      DOUBLE PRECISION AINORM,ANORM,SMACH,COND,COND1,EN,ENORM,EPS
      DOUBLE PRECISION ETNORM,FNI,FNORM,ONEPX,RCOND,RCONDB,RNORM
      DOUBLE PRECISION RTNORM,DASUM,XNORM
      DOUBLE PRECISION ABSAVE(43,15),AINVB(15,15),DENOM
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
      WRITE (IWRITE,460)
      WRITE (IWRITE,880)
C
      DO 10 I = 1, 8
         KSUSP(I) = 0
   10 CONTINUE
      KSING = 0
      KBFAIL = 0
C
C     SET EPS TO ROUNDING UNIT
C
      EPS =D1MACH(4)
      WRITE (IWRITE,470) EPS
      WRITE (IWRITE,450)
C
C     START MAIN LOOP
C
      KASE = 1
   20 CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL SGEXX(A,LDA,N,KASE,IWRITE)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 440
         ANORM = 0.0D0
         DO 30 J = 1, N
            ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   30    CONTINUE
         WRITE (IWRITE,650) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (IWRITE,450)
            DO 40 I = 1, N
               WRITE (IWRITE,700) (A(I,J), J = 1, N)
   40       CONTINUE
            WRITE (IWRITE,450)
   50    CONTINUE
C
C        GENERATE EXACT SOLUTION
C
         XEXACT(1) = 1.0D0
         IF (N .GE. 2) XEXACT(2) = 0.0D0
         IF (N .LE. 2) GO TO 70
            DO 60 I = 3, N
               XEXACT(I) = -XEXACT(I-2)
   60       CONTINUE
   70    CONTINUE
C
C        SAVE MATRIX AND GENERATE R.H.S.
C
         DO 90 I = 1, N
            DO 80 J = 1, N
               ASAVE(I,J) = A(I,J)
   80       CONTINUE
   90    CONTINUE
         CALL DGEML(N,A,LDA,XEXACT,B)
         CALL MOVEFD(N,B,X)
C
C        FACTOR AND ESTIMATE CONDITION
C
         CALL DGECE(N,A,LDA,IPVT,RCOND)
         IF (NERROR(INE).NE.0) CALL  ERROFF
C
C
C        FACTOR BAND FORM AND COMPARE
C
         KBF = .FALSE.
         ML = 0
         MU = 0
         DO 140 J = 1, N
            DO 130 I = 1, N
               IF (ASAVE(I,J) .EQ. 0.0D0) GO TO 120
                  IF (I .LT. J) MU = MAX0(MU,J-I)
                  IF (I .GT. J) ML = MAX0(ML,I-J)
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
         WRITE (IWRITE,790) ML,MU
            M = ML + MU + 1
            DO 170 J = 1, N
               I1 = MAX0(1,J-MU)
               I2 = MIN0(N,J+ML)
               DO 160 I = I1, I2
                  K=ML+1+J-I
                  AB(K,I) = ASAVE(I,J)
                  ABSAVE(K,I)=ASAVE(I,J)
  160          CONTINUE
  170       CONTINUE
C
            CALL DBACE(N,ML+1,M,AB,LDAB,AL,LDAB,IPVTB,MU,RCONDB)
          WRITE(IWRITE,820)RCOND,RCONDB
C
C
C           TEST FOR SINGULARITY
C
           IF (INE+NERROR(IRE).EQ.0) GO TO 210
             IF (IRE.NE.0) CALL ERROFF
             IF (IRE.NE.INE) WRITE(IWRITE,118)
 118         FORMAT(35H BAND AND GENERAL ROUTINES DISAGREE)
               WRITE (IWRITE,480)
               KSING = KSING + 1
            GO TO 420
  210       CONTINUE
C
C              COMPUTE INVERSE AND COND1 = TRUE CONDITION
C
               DO 230 J = 1, N
                  DO 220 I = 1, N
                     AINV(I,J) = 0.D0
                     AINVB(I,J)=0.D0
  220             CONTINUE
               AINV(J,J)=1.D0
               AINVB(J,J)=1.D0
  230          CONTINUE
               CALL DGEFS(N,A,LDA,AINV,LDA,N,IPVT)
               CALL DGEBS(N,A,LDA,AINV,LDA,N)
               CALL DBAFS(N,ML+1,AL,LDAB,IPVTB,AINVB,LDA,N)
               CALL DBABS(N,AB,LDAB,AINVB,LDA,N,MU)
               AINORM = 0.0D0
               DO 240 J = 1, N
                  AINORM = DMAX1(AINORM,DASUM(N,AINV(1,J),1))
  240          CONTINUE
               COND1 = ANORM*AINORM
               WRITE (IWRITE,510) COND1
C
C              SOLVE  A*X = B
C
               CALL DGEFS(N,A,LDA,X,N,1,IPVT)
               CALL DGEBS(N,A,LDA,X,N,1)
C
C              MORE BAND COMPARE
C
C              TEST CONSISTENCY OF DBAML AND DBALE
C
               CALL DBAML(N,ML+1,M,ABSAVE,LDAB,XEXACT,XB)
               CALL DBALE(N,ML+1,M,ABSAVE,LDAB,XB,N,1)
               IF (NERROR(IRE).EQ.0) GO TO 245
                 WRITE(IWRITE,490)
                 CALL ERROFF
                 GO TO 420
 245           CONTINUE
               IF (N .GT. NPRINT) GO TO 270
                  WRITE (IWRITE,520)
                  DO 250 I = 1, N
                     WRITE (IWRITE,740) X(I),XB(I)
  250             CONTINUE
                  WRITE (IWRITE,450)
  270          CONTINUE
C
C              RECONSTRUCT  A  FROM TRIANGULAR FACTORS , L AND U
C
               NM1 = N - 1
               IF (NM1 .LT. 1) GO TO 330
               DO 320 KB = 1, NM1
                  K = N - KB
                  KP1 = K + 1
                  L = IPVT(K)
                  DO 310 J = KP1, N
                     T = -A(K,J)
                     CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
                     T = A(L,J)
                     A(L,J) = A(K,J)
                     A(K,J) = T
  310             CONTINUE
                  T = -A(K,K)
                  CALL DSCAL(N-K,T,A(K+1,K),1)
                  T = A(L,K)
                  A(L,K) = A(K,K)
                  A(K,K) = T
  320          CONTINUE
  330          CONTINUE
C
C              COMPUTE ERRORS AND RESIDUALS
C                 E  =  X - XEXACT
C                 EB =  XB - XEXACT
C                 R  =  B - A*X
C                 F  =  A - L*U
C                 AI =  A*INV(A) - I
C                 AIB = A(BAND)*INV(A(BAND)) - I
C
               XNORM = DASUM(N,X,1)
               ENORM = 0.0D0
               EBNORM=0.D0
               FNORM = 0.0D0
               DO 350 J = 1, N
                  ENORM = ENORM + DABS(X(J)-XEXACT(J))
                  EBNORM = EBNORM + DABS(XB(J) - XEXACT(J))
                  T = -X(J)
                  CALL DAXPY(N,T,ASAVE(1,J),1,B,1)
                  FNI = 0.0D0
                  DO 340 I = 1, N
                     FNI = FNI + DABS(ASAVE(I,J)-A(I,J))
  340             CONTINUE
                  IF (FNI .GT. FNORM) FNORM = FNI
  350          CONTINUE
               RNORM = DASUM(N,B,1)
C
C              A*INV(A) - I
C
               AINORM = 0.0D0
               AIBNO=0.0D0
               DO 380 J = 1, N
                  DO 360 I = 1, N
                     B(I) = 0.0D0
                     XB(I) = 0.D0
  360             CONTINUE
                  DO 370 K = 1, N
                     T = AINV(K,J)
                     CALL DAXPY(N,T,ASAVE(1,K),1,B,1)
                     T=AINVB(K,J)
                     CALL DAXPY(N,T,ASAVE(1,K),1,XB,1)
  370             CONTINUE
                  B(J) = B(J) - 1.0D0
                  XB(J) = XB(J) -1.0D0
                  AIBNO=DMAX1(AIBNO,DASUM(N,XB,1))
                  AINORM = DMAX1(AINORM,DASUM(N,B,1))
  380          CONTINUE
C
               WRITE (IWRITE,540) ENORM,EBNORM
               WRITE (IWRITE,550) RNORM
               WRITE (IWRITE,660) FNORM
               WRITE (IWRITE,670) AINORM,AIBNO
C
C              COMPUTE TEST RATIOS
C
               Q(1) = RCOND/COND1
               Q(2) = RCONDB/COND1
               Q(3) = COND1/RCOND
               Q(4)=COND1/RCONDB
               Q(5) = ENORM/(EPS*RCOND*XNORM)
               Q(6) = EBNORM/(EPS*RCOND*XNORM)
               DENOM=DMAX1(1.0D2*D1MACH(1),EPS*ANORM*XNORM)
               Q(7) = RNORM/DENOM
               DENOM=DMAX1(1.0D2*D1MACH(1),EPS*ANORM)
               Q(8) = FNORM/DENOM
               Q(9) = AINORM/(EPS*RCOND)
               Q(10)=AIBNO/(EPS*RCOND)
               WRITE (IWRITE,450)
               WRITE (IWRITE,560)
               WRITE (IWRITE,450)
               WRITE (IWRITE,620)
               WRITE (IWRITE,630)
               WRITE (IWRITE,640)
               WRITE (IWRITE,450)
               WRITE (IWRITE,690) (Q(I), I = 1, 10)
               WRITE (IWRITE,450)
C
C              LOOK FOR SUSPICIOUS RATIOS
C
               QS(1) = 1.0D0 + 4.0D0*EPS
               QS(2)=QS(1)
               QS(4) = 10.0D0
               QS(3)=10.0D0
               EN = DBLE(FLOAT(N))
               IF (N .EQ. 1) EN = 2.0D0
               DO 390 I = 3, 10
                  QS(I) = EN
  390          CONTINUE
               KOUNT = 0
               DO 410 I = 1, 10
                  IQ(I) = 0
                  IF (Q(I) .LE. QS(I)) GO TO 400
                     IQ(I) = 1
                     KSUSP(I) = KSUSP(I) + 1
                     KOUNT = KOUNT + 1
  400             CONTINUE
  410          CONTINUE
               IF (KOUNT .EQ. 0) WRITE (IWRITE,860)
               IF (KOUNT .NE. 0) WRITE (IWRITE,870) (IQ(I), I = 1, 10)
               WRITE (IWRITE,450)
  420       CONTINUE
  430    CONTINUE
C
         WRITE (IWRITE,570)
         KASE = KASE + 1
      GO TO 20
  440 CONTINUE
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
      WRITE (IWRITE,810)
      CALL LEAVE
      RETURN
C
C     MOST FORMATS, ALSO SOME IN SGEXX
C
  450 FORMAT (1H )
  460 FORMAT (29H1  PORT  TESTER, DGE**, DBA**)
  470 FORMAT ( / 14H EPSILON     =, 1PD15.5)
  480 FORMAT ( / 19H EXACT SINGULARITY. /)
  490 FORMAT ( / 16H MAYBE SINGULAR. /)
  510 FORMAT (14H ACTUAL COND =, 1PD15.5)
  520 FORMAT(/14H X AND XBAND =)
  540 FORMAT (14H ERROR NORMS =, 2(1PD15.5))
  550 FORMAT (14H RESID NORMS =, 2(1PD15.5))
  560 FORMAT (26H TEST RATIOS.. E = EPSILON)
  570 FORMAT ( / 14H ************* /)
  580 FORMAT (8H1SUMMARY)
  590 FORMAT (18H NUMBER OF TESTS =, I4)
  600 FORMAT (30H NUMBER OF SINGULAR MATRICES =, I4)
  610 FORMAT (30H NUMBER OF SUSPICIOUS RATIOS =, 10I4)
  620 FORMAT(42H    COND  COND(B)  ACTUAL  ACTUAL  ERROR  ,
     1       40HERROR(B)  RESID  A-LU  A*AI-I  A*AI-I(B))
  630 FORMAT (10(8H   -----))
  640 FORMAT(42H    ACTUAL ACTUAL   COND  COND(B) E*COND*X,
     1       40H E*COND*X E*A*X    E*A   E*COND  E*COND )
  650 FORMAT (14H NORM(A)     =, 1PD15.5)
  660 FORMAT (14H NORM(A - LU)=, 1PD15.5)
  670 FORMAT (14H NORM(A*AI-I)=, 2(1PD15.5))
  690 FORMAT (10(1X, F7.2))
  700 FORMAT (1H , 6D15.4)
  710 FORMAT (14H 1/COND      =, 1PD15.5)
  740 FORMAT (2D18.6)
  790 FORMAT (5H ML =, I2, 6H  MU =, I2)
  810 FORMAT ( / 12H END OF TEST)
  820 FORMAT(7H COND =,1PD15.5,13H COND(BAND) =,1PD15.5 /)
  860 FORMAT (21H NO SUSPICIOUS RATIOS)
  870 FORMAT (I8, 9I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
  880 FORMAT (29H THIS VERSION DATED 03/11/78.)
      END
      SUBROUTINE SGEXX(A,LDA,N,KASE,IWRITE)
C
C     GENERATES DOUBLE PRECISION GENERAL TEST MATRICES
C
C     EXTERNAL SMACH
C     FORTRAN FLOAT,MAX0
      INTEGER LDA,N,KASE,IWRITE
      INTEGER I,J
      DOUBLE PRECISION A(LDA,1)
      DOUBLE PRECISION T1,T2
      DOUBLE PRECISION D1MACH,HUGE,TINY
C
      GO TO (10, 10, 10, 60, 60, 80, 80, 80, 120, 160, 200, 240, 280,
     *       320, 360, 410, 460), KASE
C
C     KASE 1, 2 AND 3
C
   10 CONTINUE
         N = 3*KASE
         WRITE (IWRITE,20) KASE,N
   20    FORMAT (5H KASE, I3, 3X, 16HHILBERT SLICE    / 4H N =, I4)
         DO 50 J = 1, N
            DO 40 I = 1, N
               A(I,J) = 0.0D0
               IF (I .GT. J + 2) GO TO 30
               IF (I .LT. J - 3) GO TO 30
                  A(I,J) = 1.0D0/DBLE(FLOAT(I+J-1))
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      GO TO 470
C
C     KASE 4 AND 5
C
   60 CONTINUE
         N = 1
         WRITE (IWRITE,70) KASE,N
   70    FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = 3.0D0
         IF (KASE .EQ. 5) A(1,1) = 0.0D0
      GO TO 470
C
C     KASE 6, 7 AND 8
C
   80 CONTINUE
         N = 15
         WRITE (IWRITE,90) KASE,N
   90    FORMAT (5H KASE, I3, 3X, 16HTRIDIAGONAL      / 4H N =, I4)
         T1 = 1.0D0
         T2 = 1.0D0
         IF (KASE .EQ. 7) T1 = 100.0D0
         IF (KASE .EQ. 8) T2 = 100.0D0
         DO 110 I = 1, N
            DO 100 J = 1, N
               A(I,J) = 0.0D0
               IF (I .EQ. J) A(I,I) = 4.0D0
               IF (I .EQ. J - 1) A(I,J) = T1
               IF (I .EQ. J + 1) A(I,J) = T2
  100       CONTINUE
  110    CONTINUE
      GO TO 470
C
C     KASE 9
C
  120 CONTINUE
         N = 5
         WRITE (IWRITE,130) KASE,N
  130    FORMAT (5H KASE, I3, 3X, 16HRANK ONE         / 4H N =, I4)
         DO 150 I = 1, N
            DO 140 J = 1, N
               A(I,J) = 10.0D0**(I - J)
  140       CONTINUE
  150    CONTINUE
      GO TO 470
C
C     KASE 10
C
  160 CONTINUE
         N = 4
         WRITE (IWRITE,170) KASE,N
  170    FORMAT (5H KASE, I3, 3X, 16HZERO COLUMN      / 4H N =, I4)
         DO 190 I = 1, N
            DO 180 J = 1, N
               T1 = DBLE(FLOAT(J-3))
               T2 = DBLE(FLOAT(I))
               A(I,J) = T1/T2
  180       CONTINUE
  190    CONTINUE
      GO TO 470
C
C     KASE 11
C
  200 CONTINUE
         N = 5
         WRITE (IWRITE,210) KASE,N
  210    FORMAT (5H KASE, I3, 3X, 16HTEST COND        / 4H N =, I4)
         DO 230 I = 1, N
            DO 220 J = 1, N
               IF (I .EQ. J) A(I,J) = DBLE(FLOAT(I))
               IF (I .GT. J) A(I,J) = DBLE(FLOAT(J-2))
               IF (I .LT. J) A(I,J) = DBLE(FLOAT(I-2))
  220       CONTINUE
  230    CONTINUE
      GO TO 470
C
C     KASE 12
C
  240 CONTINUE
         N = 3
         WRITE (IWRITE,250) KASE,N
  250    FORMAT (5H KASE, I3, 3X, 16HIDENTITY         / 4H N =, I4)
         DO 270 I = 1, N
            DO 260 J = 1, N
               IF (I .EQ. J) A(I,I) = 1.0D0
               IF (I .NE. J) A(I,J) = 0.0D0
  260       CONTINUE
  270    CONTINUE
      GO TO 470
C
C     KASE 13
C
  280 CONTINUE
         N = 6
         WRITE (IWRITE,290) KASE,N
  290    FORMAT (5H KASE, I3, 3X, 16HUPPER TRIANGULAR / 4H N =, I4)
         DO 310 I = 1, N
            DO 300 J = 1, N
               IF (I .GT. J) A(I,J) = 0.0D0
               IF (I .LE. J) A(I,J) = DBLE(FLOAT(J-I+1))
  300       CONTINUE
  310    CONTINUE
      GO TO 470
C
C     KASE 14
C
  320 CONTINUE
         N = 6
         WRITE (IWRITE,330) KASE,N
  330    FORMAT (5H KASE, I3, 3X, 16HLOWER TRIANGULAR / 4H N =, I4)
         DO 350 I = 1, N
            DO 340 J = 1, N
               IF (I .LT. J) A(I,J) = 0.0D0
               IF (I .GE. J) A(I,J) = DBLE(FLOAT(I-J+1))
  340       CONTINUE
  350    CONTINUE
      GO TO 470
C
C     KASE 15
C
  360 CONTINUE
         N = 5
         WRITE (IWRITE,370) KASE,N
  370    FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY =D1MACH(1)*DBLE(FLOAT(N*N))*100.D0
         WRITE (IWRITE,380) TINY
  380    FORMAT (14H TINY        =, 1PD15.5)
         DO 400 I = 1, N
            DO 390 J = 1, N
               A(I,J) = TINY*DBLE(FLOAT(J))/DBLE(FLOAT(MAX0(I,J)))
  390       CONTINUE
  400    CONTINUE
      GO TO 470
C
C     KASE 16
C
  410 CONTINUE
         N = 5
         WRITE (IWRITE,420) KASE,N
  420    FORMAT (5H KASE, I3, 3X, 16HNEAR OVERFLOW    / 4H N =, I4)
         HUGE =D1MACH(2)/DBLE(FLOAT(N*N))
         WRITE (IWRITE,430) HUGE
  430    FORMAT (14H HUGE        =, 1PD15.5)
         DO 450 I = 1, N
            DO 440 J = 1, N
               A(I,J) = HUGE*(DBLE(FLOAT(J))/DBLE(FLOAT(MAX0(I,J))))
  440       CONTINUE
  450    CONTINUE
      GO TO 470
C
  460 CONTINUE
         N = 0
  470 CONTINUE
      RETURN
C
      END
