C$TEST LYRC
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LYRC
C***********************************************************************
C
C  TEST OF THE PORT PROGRAMS CHECE AND FRIENDS
C
C***********************************************************************
C     MAIN PROGRAM
      INTEGER IWRITE
C     ALLOW 5000 UNDERFLOWS.
C
C     OUTPUT UNIT NUMBER
C
      IWRITE = I1MACH(2)
C
      CALL SPOTS(IWRITE)
      STOP
      END
      SUBROUTINE SPOTS(IWRITE)
C     IWRITE IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        CHECE;CHEFBS;CHEML;CHELE;CBACE;CBAFS;CBABS;CBAML;CBALE
C
C
C     SUBROUTINES AND FUNCTIONS
C
C     PORT CHECE,CHEFBS,CHEML,CHELE,CBACE,CBAFS,CBABS,CBAML,CBALE
C     EXTERNAL CSIXX,R1MACH
C     BLAS CAXPY,SCASUM
C     FORTRAN CABS1,AMAX1,FLOAT,MAX0
C
C     INTERNAL VARIABLES
C
      COMPLEX APSAVE(120),AB(35,15),AINV(15,15),ASAVE(15,15)
      COMPLEX AP(120),B(15),X(15),XB(15),XEXACT(15)
      COMPLEX ABSAVE(35,15),AL(15,15)
      COMPLEX XP(15),T,Z(15)
      REAL ANORM,AINORM,COND,COND1,AIBNO,EBNORM
      REAL EN,ENORM,EPS,FNORM,Q(10),QS(10),RCOND,RCONDB
      REAL RCONDP,RNORM,SCASUM,R1MACH,XNORM
      COMPLEX AINVB(15,15)
      INTEGER IPVT(15),IPVTS(15)
      INTEGER I,IQ(10),I1,J,JB
      INTEGER K,KASE,KB,KBFAIL,KNPD,KOUNT,KPFAIL
      INTEGER KSUSP(10),LDA,IWRITE,M,N,NPRINT
      LOGICAL KBF,KPF
C
      LDA = 15
      LDAB= 35
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT = 3
C
      WRITE (IWRITE,560)
      WRITE (IWRITE,1000)
C
      DO 10 I = 1,10
         KSUSP(I) = 0
   10 CONTINUE
      KNPD = 0
      KPFAIL = 0
      KBFAIL = 0
C
C     SET EPS TO ROUNDING UNIT FOR REAL ARITHMETIC
C
      EPS = R1MACH(4)
      WRITE (IWRITE,570) EPS
      WRITE (IWRITE,550)
C
        CALL ENTER(1)
C     START MAIN LOOP
C
      KASE = 1
   20 CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL CSIXX(ASAVE,LDA,N,KASE,IWRITE)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 540
         ANORM = 0.0E0
         DO 30 J = 1, N
            ANORM = AMAX1(ANORM,SCASUM(N,ASAVE(1,J),1))
   30    CONTINUE
         WRITE (IWRITE,720) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (IWRITE,550)
            DO 40 I = 1, N
               WRITE (IWRITE,760) (ASAVE(I,J), J = 1, N)
   40       CONTINUE
            WRITE (IWRITE,550)
   50    CONTINUE
C
C        GENERATE EXACT SOLUTION
C
         XEXACT(1) = (1.0E0,0.0E0)
         IF (N .GE. 2) XEXACT(2) = (0.0E0,1.0E0)
         IF (N .LE. 2) GO TO 70
            DO 60 I = 3, N
               XEXACT(I) = -XEXACT(I-2)
   60       CONTINUE
   70    CONTINUE
C
C
C      PUT INTO PACKED FORM
         K = 0
         DO 130 J = 1, N
            DO 120 I = J,N
               K = K + 1
               AP(K) = ASAVE(I,J)
               APSAVE(K)=AP(K)
  120       CONTINUE
  130    CONTINUE
         CALL CHECE(N,AP,IPVTS,RCONDP)
         IF (NERROR(IERS).NE.0) CALL ERROFF
C        FACTOR BAND FORM AND COMPARE
C
         KBF = .FALSE.
         M = 0
         DO 200 J = 1, N
            DO 190 I = 1, J
               IF (ASAVE(I,J) .NE. 0.0E0) M = MAX0(M,J-I)
  190       CONTINUE
  200    CONTINUE
C
         ML=M+1
         DO 220 J = 1, N
             I1=MIN0(N,J+M)
             I2=MAX0(1,J-M)
            DO 210 I = I2, I1
               K =ML+I-J
               AB(K,J) = ASAVE(J,I)
               ABSAVE(K,J)=AB(K,J)
  210       CONTINUE
  220    CONTINUE
        M=2*ML-1
         WRITE (IWRITE,840) M
        CALL CBACE(N,ML,M,AB,LDAB,AL,LDA,IPVT,MU,RCONDB)
           IF (NERROR(IERR).NE.0) CALL ERROFF
           IF(IERR+IERS.EQ.0) GO TO 230
             WRITE(IWRITE,580)
            WRITE (IWRITE,930) RCONDP,RCONDB
             GO TO 530
 230       CONTINUE
            WRITE(IWRITE,930)RCONDP,RCONDB
C
C           COMPUTE INVERSE AND COND1 = TRUE CONDITION
C
            DO 290 J = 1, N
               DO 280 I = 1, N
                   AINV(I,J)=(0.0,0.0)
                   AINVB(I,J)=(0.0,0.0)
  280          CONTINUE
               AINV(J,J)=(1.0,0.0)
               AINVB(J,J)=(1.0,0.0)
  290       CONTINUE
           CALL CBAFS(N,ML,AL,LDA,IPVT,AINVB,LDA,N)
           CALL CBABS(N,AB,LDAB,AINVB,LDA,N,MU)
           CALL CHEFBS(N,AP,AINV,LDA,N,IPVTS)
           AINORM=0.0
           DO 310 J=1,N
              AIS=SCASUM(N,AINV(1,J),1)
               AINORM = AMAX1(AINORM,AIS)
  310       CONTINUE
            COND1 = ANORM*AINORM
            WRITE (IWRITE,600) COND1
C
C           GENERATE RIGHT HAND SIDE FOR BOTH SYMMETRIC AND BAND
C
            CALL CHEML(N,APSAVE,XEXACT,B)
            CALL MOVEFD(N,B,X)
            CALL CBAML(N,ML,M,ABSAVE,LDAB,XEXACT,XB)
C           SOLVE A*X = B
C
            CALL CHELE(N,APSAVE,X,N,1)
            IF (NERROR(IRE).NE.0) CALL ERROFF
            CALL CBALE(N,ML,M,ABSAVE,LDAB,XB,N,1)
            IF (IRE+NERROR(IRB).EQ.0) GO TO 311
               IF (IRB.NE.0) CALL ERROFF
               WRITE(IWRITE,580)
               GO TO 530
  311       CONTINUE
C
            IF (N .GT. NPRINT) GO TO 330
               WRITE (IWRITE,610)
               DO 320 I = 1, N
                  WRITE (IWRITE,790) X(I), XB(I)
  320          CONTINUE
               WRITE (IWRITE,550)
  330       CONTINUE
C
C
C           COMPUTE ERRORS AND RESIDUALS
C              E  =  X - XEXACT
C              EB =  XB - XEXACT
C              R  =  B - A*X
C
            XNORM = SCASUM(N,X,1)
            ENORM = 0.0E0
            EBNORM = 0.E0
            DO 460 J = 1, N
               ENORM = ENORM + CABS1(X(J)-XEXACT(J))
               EBNORM = EBNORM + CABS1(XB(J)-XEXACT(J))
               T = -X(J)
               CALL CAXPY(N,T,ASAVE(1,J),1,B,1)
  460       CONTINUE
            RNORM = SCASUM(N,B,1)
C
C           A*INV(A) - I
C
            AINORM = 0.0E0
            AIBNO = 0.E0
            DO 490 J = 1, N
               DO 470 I = 1, N
                  B(I) =(0.0E0, 0.0E0)
                 XB(I) =(0.0E0, 0.0E0)
  470          CONTINUE
               DO 480 K = 1, N
                  T = AINV(K,J)
                  CALL CAXPY(N,T,ASAVE(1,K),1,B,1)
                  T=AINVB(K,J)
                  CALL CAXPY(N,T,ASAVE(1,K),1,XB,1)
  480          CONTINUE
               B(J) = B(J) - (1.0E0,0.0E0)
               XB(J) = XB(J) - (1.0E0,0.0E0)
               AINORM = AMAX1(AINORM,SCASUM(N,B,1))
               AIBNO = AMAX1(AIBNO,SCASUM(N,XB,1))
  490       CONTINUE
C
            WRITE (IWRITE,620) ENORM, EBNORM
            WRITE (IWRITE,630) RNORM
            WRITE (IWRITE,740) AINORM,AIBNO
C
C           COMPUTE TEST RATIOS
C
            Q(1) = RCONDP/COND1
            Q(2) = RCONDB/COND1
            Q(3) = COND1/RCONDP
            Q(4) = COND1/RCONDB
            Q(5) = ENORM/(EPS*RCONDP*XNORM)
            Q(6) = EBNORM/(EPS*RCONDP*XNORM)
            DENOM=AMAX1(100.E0*R1MACH(1),EPS*ANORM*XNORM)
            Q(7) = RNORM/DENOM
            Q(8) = AINORM/(EPS*RCONDP)
            Q(9) = AIBNO/(EPS*RCONDP)
            WRITE (IWRITE,550)
            WRITE (IWRITE,640)
            WRITE (IWRITE,550)
            WRITE (IWRITE,690)
            WRITE (IWRITE,700)
            WRITE (IWRITE,710)
            WRITE (IWRITE,550)
            WRITE (IWRITE,750) (Q(I), I = 1, 9)
            WRITE (IWRITE,550)
C
C           LOOK FOR SUSPICIOUS RATIOS
C
            QS(1) = 1.0E0 + 4.0E0*EPS
            QS(2) = QS(1)
            QS(3) = 10.0E0
            QS(4) =QS(3)
            EN = FLOAT(N)
            IF (N .EQ. 1) EN = 2.0E0
            DO 500 I=5,9
               QS(I) = EN
  500       CONTINUE
            KOUNT = 0
            DO 520 I = 1, 9
               IQ(I) = 0
               IF (Q(I) .LE. QS(I)) GO TO 510
                  IQ(I) = 1
                  KSUSP(I) = KSUSP(I) + 1
                  KOUNT = KOUNT + 1
  510          CONTINUE
  520       CONTINUE
            IF (KOUNT .EQ. 0) WRITE (IWRITE,980)
            IF (KOUNT .NE. 0) WRITE (IWRITE,990) (IQ(I), I = 1,9)
            WRITE (IWRITE,550)
  530    CONTINUE
C
         WRITE (IWRITE,650)
         KASE = KASE + 1
      GO TO 20
  540 CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (IWRITE,660)
      KASE = KASE - 1
      WRITE (IWRITE,670) KASE
      WRITE (IWRITE,680) KSUSP
      WRITE (IWRITE,910)
      RETURN
C
C     MOST FORMATS, ALSO SOME IN CSIXX
C
  550 FORMAT (1H )
 560  FORMAT(22HPORT TESTER,CHE**CBA**)
  570 FORMAT ( / 14H EPSILON     =, 1PE17.5)
  580 FORMAT ( / 16H MAYBE SINGULAR. /)
  600 FORMAT (14H ACTUAL COND =, 1PE17.5)
  610 FORMAT ( / 4H X =)
  620 FORMAT (14H ERROR NORM  =, 1P2E17.5)
  630 FORMAT (14H RESID NORM  =, 1P1E17.5)
  640 FORMAT (26H TEST RATIOS.. E = EPSILON)
  650 FORMAT ( / 14H ************* /)
  660 FORMAT (8H1SUMMARY)
  670 FORMAT (18H NUMBER OF TESTS =, I4)
  680 FORMAT ( / 30H NUMBER OF SUSPICIOUS RATIOS =, 10I4)
  690 FORMAT( 42H    COND  COND(B)  ACTUAL  ACTUAL  ERROR  ,
     1        34HERROR(B)  RESID   A*AI-I A*AI-I(B))
  700 FORMAT (10(8H   -----))
 710  FORMAT(42H    ACTUAL ACTUAL   COND  COND(B) E*COND*X,
     1       34H E*COND*X E*A*X    E*COND  E*COND )
  720 FORMAT (14H NORM(A)     =, 1PE17.5)
  730 FORMAT (14H NORM(A-RT*R)=, 1PE17.5)
  740 FORMAT (14H NORM(A*AI-I)=, 1P2E17.5)
  750 FORMAT (10(1X, F7.2))
 760  FORMAT(1H ,E17.4)
  780 FORMAT (2G14.6)
  790 FORMAT (2G14.6)
  830 FORMAT ( / 28H BAND ROUTINES DO NOT AGREE,)
  840 FORMAT (5H M  =, I2)
  910 FORMAT ( / 12H END OF TEST)
  930 FORMAT (8H RCOND =, 1P3E17.5)
  980 FORMAT (21H NO SUSPICIOUS RATIOS)
  990 FORMAT (I8, 5I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
 1000 FORMAT (29H THIS VERSION DATED 09/21/78.)
      END
      SUBROUTINE CSIXX(A,LDA,N,KASE,IWRITE)
      INTEGER LDA,N,KASE,IWRITE
      COMPLEX A(LDA,1)
C
C     GENERATES COMPLEX SYMMETRIC INDEFINITE TEST MATRICES
C
C     EXTERNAL R1MACH
C     FORTRAN ABS,AIMAG,CABS,CMPLX,FLOAT,IABS,REAL
      COMPLEX T
      REAL TINY,HUGE,R1MACH
      INTEGER I,J
      REAL CABS1
      COMPLEX ZDUM
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
      GO TO (10, 10, 10, 50, 50, 70, 70, 70, 110, 150, 190, 230, 250,
     *       290, 330, 330, 380, 430, 480), KASE
C
C     KASE 1, 2 AND 3
C
   10 CONTINUE
         N = 5*KASE
         WRITE (IWRITE,20) KASE,N
   20    FORMAT (5H KASE, I3, 3X, 16HCOMPLEX HILBERT  / 4H N =, I4)
         T = (1.0E0,1.0E0)
         T = T/CABS(T)
         DO 40 J = 1, N
            DO 30 I = 1, J
               A(I,J) = T**(J - I)/CMPLX(FLOAT(I+J-1),0.0E0)
              A(J,I) = CONJG( A(I,J) )
C              FOR NON-C0MPLEX MATRICES, A(I,J) = 1.0/FLOAT(I+J-1)
   30       CONTINUE
   40    CONTINUE
      GO TO 490
C
C     KASE 4 AND 5
C
   50 CONTINUE
         N = 1
         WRITE (IWRITE,60) KASE,N
   60    FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = (3.0E0,0.0E0)
         IF (KASE .EQ. 5) A(1,1) = (0.0E0,0.0E0)
      GO TO 490
C
C     KASE 6, 7 AND 8
C
   70 CONTINUE
         N = 15
         WRITE (IWRITE,80) KASE,N
   80    FORMAT (5H KASE, I3, 3X, 16HTRIDIAGONAL      / 4H N =, I4)
         T = (1.0E0,0.0E0)
         IF (KASE .EQ. 7) T = (3.0E0,1.0E0)
         IF (KASE .EQ. 8) T = (100.0E0,100.0E0)
         DO 100 J = 1, N
            DO 90 I = 1, J
               A(I,J) = (0.0E0,0.0E0)
               IF (I .EQ. J) A(I,I) = (4.0E0,0.0E0)
               IF (I .EQ. J - 1) A(I,J) = T
              A(J,I) = CONJG( A(I,J) )
   90       CONTINUE
  100    CONTINUE
      GO TO 490
C
C     KASE 9
C
  110 CONTINUE
         N = 5
         WRITE (IWRITE,120) KASE,N
  120    FORMAT (5H KASE, I3, 3X, 16HPENTADIAGONAL    / 4H N =, I4)
         DO 140 J = 1, N
            DO 130 I = 1, J
               A(I,J) = (0.0E0,0.0E0)
               IF (I .GE. J - 2)
     *            A(I,J) = CMPLX((5.0E0-FLOAT(IABS(I-J)))**(10-I-J),
     *                           0.0E0)
              A(J,I) = CONJG( A(I,J) )
  130       CONTINUE
  140    CONTINUE
      GO TO 490
C
C     KASE 10
C
  150 CONTINUE
         N = 6
         WRITE (IWRITE,160) KASE,N
  160    FORMAT (5H KASE, I3, 3X, 16HTRIDIAG INVERSE  / 4H N =, I4)
         DO 180 J = 1, N
            DO 170 I = 1, J
               A(I,J) = CMPLX(FLOAT(N+1-J),0.0E0)
              A(J,I) = CONJG( A(I,J) )
  170       CONTINUE
  180    CONTINUE
      GO TO 490
C
C     KASE 11
C
  190 CONTINUE
         N = 10
         WRITE (IWRITE,200) KASE,N
  200    FORMAT (5H KASE, I3, 3X, 16HZERO DIAGONAL    / 4H N =, I4)
         DO 220 J = 1, N
            DO 210 I = 1, J
               A(I,J) = (0.0E0,0.0E0)
               IF (I .EQ. J - 1) A(I,J) = (1.0E0,1.0E0)
              A(J,I) = CONJG( A(I,J) )
  210       CONTINUE
  220    CONTINUE
      GO TO 490
C
C     KASE 12
C
  230 CONTINUE
         N = 2
         WRITE (IWRITE,240) KASE,N
  240    FORMAT (5H KASE, I3, 3X, 16HTWO BY TWO       / 4H N =, I4)
         A(1,1) = (4.0E0,0.0E0)
         A(1,2) = (1.0E0,2.0E0)
         A(2,1) = CONJG(A(1,2))
         A(2,2) = (0.0E0,0.0E0)
      GO TO 490
C
C     KASE 13
C
  250 CONTINUE
         N = 6
         WRITE (IWRITE,260) KASE,N
  260    FORMAT (5H KASE, I3, 3X, 16H ZERO MATRIX     / 4H N =, I4)
         DO 280 J = 1, N
            DO 270 I = 1, J
               A(I,J) = (0.0E0,0.0E0)
              A(J,I) = CONJG( A(I,J) )
  270       CONTINUE
  280    CONTINUE
      GO TO 490
C
C     KASE 14
C
  290 CONTINUE
         N = 3
         WRITE (IWRITE,300) KASE,N
  300    FORMAT (5H KASE, I3 / 4H N =, I4)
         DO 320 I = 1, N
            DO 310 J = 1, N
               A(I,J) = (0.0E0,0.0E0)
  310       CONTINUE
  320    CONTINUE
         A(1,3) = (1.0E0,0.0E0)
         A(3,1) = (1.0E0,0.0E0)
      GO TO 490
C
C     KASE 15 AND 16
C
  330 CONTINUE
         N = 15
         WRITE (IWRITE,340) KASE,N
  340    FORMAT (5H KASE, I3, 3X, 16H                 / 4H N =, I4)
         DO 370 J = 1, N
            DO 360 I = 1, J
               A(I,J) = (-1.0E0,1.0E0)
              A(J,I) = CONJG( A(I,J) )
               IF (I .NE. J) GO TO 350
                  IF (KASE .EQ. 15) A(I,I) = (26.0E0,0.0E0)
                  IF (KASE .EQ. 16) A(I,I) = CMPLX(FLOAT(I),0.0E0)
  350          CONTINUE
  360       CONTINUE
  370    CONTINUE
      GO TO 490
C
C     KASE 17
C
  380 CONTINUE
         N = 5
         WRITE (IWRITE,390) KASE,N
  390    FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY = R1MACH(1)*FLOAT(100*N*N)
         WRITE (IWRITE,400) TINY
  400    FORMAT (14H TINY        =, 1PE17.5)
         DO 420 J = 1, N
            DO 410 I = 1, J
               A(I,J) = TINY
     *                  *CMPLX(FLOAT(IABS(I-J))/FLOAT(I+J),FLOAT(I-J))
              A(J,I) = CONJG( A(I,J) )
  410       CONTINUE
  420    CONTINUE
      GO TO 490
C
C     KASE 18
C
  430 CONTINUE
         N = 5
         WRITE (IWRITE,440) KASE,N
  440    FORMAT (5H KASE, I3, 3X, 16HNEAR OVERFLOW    / 4H N =, I4)
         HUGE = SQRT(R1MACH(2))/FLOAT(100*N*N)
         WRITE (IWRITE,450) HUGE
  450    FORMAT (14H HUGE        =, 1PE17.5)
         DO 470 J = 1, N
            DO 460 I = 1, J
               A(I,J) = HUGE
     *                  *CMPLX(FLOAT(IABS(I-J))/FLOAT(I+J),FLOAT(I-J))
              A(J,I) = CONJG( A(I,J) )
  460       CONTINUE
  470    CONTINUE
      GO TO 490
C
  480 CONTINUE
         N = 0
  490 CONTINUE
      RETURN
C
      END
