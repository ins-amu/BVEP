C$TEST LYSC
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LYSC
C***********************************************************************
C
C  SECOND TEST OF THE PORT PROGRAMS CHECE AND FRIENDS
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
C        CHECE;CHEFBS;CHEML;CHELE;CBPCE;CBPFS;CBPBS;CBPML;CBPLE
C
C
C     SUBROUTINES AND FUNCTIONS
C
C     PORT CHECE,CHEFBS,CHEML,CHELE,CBPCE,CBPFS,CBPBS,CBPML,CBPLE
C     EXTERNAL CPOXX,R1MACH
C     BLAS CAXPY,SCASUM
C     FORTRAN CABS1,AMAX1,FLOAT,MAX0
C
C     INTERNAL VARIABLES
C
      COMPLEX APSAVE(120),AB(15,15),AINV(15,15),ASAVE(15,15)
      COMPLEX AP(120),B(15),X(15),XB(15),XEXACT(15)
      COMPLEX ABSAVE(15,15)
      COMPLEX XP(15),T,Z(15)
      REAL ANORM,AINORM,COND,COND1,AIBNO,EBNORM
      REAL EN,ENORM,EPS,FNORM,Q(10),QS(10),RCOND,RCONDB
      REAL RCONDP,RNORM,SCASUM,R1MACH,XNORM
      COMPLEX AINVB(15,15),SC
      REAL AIS
      INTEGER IPVT(15),IPVTS(15)
      INTEGER I,IQ(10),I1,J,JB
      INTEGER K,KASE,KB,KBFAIL,KNPD,KOUNT,KPFAIL
      INTEGER KSUSP(10),LDA,IWRITE,M,N,NPRINT
      LOGICAL KBF,KPF
C
      LDA = 15
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
         CALL CPOXX(ASAVE,LDA,N,KASE,IWRITE)
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
         XEXACT(1) = (1.0,0.0)
         IF (N .GE. 2) XEXACT(2) = (0.0,1.0)
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
         IF (NERROR(IERR).NE.0) CALL ERROFF
C        FACTOR BAND FORM AND COMPARE
C
         KBF = .FALSE.
         M = 0
         DO 200 J = 1, N
            DO 190 I = 1, J
               IF (CABS1(ASAVE(I,J)) .NE. 0.0E0) M = MAX0(M,J-I)
  190       CONTINUE
  200    CONTINUE
         MP1 =M+1
C
         DO 220 J = 1, N
             I1=MIN0(N,J+M)
            DO 210 I = J, I1
               K =I-J+1
               AB(K,J) = ASAVE(I,J)
               ABSAVE(K,J)=AB(K,J)
  210       CONTINUE
  220    CONTINUE
         WRITE (IWRITE,840) M
        CALL CBPCE(N,MP1,AB,LDA,RCONDB)
           IF (NERROR(IERR).EQ.0) GO TO 230
          CALL ERROFF
          IF ((IERR).LT.10+N)GO TO 226
             WRITE(IWRITE,860)IERR
             J = IERR - N - 10
             WRITE(IWRITE,229)AB(1,J)
 229          FORMAT(19H OFFENDING DIAGONAL,E15.7)
             GO TO 530
 226         WRITE(IWRITE,580)
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
           CALL CBPFS(N,MP1,AB,LDA,AINVB,LDA,N)
            CALL CBPBS(N,MP1,AB,LDA,AINVB,LDA,N)
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
            CALL CBPML(N,MP1,ABSAVE,LDA,XEXACT,XB)
C           SOLVE A*X = B
C
            CALL CHELE(N,APSAVE,X,N,1)
            IF (NERROR(IRE).NE.0) CALL ERROFF
            CALL CBPLE(N,MP1,ABSAVE,LDA,XB,N,1)
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
C              F  =  A - TRANS(R)*R
C
            XNORM = SCASUM(N,X,1)
            ENORM = 0.0E-0
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
                  B(I) = (0.0,0.0)
                 XB(I) = (0.0,0.0)
  470          CONTINUE
               DO 480 K = 1, N
                  T = AINV(K,J)
                  CALL CAXPY(N,T,ASAVE(1,K),1,B,1)
                  T=AINVB(K,J)
                  CALL CAXPY(N,T,ASAVE(1,K),1,XB,1)
  480          CONTINUE
               B(J) = B(J) - (1.0,0.0)
               XB(J) = XB(J) - (1.0,0.0)
               AINORM = AMAX1(AINORM,SCASUM(N,B,1))
               AIBNO = AMAX1(AIBNO,SCASUM(N,XB,1))
  490       CONTINUE
            FNORM = 0.0E0
            ML=M+1
            NP1=N+1
            MLP1=ML+1
            DO 495 J=1,N
               NUMAX= MIN0(ML,NP1-J)
               JM1=J-1
               IEND=JM1+NUMAX
               DO 491 I=J,IEND
                  B(I)=-ASAVE(I,J)
 491           CONTINUE
               KBEGIN=MAX0(1,J-M)
               L=MIN0(J,ML)
               NUM=MLP1-L
               IF (JM1.LT.KBEGIN) GO TO 493
               DO 492 K=KBEGIN,JM1
                  SC=(AB(L,K))*AB(1,K)
                  IEND=JM1+NUM
                  LL=L
                  IF(IEND.LT.J) GO TO 4912
                  DO 4911 I=J,IEND
                      B(I)=B(I)+SC*CONJG(AB(LL,K))
                      LL=LL+1
 4911             CONTINUE
 4912             L=L-1
                  NUM=MIN0(NUM+1,NUMAX)
 492           CONTINUE
 493           SC=AB(1,J)
                  IEND=JM1+NUM
                  JP1=J+1
                  IF (IEND.LT.JP1) GO TO 4941
                  L=2
                  DO 494 I=JP1,IEND
                     B(I)=B(I)+SC*CONJG(AB(L,J))
                     L=L+1
 494              CONTINUE
 4941          B(J)=B(J)+AB(1,J)
               FNORM=AMAX1(FNORM,SCASUM(NUMAX,B(J),1))
 495       CONTINUE
C
            WRITE (IWRITE,620) ENORM, EBNORM
            WRITE (IWRITE,630) RNORM
            WRITE (IWRITE,730) FNORM
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
            DENOM=AMAX1(EPS*ANORM*XNORM,100.0*R1MACH(1))
            Q(7) = RNORM/DENOM
            DENOM=AMAX1(EPS*ANORM,100.0*R1MACH(1))
            Q(8) = FNORM/DENOM
            Q(9) = AINORM/(EPS*RCONDP)
            Q(10) = AIBNO/(EPS*RCONDP)
            WRITE (IWRITE,550)
            WRITE (IWRITE,640)
            WRITE (IWRITE,550)
            WRITE (IWRITE,690)
            WRITE (IWRITE,700)
            WRITE (IWRITE,710)
            WRITE (IWRITE,550)
            WRITE (IWRITE,750) (Q(I), I = 1, 10)
            WRITE (IWRITE,550)
C
C           LOOK FOR SUSPICIOUS RATIOS
C
            QS(1) = 1.0E0 + 4.0E0*EPS
            QS(2) = QS(1)
            QS(3) = 10.0E0
            QS(4) =QS(3)
            EN = DBLE(FLOAT(N))
            IF (N .EQ. 1) EN = 2.0E0
            DO 500 I=5,10
               QS(I) = EN
  500       CONTINUE
            KOUNT = 0
            DO 520 I = 1, 10
               IQ(I) = 0
               IF (Q(I) .LE. QS(I)) GO TO 510
                  IQ(I) = 1
                  KSUSP(I) = KSUSP(I) + 1
                  KOUNT = KOUNT + 1
  510          CONTINUE
  520       CONTINUE
            IF (KOUNT .EQ. 0) WRITE (IWRITE,980)
            IF (KOUNT .NE. 0) WRITE (IWRITE,990) (IQ(I), I = 1,10)
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
      WRITE (IWRITE,900) KNPD
      WRITE (IWRITE,680) KSUSP
      WRITE (IWRITE,910)
      RETURN
C
C     MOST FORMATS, ALSO SOME IN CPOXX
C
  550 FORMAT (1H )
 560  FORMAT(22HPORT TESTER,CHE**CBP**)
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
     1        40HERROR(B)  RESID  A-RT*R A*AI-I A*AI-I(B))
  700 FORMAT (10(8H   -----))
 710  FORMAT(42H    ACTUAL ACTUAL   COND  COND(B) E*COND*X,
     1       40H E*COND*X E*A*X    E*A   E*COND  E*COND )
  720 FORMAT (14H NORM(A)     =, 1PE17.5)
  730 FORMAT (14H NORM(A-RT*R)=, 1PE17.5)
  740 FORMAT (14H NORM(A*AI-I)=, 1P2E17.5)
  750 FORMAT (10(1X, F7.2))
  760 FORMAT(1H ,3E17.4)
  780 FORMAT (2G14.6)
  790 FORMAT (2G14.6)
  830 FORMAT ( / 28H BAND ROUTINES DO NOT AGREE,)
  840 FORMAT (5H M  =, I2)
  860 FORMAT (30H NOT POSITIVE DEFINITE, INFO =, I2)
  900 FORMAT (34H NUMBER OF NOT POSITIVE DEFINITE =, I4)
  910 FORMAT ( / 12H END OF TEST)
  930 FORMAT (8H RCOND =, 1P3E17.5)
  980 FORMAT (21H NO SUSPICIOUS RATIOS)
  990 FORMAT (I8, 5I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
 1000 FORMAT (29H THIS VERSION DATED 09/21/78.)
      END
      SUBROUTINE CPOXX(A,LDA,N,KASE,IWRITE)
      INTEGER LDA,N,KASE,IWRITE
      COMPLEX A(LDA,1)
C
C     GENERATES COMPLEX POSITIVE DEFINITE TEST MATRICES
C
C     EXTERNAL CMACH
C     FORTRAN CABS,CMPLX,CONJG,FLOAT,IABS,MAX0,MIN0
      COMPLEX T
      REAL TINY,HUGE,CMACH
      INTEGER I,J
C
      GO TO (10, 10, 10, 50, 50, 70, 70, 70, 120, 160, 200, 240, 290,
     *       340), KASE
C
C     KASE 1, 2 AND 3
C
   10 CONTINUE
         N = 5*KASE
         WRITE (IWRITE,20) KASE,N
   20    FORMAT (5H KASE, I3, 3X, 16HHILBERT          / 4H N =, I4)
         T = (1.0E0,1.0E0)
         T = T/CABS(T)
         DO 40 J = 1, N
            DO 30 I = 1, N
               A(I,J) = T**(J - I)/CMPLX(FLOAT(I+J-1),0.0E0)
C              FOR REAL MATRICES, A(I,J) = 1.0/FLOAT(I+J-1)
   30       CONTINUE
   40    CONTINUE
      GO TO 350
C
C     KASE 4 AND 5
C
   50 CONTINUE
         N = 1
         WRITE (IWRITE,60) KASE,N
   60    FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = (3.0E0,0.0E0)
         IF (KASE .EQ. 5) A(1,1) = (0.0E0,0.0E0)
      GO TO 350
C
C     KASE 6, 7 AND 8
C
   70 CONTINUE
         N = 15
         IF (KASE .NE. 8) WRITE (IWRITE,80) KASE,N
   80    FORMAT (5H KASE, I3, 3X, 16HTRIDIAGONAL      / 4H N =, I4)
         IF (KASE .EQ. 8) WRITE (IWRITE,90) KASE,N
   90    FORMAT (5H KASE, I3, 3X, 16HDIAGONAL         / 4H N =, I4)
         T = (1.0E0,0.0E0)
         IF (KASE .EQ. 7) T = (2.0E0,1.0E0)
         IF (KASE .EQ. 8) T = (0.0E0,0.0E0)
         DO 110 J = 1, N
            DO 100 I = 1, J
               A(I,J) = (0.0E0,0.0E0)
               IF (I .EQ. J) A(I,I) = (4.0E0,0.0E0)
               IF (I .EQ. J - 1) A(I,J) = T
               A(J,I) = CONJG(A(I,J))
  100       CONTINUE
  110    CONTINUE
      GO TO 350
C
C     KASE 9
C
  120 CONTINUE
         N = 5
         WRITE (IWRITE,130) KASE,N
  130    FORMAT (5H KASE, I3, 3X, 16HPENTADIAGONAL    / 4H N =, I4)
         DO 150 J = 1, N
            DO 140 I = 1, N
               A(I,J) = (0.0E0,0.0E0)
               IF (IABS(I-J) .LE. 2)
     *            A(I,J) = CMPLX((5.0E0-FLOAT(IABS(I-J)))**(10-I-J),
     *                           0.0E0)
  140       CONTINUE
  150    CONTINUE
      GO TO 350
C
C     KASE 10
C
  160 CONTINUE
         N = 6
         WRITE (IWRITE,170) KASE,N
  170    FORMAT (5H KASE, I3, 3X, 16HTRIDIAG INVERSE  / 4H N =, I4)
         DO 190 J = 1, N
            DO 180 I = 1, J
               A(I,J) = CMPLX(FLOAT(N+1-J),0.0E0)
               A(J,I) = A(I,J)
  180       CONTINUE
  190    CONTINUE
      GO TO 350
C
C     KASE 11
C
  200 CONTINUE
         N = 15
         WRITE (IWRITE,210) KASE,N
  210    FORMAT (5H KASE, I3, 3X, 16HTEST COND        / 4H N =, I4)
         DO 230 J = 1, N
            DO 220 I = 1, N
               IF (I .EQ. J) A(I,J) = CMPLX(FLOAT(I),0.0E0)
               IF (I .GT. J) A(I,J) = CMPLX(FLOAT(J-2),0.0E0)
               IF (I .LT. J) A(I,J) = CMPLX(FLOAT(I-2),0.0E0)
  220       CONTINUE
  230    CONTINUE
      GO TO 350
C
C     KASE 12
C
  240 CONTINUE
         N = 5
         WRITE (IWRITE,250) KASE,N
  250    FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY = 100. * R1MACH(1)
         WRITE (IWRITE,260) TINY
  260    FORMAT (14H TINY        =, 1PE17.5)
         DO 280 I = 1, N
            DO 270 J = 1, N
               A(I,J) = CMPLX(TINY*FLOAT(MIN0(I,J))/FLOAT(MAX0(I,J)),
     *                        0.0E0)
  270       CONTINUE
  280    CONTINUE
      GO TO 350
C
C     KASE 13
C
  290 CONTINUE
         N = 5
         WRITE (IWRITE,300) KASE,N
  300    FORMAT (5H KASE, I3, 3X, 16HNEAR OVERFLOW    / 4H N =, I4)
         HUGE= SQRT(R1MACH(2))/(FLOAT(100*N*N))
         WRITE (IWRITE,310) HUGE
  310    FORMAT (14H HUGE        =, 1PE17.5)
         DO 330 I = 1, N
            DO 320 J = 1, N
               A(I,J) = CMPLX(HUGE*FLOAT(MIN0(I,J))/FLOAT(MAX0(I,J)),
     *                        0.0E0)
  320       CONTINUE
  330    CONTINUE
      GO TO 350
C
  340 CONTINUE
         N = 0
  350 CONTINUE
      RETURN
C
      END
