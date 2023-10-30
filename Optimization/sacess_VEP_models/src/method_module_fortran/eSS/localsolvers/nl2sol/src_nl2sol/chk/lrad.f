C$TEST LRAD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LRAD
C***********************************************************************
C
C  TEST OF THE PORT LINEAR PROGRAMMING - DLINPR, ETC.
C
C***********************************************************************
C  MAIN PROGRAM
      DOUBLE PRECISION DSTAK(8000)
      DOUBLE PRECISION ANULL, FEASIN, SN
      COMMON /CSTAK/ DSTAK
      EXTERNAL DLPMAN,  LPT,DLPRNT
      INTEGER IAG, IAS,  NERROR, I, J
      INTEGER K, L, M, N, KASE
      INTEGER IERR, IPTG(150), IA, IE, IL, IS
      INTEGER IROLD, ISIMP(200), ITMAX, ITYPE, IWRITE
      COMMON /KA/ KASE
      DOUBLE PRECISION CHANGE,  CTX, SUM, CHAGE2, CHAGE3
      DOUBLE PRECISION CHAGE4, CHAGE5, A(150, 50), B(150), C(150), X(
     1   100, 5)
      DOUBLE PRECISION Y(100), DDOT, CTXD, SIMP(200),  DASUM
      REAL UNI
      DOUBLE PRECISION D1MACH,Q(8),EPS
C  SET UP OUTPUT WRITE UNIT
C
       IWRITE = I1MACH(2)
      CALL ISTKIN(8000, 4)
      CALL ENTSRC(IROLD, 1)
      KASE = 0
          EPS=DSQRT(D1MACH(4))
         NOBAD=0
      IA = 150
 10   CONTINUE
      KASE = KASE+1
      CALL SETUP(A, IA, M, N, IE, IS, KASE, ISIMP, SIMP, ITYPE, B, C)
          DO 20 I=1,8
                 Q(I)=0.0D0
 20         CONTINUE
       IF (N.GT.0) GO TO 40
            WRITE(IWRITE,30)NOBAD
 30          FORMAT(26H NUMBER OF BAD EXAMPLES IS, I5)
                STOP
 40               CONTINUE
      WRITE(IWRITE,60)KASE
      WRITE (IWRITE,  50) M, N, IE, IS, ITYPE
 50   FORMAT (3H M=, I10, 2HN=, I10, 13H NO. EQUALITY, I10,
     1   11H NO. SIMPLE, I10, 6H ITYPE, I10)
 60   FORMAT(//6H KASE=,I5)
      DO  70 I = 1, N
         X(I, 1) = DBLE(UNI(0))
         X(I, 2) = X(I, 1)
 70      CONTINUE
      ITMAX = 5*N
         CALL DLINPA(A,M,N,DLPMAN,IA,B,C,X,ITMAX,CTX,IS,SIMP,
     1      ISIMP, IE, DLPRNT, IAG, IAS, IPTG, X(I, 5))
      IF (NERROR(IERR) .EQ. 0) GOTO 90
         WRITE (IWRITE,  80) IERR
 80      FORMAT (4H ERR, I10)
         CALL ERROFF
         IF (IERR.EQ.8) GO TO 170
         GOTO  10
 90       CONTINUE
           FEASIN=0.0D0
           IAGIE=IAG+IE
           ANULL=0.0D0
           DO 100 I=1,M
              JJ=IPTG(I)
              SUM=DDOT(N,A(JJ,1),IA,X,1)-B(JJ)
              FEASIN=DMAX1(FEASIN,-SUM)
              IF (I.LE.IAGIE)ANULL=DMAX1(ANULL,DABS(SUM))
 100       CONTINUE
      SN= DASUM(N,X,1)
      SN=DMAX1(SN,1.D0)
           FEASIN=FEASIN/SN
            ANULL=ANULL/SN
           Q(1)=FEASIN
           Q(2)=ANULL
           WRITE(IWRITE,110)FEASIN,ANULL
 110       FORMAT(22H MAX. INFEASABILIBITY=,E15.5,18H NULLITY VIOLATION,
     1     E15.5)
C SOLVE THE PROBLEM AGAIN
C
      DO  120 I = 1, N
         X(I, 3) = X(I, 1)
 120     CONTINUE
      CALL DLINPR(A, M, N, IA, B, C, X(1, 3), ITMAX, CTX, IS,
     1   SIMP, ISIMP, IE)
      IF (NERROR(IERR) .EQ. 0) GOTO 140
         WRITE (IWRITE,  130) IERR
 130     FORMAT (4H ERR, I10)
         CALL ERROFF
         GOTO  10
 140  DO  150 I = 1, N
         Y(I) = X(I, 3)-X(I, 1)
 150     CONTINUE
      CHANGE = DASUM(N, Y, 1)/SN
            Q(3)=CHANGE
      WRITE (IWRITE,  160) CHANGE
 160  FORMAT(51H SOLVING SECOND TIME WITH INITIAL GUESS AS SOLUTION/
     1  9H CHANGE= ,1PD23.15)
C
C DIFFERENT INITIAL POINT
C
 170      CONTINUE
      DO  180 I = 1, N
         X(I, 3) = -X(I, 2)
 180     CONTINUE
      CALL DLINPR(A, M, N, IA, B, C, X(1, 3), ITMAX, CTX, IS,
     1   SIMP, ISIMP, IE)
      IF (NERROR(IERR) .EQ. 0) GOTO 200
         WRITE (IWRITE,  190) IERR
 190     FORMAT (4H ERR, I10)
         CALL ERROFF
         GOTO  10
 200  DO  210 I = 1, N
         Y(I) = X(I, 3)-X(I, 1)
 210     CONTINUE
      CHAGE2=DASUM(N,Y,1)/SN
            Q(4)=CHAGE2
      WRITE (IWRITE,  220) CHAGE2
 220  FORMAT(44H INITIAL POINT IS NEGATED AND PROBLEM SOLVED/
     *9H CHANGE= ,1PD23.15)
C
C MAKE SURE HAVE INFEASIBLE STARTING POINT
C
      IF (M .LE. 0) GOTO 330
         SUM = B(1)-DDOT(N, A(1, 1), IA, X(1, 2), 1)
         IF (SUM .LE. 0D0) GOTO 240
            WRITE (IWRITE,  230)
 230        FORMAT (29H INITIAL GUESS WAS INFEASIBLE)
            GOTO  320
 240        I = 1
 250        IF (A(I, 1) .NE. 0.0D0) GOTO  260
               I = I+1
               GOTO  250
 260        DO  270 J = 1, N
               X(J, 3) = X(J, 2)
 270           CONTINUE
            X(I, 3) = X(I, 3)+(SUM+1.D0)/A(I, 1)
            CALL DLINPR(A, M, N, IA, B, C, X(1, 3), ITMAX, CTX,
     1    IS, SIMP, ISIMP, IE)
            IF (NERROR(IERR) .EQ. 0) GOTO 290
               WRITE (IWRITE,  280) IERR
 280           FORMAT (4H ERR, I10)
               CALL ERROFF
               GOTO  10
 290        DO  300 I = 1, N
               Y(I) = X(I, 3)-X(I, 1)
 300           CONTINUE
            CHAGE3 = DASUM(N, Y, 1)/SN
             Q(5)=CHAGE3
            WRITE (IWRITE,  310) CHAGE3
 310     FORMAT(43H INITIAL POINT MADE INFEASIBLE AND RESOLVED/
     *10H CHANGE = ,1PD23.15)
 320     CONTINUE
C
C SOLVE WITH SIMPLE CONVERTED INTO GENERAL
C
 330  IF (IS .EQ. 0) GOTO 430
         L = M
         DO  370 I = 1, IS
            L = L+1
            DO  340 J = 1, N
               A(L, J) = 0.0
 340           CONTINUE
            K = IABS(ISIMP(I))
            A(L, K) = -1.D0
            B(L) = SIMP(I)
            IF (ISIMP(I) .LE. 0) GOTO 350
               A(L, K) = 1.D0
               GOTO  360
 350           B(L) = -SIMP(I)
 360        CONTINUE
 370        CONTINUE
         DO  380 I = 1, N
            X(I, 3) = X(I, 2)
 380        CONTINUE
         IL = 0
         MM=M+IS
         CALL DLINPR(A,MM, N, IA, B, C, X(1, 3), ITMAX, CTX, IL
     1 , SIMP, ISIMP, IE)
         IF (NERROR(IERR) .EQ. 0) GOTO 400
            WRITE (IWRITE,  390) IERR
 390        FORMAT (4H ERR, I10)
            CALL ERROFF
            GOTO  10
 400     DO  410 I = 1, N
            Y(I) = X(I, 3)-X(I, 1)
 410        CONTINUE
         CHAGE4 = DASUM(N, Y, 1)/SN
             Q(6)=CHAGE4
         WRITE (IWRITE,  420) CHAGE4
 420    FORMAT(38H SIMPLE CONSTRAINTS WRITTEN AS GENERAL/
     1 9H CHANGE= ,1PD23.15)
C
C TEST DUALITY IF APPLICABLE
C
 430  IF (ITYPE .NE. 1) GOTO 520
       WRITE(IWRITE,440)
 440   FORMAT(17H TEST FOR DUALITY)
         DO  460 I = 1, M
             DO 450 J=1,N
                A(I,J)=-A(I,J)
 450          CONTINUE
            X(I, 4) = DBLE(UNI(0))
 460        CONTINUE
         CALL DLINPA(A,N,M,LPT,IA,C,B,X(1,4),ITMAX,CTXD,0,SIMP,
     1      ISIMP, 0, DLPRNT, IAG, IAS, IPTG, X(1, 3))
         IF (NERROR(IERR) .EQ. 0) GOTO 480
            WRITE (IWRITE,  470) IERR
 470        FORMAT (4H ERR, I10)
            CALL ERROFF
            GOTO  10
 480     DO  490 I = 1, M
            II=IPTG(I)
            Y(I) = X(I, 3)+X(II, 1)
 490        CONTINUE
         CHAGE5 = DASUM(M, Y, 1)/SN
         CTXD=-CTXD
              Q(7)=DABS(CTX-CTXD)/DABS(CTX)
             Q(8)=CHAGE5
         WRITE (IWRITE,  500) CTX, CTXD
 500     FORMAT (4H CTX, 1PD23.15, 4HCTXX, 1PD23.15)
         WRITE (IWRITE,  510) CHAGE5
 510     FORMAT (17H CHANGE FROM DUAL, 1PD23.15)
C
C
 520         CONTINUE
             IBAD=0
             DO 530 I=1,8
                  IF (Q(I).GT.EPS)IBAD=1
 530          CONTINUE
             IF (IBAD.EQ.0)WRITE(IWRITE,540)
 540        FORMAT(12H NO PROBLEMS)
              IF (IBAD.EQ.1)WRITE(IWRITE,550)KASE
 550         FORMAT(18H PROBLEM WITH CASE,I5)
             IF(IBAD.EQ.1)NOBAD=NOBAD+1
      GOTO  10
      END
      SUBROUTINE SETUP(A, IA, M, N, IE, IS, KASE, ISIMP, SIMP,
     1   ITYPE, B, C)
      INTEGER IA
      INTEGER M, N, IE, IS, KASE, ISIMP(1)
      INTEGER ITYPE
      DOUBLE PRECISION A(IA, 1), SIMP(1), B(1), C(1)
      INTEGER KAS, I, J, KA, JJ
      REAL UNI
      INTEGER TEMP
C THIS SUBROUTINE SETS UP VARIOUS PROBLEMS
      KAS = (KASE+1)/2
      N = 0
      IF (KAS.EQ.6)RETURN
      KA = 5
      IF (2*KAS .NE. KASE) KA = 1
      N = 5*KA
      M = N/2
       IF (KAS.EQ.2)M=M/2
C SET UP CONSTRAINT MATRIX
C
      DO  20 I = 1, M
          SUM=0.0
         DO  10 J = 1, N
            A(I, J) = 2.D0*DBLE(UNI(0)) - 1.D0
             SUM=SUM+A(I,J)
 10         CONTINUE
           B(I)=0.5D0*SUM
 20      CONTINUE
C
C SET UP COST FUNCTIONAL
C
      DO  30 I = 1, N
         C(I) = -DBLE(UNI(0))
 30      CONTINUE
C
C
C SET UP SIMPLE CONSTRAINTS
C
      JJ = N
      DO  40 I = 1, JJ
         SIMP(I) = 0.0D0
         ISIMP(I) = I
         TEMP = I+JJ
         SIMP(TEMP) =0.6D0
         TEMP = I+JJ
         ISIMP(TEMP) = -I
 40      CONTINUE
      GO TO (50,60,70,80,90),KAS
 50      ITYPE = 1
         IS = JJ
         IE = M
         GO TO 100
 60      ITYPE = 10
         IS = 2*JJ
         IE = M
         GO TO 100
 70      ITYPE = 2
         IS = JJ
         IE = 0
         GO TO 100
 80      ITYPE = 10
         IS = 2*JJ
         IE = 0
         GO TO 100
 90      ITYPE = 10
         IS = JJ
         IE = M/8
         GOTO  110
 100     CONTINUE
 110  RETURN
      END
       SUBROUTINE LPT(L,AA,IA,N,I,TVEC,T)
       LOGICAL L
       DOUBLE PRECISION AA(IA,N),TVEC(N),T,DDOT
       IF (L)GO TO 10
       CALL DCOPY(N,AA(1,I),1,TVEC,1)
       RETURN
 10    T=DDOT(N,AA(1,I),1,TVEC,1)
        RETURN
       END
