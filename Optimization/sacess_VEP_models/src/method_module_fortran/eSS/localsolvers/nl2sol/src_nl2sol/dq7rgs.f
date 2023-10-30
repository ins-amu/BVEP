      SUBROUTINE DQ7RGS(IERR, IPIVOT, L, N, NN, NOPIVK, P, Q, R, W)
C
C  ***  COMPUTE QR FACTORIZATION VIA MODIFIED GRAM-SCHMIDT PROCEDURE
C  ***  WITH COLUMN PIVOTING  ***
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER IERR, L, N, NN, NOPIVK, P
      INTEGER IPIVOT(P)
      DOUBLE PRECISION Q(NN,P), R(1), W(P)
C     DIMENSION R(P*(P+1)/2)
C
C----------------------------  DESCRIPTION  ----------------------------
C
C        THIS ROUTINE COMPUTES COLUMNS  L  THROUGH  P  OF A QR FACTORI-
C     ZATION OF THE MATRIX  A  THAT IS ORIGINALLY STORED IN COLUMNS  L
C     THROUGH  P  OF  Q.  IT IS ASSUMED THAT COLUMNS 1 THROUGH  L-1  OF
C     THE FACTORIZATION HAVE ALREADY BEEN STORED IN  Q  AND  R.  THIS
C     CODE USES THE MODIFIED GRAM-SCHMIDT PROCEDURE WITH REORTHOGONALI-
C     ZATION AND, IF  NOPIVK  ALLOWS IT, WITH COLUMN PIVOTING -- IF
C     K .GT. NOPIVK,  THEN ORIGINAL COLUMN  K  IS ELIGIBLE FOR PIVOTING.
C     IF  IPIVOT(L) = 0  ON INPUT, THEN  IPIVOT  IS INITIALIZED SO THAT
C     IPIVOT(I) = I  FOR  I = L,...,P.  WHATEVER THE ORIGINAL VALUE OF
C     IPIVOT(L), THE CORRESPONDING ELEMENTS OF  IPIVOT  ARE INTERCHANGED
C     WHENEVER COLUMN PIVOTING OCCURS.  THUS IF  IPIVOT(L) = 0  ON IN-
C     PUT, THEN THE  Q  AND  R  RETURNED ARE SUCH THAT COLUMN  I  OF
C     Q*R  EQUALS COLUMN  IPIVOT(I)  OF THE ORIGINAL MATRIX  A.  THE UP-
C     PER TRIANGULAR MATRIX  R  IS STORED COMPACTLY BY COLUMNS, I.E.,
C     THE OUTPUT VECTOR  R  CONTAINS  R(1,1), R(1,2), R(2,2), R(1,3),
C     R(2,3), ..., R(P,P) (IN THAT ORDER).  IF ALL GOES WELL, THEN THIS
C     ROUTINE SETS  IERR = 0.  BUT IF (PERMUTED) COLUMN  K  OF  A  IS
C     LINEARLY DEPENDENT ON (PERMUTED) COLUMNS 1,2,...,K-1, THEN  IERR
C     IS SET TO  K AND THE R MATRIX RETURNED HAS  R(I,J) = 0  FOR
C     I .GE. K  AND  J .GE. K.  IN THIS CASE COLUMNS  K  THROUGH  P
C     OF THE  Q  RETURNED ARE NOT ORTHONORMAL.  W IS A SCRATCH VECTOR.
C        THE ORIGINAL MATRIX  A  AND THE COMPUTED ORTHOGONAL MATRIX  Q
C     ARE  N BY P  MATRICES.  NN  IS THE LEAD DIMENSION OF THE ARRAY  Q
C     AND MUST SATISFY  NN .GE. N.  NO PARAMETER CHECKING IS DONE.
C
C        CODED BY DAVID M. GAY (FALL 1979, SPRING 1984).
C
C--------------------------  LOCAL VARIABLES  --------------------------
C
      INTEGER I, II, J, K, KK, KM1, KP1, LM1
      LOGICAL IPINIT
      DOUBLE PRECISION AK, SINGTL, T, T1, T2, WK
      EXTERNAL DD7TPR, DR7MDC,DV2AXY, DV7SCP,DV7SWP, DV2NRM, DV7SCL
      DOUBLE PRECISION DD7TPR, DR7MDC, DV2NRM
C/+
      DOUBLE PRECISION DSQRT
C/
      DOUBLE PRECISION BIG, MEPS10, ONE, REOTOL, TEN, TINY, WTOL, ZERO
C/6
C     DATA ONE/1.0D+0/, REOTOL/0.25D+0/, TEN/1.D+1/, WTOL/0.75D+0/,
C    1     ZERO/0.0D+0/
C/7
      PARAMETER (ONE=1.0D+0, REOTOL=0.25D+0, TEN=1.D+1, WTOL=0.75D+0,
     1           ZERO=0.0D+0)
      SAVE MEPS10, TINY
C/
      DATA MEPS10/0.0D+0/, TINY/0.0D+0/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      IERR = 0
      IF (MEPS10 .GT. ZERO) GO TO 10
         MEPS10 = TEN * DR7MDC(3)
         TINY = DR7MDC(1)
         BIG = DR7MDC(6)
         IF (TINY*BIG .LT. ONE) TINY = ONE / BIG
 10   SINGTL = FLOAT(MAX0(N,P)) * MEPS10
      LM1 = L - 1
      J = L*LM1/2
      KK = J
      IPINIT = IPIVOT(L) .EQ. 0
C
C  ***  INITIALIZE W, IPIVOT, DIAG(R), AND R(I,J) FOR I = 1,2,...,L-1
C  ***  AND J = L,L+1,...,P.
C
      DO 50 I = L, P
         IF (IPINIT) IPIVOT(I) = I
         T = DV2NRM(N, Q(1,I))
         IF (T .GT. ZERO) GO TO 20
              W(I) = ONE
              J = J + LM1
              GO TO 40
 20      W(I) = ZERO
         IF (LM1 .EQ. 0) GO TO 40
              DO 30 K = 1, LM1
                   J = J + 1
                   T1 = DD7TPR(N, Q(1,K), Q(1,I))
                   R(J) = T1
                   CALL DV2AXY(N, Q(1,I), -T1, Q(1,K), Q(1,I))
                   W(I) = W(I) + (T1/T)**2
 30                CONTINUE
 40      J = J + I - LM1
         R(J) = T
 50      CONTINUE
C
C  ***  MAIN LOOP  ***
C
      DO 140 K = L, P
         KK = KK + K
         KP1 = K + 1
         IF (K .LE. NOPIVK) GO TO 70
         IF (K .GE. P) GO TO 70
C
C        ***  FIND COLUMN WITH MINIMUM WEIGHT LOSS  ***
C
              T = W(K)
              IF (T .LE. ZERO) GO TO 70
              J = K
              DO 60 I = KP1, P
                   IF (W(I) .GE. T) GO TO 60
                        T = W(I)
                        J = I
 60                CONTINUE
              IF (J .EQ. K) GO TO 70
C
C             ***  INTERCHANGE COLUMNS K AND J  ***
C
                   I = IPIVOT(K)
                   IPIVOT(K) = IPIVOT(J)
                   IPIVOT(J) = I
                   W(J) = W(K)
                   W(K) = T
                   I = J*(J+1)/2
                   T1 = R(I)
                   R(I) = R(KK)
                   R(KK) = T1
                   CALL DV7SWP(N, Q(1,K), Q(1,J))
                   IF (K .LE. 1) GO TO 70
                        I = I - J + 1
                        J = KK - K + 1
                        CALL DV7SWP(K-1, R(I), R(J))
C
C        ***  COLUMN K OF Q SHOULD BE NEARLY ORTHOGONAL TO THE PREVIOUS
C        ***  COLUMNS.  NORMALIZE IT, TEST FOR SINGULARITY, AND DECIDE
C        ***  WHETHER TO REORTHOGONALIZE IT.
C
 70      AK = R(KK)
         IF (AK .LE. ZERO) GO TO 150
         T1 = AK
         R(KK) = ONE
         T2 = ONE
         WK = W(K)
C
C        *** SET T TO THE NORM OF (Q(K,K),...,Q(N,K))
C        *** AND CHECK FOR SINGULARITY.
C
 80      IF (WK .LT. WTOL) GO TO 90
            T = DV2NRM(N, Q(1,K))
            IF (T*T2 / AK .GT. SINGTL) GO TO 100
            GO TO 150
 90      T = DSQRT(ONE - WK)
         IF (T*T2 .LE. SINGTL) GO TO 150
         T = T * AK
C
 100     IF (T .LT. TINY) GO TO 150
         R(KK) = T * R(KK)
         CALL DV7SCL(N, Q(1,K), ONE/T, Q(1,K))
         IF (T/T1 .GE. REOTOL) GO TO 120
C
C     ***  REORTHOGONALIZE COLUMN K  ***
C
              AK = ONE
              T2 = T * T2
              WK = ZERO
              J = KK - K
              KM1 = K - 1
              DO 110 I = 1, KM1
                   J = J + 1
                   T = DD7TPR(N, Q(1,I), Q(1,K))
                   WK = WK + T*T
                   R(J) = R(J) + T*R(KK)
 110               CALL DV2AXY(N, Q(1,K), -T, Q(1,I), Q(1,K))
              T1 = ONE
              GO TO 80
C
C        ***  COMPUTE R(K,I) FOR I = K+1,...,P AND UPDATE Q  ***
C
 120     IF (K .GE. P) GO TO 999
         J = KK + K
         II = KK
         DO 130 I = KP1, P
              II = II + I
              T = DD7TPR(N, Q(1,K), Q(1,I))
              R(J) = T
              J = J + I
              CALL DV2AXY(N, Q(1,I), -T, Q(1,K), Q(1,I))
              T1 = R(II)
              IF (T1 .GT. ZERO)  W(I) = W(I) + (T/T1)**2
 130          CONTINUE
 140     CONTINUE
C
C  ***  SINGULAR Q  ***
C
 150  IERR = K
      KM1 = K - 1
      J = KK
      DO 160 I = K, P
         CALL DV7SCP(I-KM1, R(J), ZERO)
         J = J + I
 160     CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DQ7RGS FOLLOWS  ***
      END
