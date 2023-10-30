      SUBROUTINE R4RPLY(TYPE, UU, VV, P, K)
C COMPUTE NEW ESTIMATES OF THE QUADRATIC COEFFICIENTS
C USING THE SCALARS COMPUTED IN R2RPLY.
      COMMON /P66PLY/ SR, SI, U,
     1 V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     2 H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
C
      INTEGER N, NN
      INTEGER TYPE
      REAL ETA, ARE, MRE
C
      REAL P(1), K(1),
     1  SR, SI, U, V, A, B, C, D,
     2 A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     3 LZR, LZI
      REAL A4, A5, B1, B2, C1, C2, C3,
     1 C4, TEMP, UU, VV
C USE FORMULAS APPROPRIATE TO SETTING OF TYPE.
      IF (TYPE.EQ.3) GO TO 30
      IF (TYPE.EQ.2) GO TO 10
      A4 = A + U*B + H*F
      A5 = C + (U+V*F)*D
      GO TO 20
   10 A4 = (A+G)*F + H
      A5 = (F+U)*C + V*D
C EVALUATE NEW QUADRATIC COEFFICIENTS.
   20 B1 = -K(N)/P(NN)
      B2 = -(K(N-1)+B1*P(N))/P(NN)
      C1 = V*B2*A1
      C2 = B1*A7
      C3 = B1*B1*A3
      C4 = C1 - C2 - C3
      TEMP = A5 + B1*A4 - C4
      IF (TEMP.EQ.0.E0) GO TO 30
      UU = U - (U*(C3+C2)+V*(B1*A1+B2*A7))/TEMP
      VV = V*(1.+C4/TEMP)
      RETURN
C IF TYPE=3 THE QUADRATIC IS ZEROED
   30 UU = 0.E0
      VV = 0.E0
      RETURN
      END
