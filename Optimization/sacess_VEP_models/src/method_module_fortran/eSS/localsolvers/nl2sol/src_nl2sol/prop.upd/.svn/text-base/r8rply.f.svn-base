      SUBROUTINE R8RPLY(NN, U, V, P, Q, A, B)
C DIVIDES P BY THE QUADRATIC  1,U,V  PLACING THE
C QUOTIENT IN Q AND THE REMAINDER IN A,B
      REAL P(NN), Q(NN), U, V, A, B, C
      INTEGER I
      B = P(1)
      Q(1) = B
      A = P(2) - U*B
      Q(2) = A
      DO 10 I=3,NN
        C = P(I) - U*A - V*B
        Q(I) = C
        B = A
        A = C
   10 CONTINUE
      RETURN
      END
