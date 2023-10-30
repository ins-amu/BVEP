      SUBROUTINE SROTG(SA, SB, SC, SS)
      REAL SA, SB, SC, SS
      REAL ONE, XR, YR, ZERO, ABS, SQRT
      REAL SIGN
      DATA ZERO/0./
      DATA ONE/1./
C COMPUTE.. MATRIX ( SC SS) SO  ( SC SS)(SA) = (SQRT(SA**2+SB**2))
C                  (-SS SC)     (-SS SC)(SB)   (   0             )
C  SQRT(SA**2+SB**2) REPLACES SA IN STORAGE.
C  ZERO REPLACES SB
      IF (ABS(SA) .LE. ABS(SB)) GOTO 1
         XR = SB/SA
         YR = SQRT(ONE+XR**2)
         SC = SIGN(ONE/YR, SA)
         SS = SC*XR
         SA = ABS(SA)*YR
         GOTO  4
   1     IF (SB .EQ. ZERO) GOTO 2
            XR = SA/SB
            YR = SQRT(ONE+XR**2)
            SS = SIGN(ONE/YR, SB)
            SC = SS*XR
            SA = ABS(SB)*YR
            GOTO  3
   2        SC = ONE
            SS = ZERO
   3  CONTINUE
   4  SB = ZERO
      RETURN
      END
