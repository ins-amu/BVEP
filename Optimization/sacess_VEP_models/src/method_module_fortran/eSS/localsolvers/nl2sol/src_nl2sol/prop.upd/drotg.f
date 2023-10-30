      SUBROUTINE DROTG(SA, SB, SC, SS)
      DOUBLE PRECISION SA, SB, SC, SS
      DOUBLE PRECISION ONE, XR, YR, ZERO, DSQRT
      DATA ZERO/0.0D0/
      DATA ONE/1.0D0/
C COMPUTE.. MATRIX ( SC SS) SO  ( SC SS)(SA) = (DSQRT(SA**2+SB**2))
C                  (-SS SC)     (-SS SC)(SB)   (   0             )
C  DSQRT(SA**2+SB**2) REPLACES SA IN STORAGE.
C  ZERO REPLACES SB
      IF (DABS(SA) .LE. DABS(SB)) GOTO 1
         XR = SB/SA
         YR = DSQRT(ONE+XR**2)
         SC = DSIGN(ONE/YR, SA)
         SS = SC*XR
         SA = DABS(SA)*YR
         GOTO  4
   1     IF (SB .EQ. ZERO) GOTO 2
            XR = SA/SB
            YR = DSQRT(ONE+XR**2)
            SS = DSIGN(ONE/YR, SB)
            SC = SS*XR
            SA = DABS(SB)*YR
            GOTO  3
   2        SC = ONE
            SS = ZERO
   3  CONTINUE
   4  SB = ZERO
      RETURN
      END
