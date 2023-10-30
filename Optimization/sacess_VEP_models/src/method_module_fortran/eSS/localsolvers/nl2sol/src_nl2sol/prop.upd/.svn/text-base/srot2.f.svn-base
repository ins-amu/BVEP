      SUBROUTINE SROT2(N,X,INCX,Y,INCY,A,B)
C
C THIS SUBROUTINE IS SROT WITH THE OPTION OF GOING
C A VARIED AMOUNT IN INCREMENTS
C
        REAL X(1),Y(1)
        REAL A,B,W,V
        IF (N.LT.1)RETURN
        IX=IABS(INCX)
        IY=IABS(INCY)
        L=1
        K=1
        DO 10 I=1,N
           W=X(L)
           V=Y(K)
           X(L)=A*W+B*V
           Y(K)=-B*W+A*V
           L=L+IX
           K=K+IY
           IF (INCX.GT.0)IX=IX+1
           IF (INCY.GT.0)IY=IY+1
 10    CONTINUE
       RETURN
       END
