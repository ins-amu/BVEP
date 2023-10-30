C$TEST LLZA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LLZA
C***********************************************************************
C
C  TEST OF THE PORT PROGRAM LZ
C
C***********************************************************************
          COMPLEX A(20,20),B(20,20),X(20,20),EIGA(20),EIGB(20)
          COMPLEX ASAVE(20,20),BSAVE(20,20),EIG(20),EI
          INTEGER IP(20), I1MACH, IWRITE
          REAL EIGBN(20),EIGN(20),ABNORM,RESERR,RESUM
          COMPLEX SUMA,SUMB
          REAL AN,BN,ANORM,BNORM
          INTEGER MR(20,20),MI(20,20),LR(20),LI(20)
C
C THIS IS TEST PROGRAM FOR LZ
C
C
C SET UP OUTPUT WRITE UNIT
C
      IWRITE = I1MACH(2)
          NM=20
          WRITE(IWRITE,1)
 1        FORMAT(12H TEST FOR LZ)
          DO 100 ICASE=1,9
C SETUP PROBLEM
              CALL SETUP(NM,A,B,MR,MI,LR,LI,N,NA,NB,ICASE)
C
C SAVE MATRICES AND FIND THEIR NORM
C
              ANORM=0.0
              BNORM=0.0
              DO 10 I=1,N
                 AN=0.0
                 BN=0.0
                 DO 5 J=1,N
                    ASAVE(I,J)=A(I,J)
                    BSAVE(I,J)=B(I,J)
                    AN=AN+CABS(A(I,J))
                    BN=BN+CABS(B(I,J))
 5               CONTINUE
              ANORM=AMAX1(ANORM,AN)
              BNORM=AMAX1(BNORM,BN)
 10           CONTINUE
              ABNORM=ANORM+BNORM
C SOLVE GENERALIZED EIGENVALUE PROBLEM
              CALL LZ(N,A,NM,B,NM,X,NM,.TRUE.,EIGA,EIGB)
C TRY TO FIND ABSOLUTE ERROR WITHOUT DIVISION BY ZERO
              DO 20 I=1,N
                 EIGBN(I)=CABS(EIGB(I))
 20           CONTINUE
              CALL SRTPDR(EIGBN,1,IP,1,N)
              DO 30 I=1,N
                 IN=IP(I)
                 EIG(I)=EIGA(IN)/EIGB(IN)
                 EIGN(I)=CABS(EIG(I))
 30           CONTINUE
              CALL SRTPAR(EIGN,1,IP,1,NB)
              ERR=0.0
              DO 40 I=1,NB
                 IN=IP(I)
                 ERR=ERR+CABS(EIG(IN)-CMPLX(FLOAT(I),FLOAT(I)))
 40           CONTINUE
C PRINT EIGENVALUES FOR SMALL PROBLEMS
              IF (N.GT.5) GO TO 60
              WRITE(IWRITE,35)
 35           FORMAT(5H EIGA,17X,5H EIGB,17X,20H COMPLEX EIGENVALUES)
              DO 50 I=1,N
                 IF (EIGB(I).NE.0.0) GO TO 45
                 WRITE(IWRITE,43)EIGA(I),EIGB(I)
 43              FORMAT(4E11.3)
                 GO TO 50
 45              CONTINUE
                 EI=EIGA(I)/EIGB(I)
                 WRITE(IWRITE,44)EIGA(I),EIGB(I),EI
 44              FORMAT(4E10.3,2E15.6)
 50           CONTINUE
 60         CONTINUE
C DETERMINE RELATIVE RESIDUAL
            RESERR=0.0
            DO 90 K=1,N
               RESUM=0.0
               DO 80 I=1,N
                  SUMA=0.0
                  SUMB=0.0
                  DO 70 J=1,N
                     SUMA=SUMA+ASAVE(I,J)*X(J,K)
                     SUMB=SUMB+BSAVE(I,J)*X(J,K)
 70               CONTINUE
                  RESUM=RESUM+CABS(SUMA*EIGB(K)-SUMB*EIGA(K))
 80           CONTINUE
              RESERR=AMAX1(RESUM,RESERR)
 90        CONTINUE
           RESERR=RESERR/(ABNORM)
           WRITE(IWRITE,95)N,NA,NB,RESERR,ERR
 95         FORMAT(36H N,RANKA,RANKB,REL.RESID.,ABS. ERROR,3I3,2E15.5)
 100       CONTINUE
           STOP
           END
            SUBROUTINE SETUP(NM,A,B,MR,MI,LR,LI,N,NA,NB,ICASE)
            COMPLEX A(NM,NM),B(NM,NM)
            INTEGER MR(NM,NM),MI(NM,NM),LR(NM),LI(NM)
            N=((ICASE+2)/3)*5
            NA=N
            NB=N
            IF (ICASE.NE.1.AND.ICASE.NE.4.AND.ICASE.NE.7)NB=N-2
            IF (ICASE.EQ.3.OR.ICASE.EQ.6.OR.ICASE.EQ.9)NA=N-1
            DO 20 I=1,N
               DO 10 J=1,N
                   A(I,J)=(0.0,0.0)
                   B(I,J)=(0.0,0.0)
                   MR(I,J)=N*UNI(0)
                   MI(I,J)=2*N*UNI(0)-N
 10            CONTINUE
 20         CONTINUE
            DO 80 I=1,N
               DO 30 K=1,N
                  LR(K)=2*N*UNI(0)-N
                  LI(K)=2*N*UNI(0)-N
 30            CONTINUE
               DO 70 J=1,N
                  DO 50 K=1,NA
                   LLR=LR(K)-LI(K)
                   LLI=LR(K)+LI(K)
                   A(I,J)=A(I,J)+CMPLX((LLR*MR(K,J)-LLI*MI(K,J))*K,(LLR*
     1   MI(K,J)+LLI*MR(K,J))*K)
 50              CONTINUE
                DO 60 K=1,NB
                   B(I,J)=B(I,J)+CMPLX(LR(K)*MR(K,J)-LI(K)*MI(K,J),
     1    LI(K)*MR(K,J)+LR(K)*MI(K,J))
 60               CONTINUE
 70            CONTINUE
 80         CONTINUE
            RETURN
            END
