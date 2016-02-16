C
C     Returns the general complete elliptic integral
C     cel(kc, p, a, b) with QQC = kc, PP = p
C
C     This is the Fortran code from the first edition
C     of "Numerical Recipies"
C
C     The Fortran code is a translation of the ALGOL
C     code from:
C     "Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions. III*" by Roland Bulirsch.
C     Numerische Mathematik 13, 305--315 (1969)
C
C     Complete elliptic integral of the 1st kind:
C      cel(kc, 1, 1, 1)
C     Complete elliptic integral of the 2nd kind:
C      cel(kc, 1, 1, kc*kc)
C     Complete elliptic integral of the 3rd kind:
C      cel(kc, p, 1, 1)
C
C
      FUNCTION CEL(QQC,PP,AA,BB)
      REAL*8 CEL,QQC,PP,AA,BB,A,B,P,E,EM,QC,F,G,Q
      PARAMETER (CA=.0003, PI02=1.5707963268)
      IF(QQC.EQ.0.)PAUSE 'failure in CEL'
      QC=ABS(QQC)
      A=AA
      B=BB
      P=PP
      E=QC
      EM=1.
      IF(P.GT.0.)THEN
              P=SQRT(P)
              B=B/P
      ELSE
              F=QC*QC
              Q=1.-F
              G=1.-P
              F=F-P
              Q=Q*(B-A*P)
              P=SQRT(F/G)
              A=(A-B)/G
              B=-Q/(G*G*P)+A*P
      ENDIF
1     F=A
      A=A+B/P
      G=E/P
      B=B+F*G
      B=B+B
      P=G+P
      G=EM
      EM=QC+EM
      IF(ABS(G-QC).GT.G*CA)THEN
              QC=SQRT(E)
              QC=QC+QC
              E=QC*EM
              GO TO 1
      ENDIF
      CEL=PI02*(B+A*EM)/(EM*(EM+P))
      RETURN
      END
