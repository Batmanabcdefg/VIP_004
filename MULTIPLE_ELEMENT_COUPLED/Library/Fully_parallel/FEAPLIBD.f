
      SUBROUTINE ELMT101(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
    1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT02(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT103(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT04(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT05(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT06(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT07(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT08(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE ELMT09(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1                  NDF,NDM,NST,NS,NQ,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DUMMY ELEMENT ROUTINE

      RETURN
      END
      SUBROUTINE SHAP3(RR,SS,TT,X,SHP,XSJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- 3D SHAPE FUNCTION ROUTINE

      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL
      DIMENSION SHP(4,*),X(3,*),R(8),S(8),T(8),XS(3,3),
     1 SX(3,3),DUM(3)
      DATA R/-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5/
      DATA S/-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5/
      DATA T/-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5/

!---- FORM 8-NODE BRICK SHAPE FUNCTIONS
      DO 10 I = 1,8
      SHPR = 0.5 + R(I)*RR
      SHPS = 0.5 + S(I)*SS
      SHPT = 0.5 + T(I)*TT
      SHP(1,I) = R(I)*SHPS*SHPT
      SHP(2,I) = SHPR*S(I)*SHPT
      SHP(3,I) = SHPR*SHPS*T(I)
   10 SHP(4,I) = SHPR*SHPS*SHPT

!---- CONSTRUCT JACOBIAN AND INVERSE
      DO 20 I = 1,3
      DO 20 J = 1,3
      XS(I,J) = 0.D00
      DO 20 K = 1,8
   20 XS(I,J) = XS(I,J) + X(I,K)*SHP(J,K)

      XSJ = XS(1,1)*XS(2,2)*XS(3,3) + XS(2,1)*XS(3,2)*XS(1,3)
     1    + XS(1,2)*XS(2,3)*XS(3,1) - XS(1,1)*XS(2,3)*XS(3,2)
     2    - XS(1,2)*XS(2,1)*XS(3,3) - XS(1,3)*XS(2,2)*XS(3,1)

      IF(XSJ.LE.0.D0) THEN
      WRITE(6,2001) NELMT,XSJ
      STOP
      END IF
      SX(1,1) = (XS(2,2)*XS(3,3) - XS(3,2)*XS(2,3))/XSJ
      SX(1,2) = (XS(3,2)*XS(1,3) - XS(1,2)*XS(3,3))/XSJ
      SX(1,3) = (XS(1,2)*XS(2,3) - XS(2,2)*XS(1,3))/XSJ
      SX(2,1) = (XS(3,1)*XS(2,3) - XS(2,1)*XS(3,3))/XSJ
      SX(2,2) = (XS(1,1)*XS(3,3) - XS(1,3)*XS(3,1))/XSJ
      SX(2,3) = (XS(2,1)*XS(1,3) - XS(1,1)*XS(2,3))/XSJ
      SX(3,1) = (XS(2,1)*XS(3,2) - XS(3,1)*XS(2,2))/XSJ
      SX(3,2) = (XS(3,1)*XS(1,2) - XS(1,1)*XS(3,2))/XSJ
      SX(3,3) = (XS(1,1)*XS(2,2) - XS(1,2)*XS(2,1))/XSJ

!---- FORM GLOBAL DERIVATIVES
      DO 30 I = 1,8
      DO 40 J = 1,3
      DUM(J) = 0.D0
      DO 40 K = 1,3
   40 DUM(J) = DUM(J) + SHP(K,I)*SX(K,J)
      DO 30 J = 1,3
   30 SHP(J,I) = DUM(J)
      RETURN
2001  FORMAT(//5X,'** FATAL ERROR ** JACOBIAN OF TRANSFORMATION',
     1 ' FOR ELEM',I5,' IS ',E15.7)
      END
      SUBROUTINE PGAUS3(L,LINT,RR,SS,TT,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- 3D GAUSS POINTS AND WEIGHTS

      DIMENSION RR(*),SS(*),TT(*),W(*)
      DIMENSION R(8),S(8),T(8)
      DATA R/-1.,1.,1.,-1.,-1.,1.,1.,-1./,
     1     S/-1.,-1.,1.,1.,-1.,-1.,1.,1./,
     2     T/-1.,-1.,-1.,-1.,1.,1.,1.,1./
      LL = L
      IF (LL.LE.0) LL = 1
      IF (LL.GT.2) LL = 2
      LINT = LL*LL*LL
      GO TO (10,20),LL

!---- 1ST ORDER QUADRATURE
   10 RR(1) = 0.D0
      SS(1) = 0.D0
      TT(1) = 0.D0
      W(1)  = 8.D0
      RETURN

!---- 2ND ORDER QUADRATURE
   20 CONST = 1./DSQRT(3.D0)
!  20 CONST = 1./SQRT(3.E0)
      DO 30 I = 1,8
      RR(I) = R(I)*CONST
      SS(I) = S(I)*CONST
      TT(I) = T(I)*CONST
   30 W(I)  = 1.D0
      RETURN

      END
      SUBROUTINE PGAUS2(L,LINT,R,Z,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- GAUSS POINTS AND WEIGHTS FOR TWO DIMENSIONS

      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
      DIMENSION LR(9),LZ(9),LW(9),R(*),Z(*),W(*),G4(4),H4(4)
      DATA LR/-1,1,1,-1,0,1,0,-1,0/,LZ/-1,-1,1,1,-1,0,1,0,0/
      DATA LW/4*25,4*40,64/
      IF(L.LT.0) GO TO 10
      LINT = L*L
      GO TO (1,2,3,4),L
!---- 1X1 INTEGRATION
1     R(1) = 0.
      Z(1) = 0.
      IF(NEL.EQ.3) Z(1) = -1./3.
      W(1) = 4.
      RETURN
!---- 2X2 INTEGRATION
2     G = 1./DSQRT(3.0D0)
!2     G = 1./SQRT(3.0E0)
      DO 21 I = 1,4
      R(I) = G*LR(I)
      Z(I) = G*LZ(I)
21    W(I) = 1.
      RETURN
!---- 3X3 INTEGRATION
3     G = DSQRT(0.6D0)
!3     G = SQRT(0.6E0)
      H = 1./81.
      DO 31 I = 1,9
      R(I) = G*LR(I)
      Z(I) = G*LZ(I)
31    W(I) = H*LW(I)
      RETURN
!---- 4X4 INTEGRATION
   4  G = DSQRT(4.8D0)
!4     G = SQRT(4.8E0)
      H = DSQRT(30.0D0)/36.
!     H = SQRT(30.0E0)/36.
      G4(1) = DSQRT((3.0D0+G)/7.)
!     G4(1) = SQRT((3.0E0+G)/7.)
      G4(4) = - G4(1)
      G4(2) = DSQRT((3.0D0-G)/7.)
!     G4(2) = SQRT((3.0E0-G)/7.)
      G4(3) = -G4(2)
      H4(1) = 0.5 - H
      H4(2) = 0.5 + H
      H4(3) = 0.5 + H
      H4(4) = 0.5 - H
      I = 0
      DO 41 J = 1,4
      DO 41 K = 1,4
      I = I + 1
      R(I) = G4(K)
      Z(I) = G4(J)
      W(I) = H4(J)*H4(K)
41    CONTINUE
      RETURN
10    LINT = 3
      G = DSQRT(1.0D0/3.0D0)
!     G = SQRT(1.0E0/3.0E0)
      H = DSQRT(2.0D0/3.0D0)
!     H = SQRT(2.0E0/3.0E0)
      R(1) = -H
      R(2) =  H
      R(3) = 0.
      Z(1) = -G
      Z(2) = -G
      Z(3) = G
      W(1) = 1.
      W(2) = 1.
      W(3) = 2.
      RETURN
      END
      SUBROUTINE SHAP2D(SS,TT,X,SHP,XSJ,NDM,NEL,IX,FLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FLAG

!---- SHAPE FUNCTION ROUTINE FOR TWO DIMENSIONAL ELEMENTS

      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL1
      DIMENSION SHP(3,*),X(NDM,*),S(4),T(4),XS(2,2),SX(2,2),IX(*)
      DATA S/-0.5,0.5,0.5,-0.5/,T/-0.5,-0.5,0.5,0.5/
!---- FORM 4-NODE QUADRILATERAL SHAPE FUNCTIONS
      DO 100 I = 1,4
      SHP(3,I) = (0.5+S(I)*SS)*(0.5+T(I)*TT)
      SHP(1,I) = S(I)*(0.5+T(I)*TT)
100   SHP(2,I) = T(I)*(0.5+S(I)*SS)
      IF(NEL.GE.4) GO TO 120
!---- FORM TRIANGLE BY ADDING THIRD AND FOURTH TOGETHER
      DO 110 I = 1,3
110   SHP(I,3) = SHP(I,3)+SHP(I,4)
!---- ADD QUADRATIC TERMS IF NECESSARY
120   IF(NEL.GT.4) CALL QUAD(SS,TT,SHP,IX,NEL)
!---- CONSTRUCT JACOBIAN AND ITS INVERSE
      DO 130 I = 1,NDM
      DO 130 J = 1,2
      XS(I,J) = 0.D0
      DO 130 K = 1,NEL
130   XS(I,J) = XS(I,J) + X(I,K)*SHP(J,K)
      XSJ = XS(1,1)*XS(2,2)-XS(1,2)*XS(2,1)
      IF(XSJ.LE.0.D0) THEN
      WRITE(6,2001) NELMT,XSJ
      STOP
      END IF
      IF (.NOT.FLAG) RETURN
      SX(1,1) = XS(2,2)/XSJ
      SX(2,2) = XS(1,1)/XSJ
      SX(1,2) =-XS(1,2)/XSJ
      SX(2,1) =-XS(2,1)/XSJ
!---- FORM GLOBAL DERIVATIVES
      DO 140 I = 1,NEL
      TP        = SHP(1,I)*SX(1,1)+SHP(2,I)*SX(2,1)
      SHP(2,I)  = SHP(1,I)*SX(1,2)+SHP(2,I)*SX(2,2)
140   SHP(1,I) = TP
      RETURN
2001  FORMAT(//5X,'** FATAL ERROR ** JACOBIAN OF TRANSFORMATION',
     1 ' FOR ELEM',I5,' IS ',E15.7)
      END
      SUBROUTINE QUAD(S,T,SHP,IX,NEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ADD QUADRATIC FUNCTIONS AS NECESSARY

      DIMENSION IX(*),SHP(3,*)
      S2 = (1.-S*S)/2.
      T2 = (1.-T*T)/2.
      DO 100 I=5,9
      DO 100 J = 1,3
100   SHP(J,I) = 0.D0
!---- MIDSIDE NODES (SERENDIPITY)
      IF(IX(5).EQ.0) GO TO 101
      SHP(1,5) = -S*(1.-T)
      SHP(2,5) = -S2
      SHP(3,5) = S2*(1.-T)
101   IF(NEL.LT.6) GO TO 107
      IF(IX(6).EQ.0) GO TO 102
      SHP(1,6) = T2
      SHP(2,6) = -T*(1.+S)
      SHP(3,6) = T2*(1.+S)
102   IF(NEL.LT.7) GO TO 107
      IF(IX(7).EQ.0) GO TO 103
      SHP(1,7) = -S*(1.+T)
      SHP(2,7) = S2
      SHP(3,7) = S2*(1.+T)
103   IF(NEL.LT.8) GO TO 107
      IF(IX(8).EQ.0) GO TO 104
      SHP(1,8) = -T2
      SHP(2,8) = -T*(1.-S)
      SHP(3,8) = T2*(1.-S)
!---- INTERIOR NODE (LAGRANGIAN)
104   IF(NEL.LT.9) GO TO 107
      IF(IX(9).EQ.0) GO TO 107
      SHP(1,9) = -4.*S*T2
      SHP(2,9) = -4.*T*S2
      SHP(3,9) = 4.*S2*T2
!---- CORRECT EDGE NODES FOR INTERIOR NODE (LAGRANGIAN)
      DO 106 J= 1,3
      DO 105 I = 1,4
105   SHP(J,I) = SHP(J,I) - 0.25*SHP(J,9)
      DO 106 I = 5,8
106   IF(IX(I).NE.0) SHP(J,I) = SHP(J,I) - .5*SHP(J,9)
!---- CORRECT CORNER NODES FOR PRESENSE OF MIDSIDE NODES
107   K = 8
      DO 109 I = 1,4
      L = I + 4
      DO 108 J = 1,3
108   SHP(J,I) = SHP(J,I) - 0.5*(SHP(J,K)+SHP(J,L))
109   K = L
      RETURN
      END

      SUBROUTINE PGAUS1(LINT,R,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- GAUSS POINTS AND WEIGHTS FOR 1 DIMENSION (BEAM ELEMENT)

      DIMENSION R(*),W(*)
      GO TO (1,2,3) LINT
!---- 1 POINT QUAD
1     CONTINUE
      R(1)=0.D0
      W(1)=2.D0
      RETURN
!---- 2 POINT QUAD
2     CONTINUE
      G=1.D0/DSQRT(3.0D0)
!     G=1.0/SQRT(3.0E0)
      R(1)=-G
      R(2)=G
      W(1)=1.D0
      W(2)=1.D0
      RETURN
!---- 3 POINT QUAD
3     CONTINUE
      G=DSQRT(0.6D0)
!     G=SQRT(0.6E0)
      H1=5.D0/9.D0
      H2=8.D0/9.D0
      R(1)=-G
      R(2)=0.D0
      R(3)=G
      W(1)=H1
      W(2)=H2
      W(3)=H1
      RETURN
      END

      SUBROUTINE SUBST(W,B,X,IPIVOT,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION W(N,*),B(*),X(*),IPIVOT(*)

      IF (N.GT.1) GO TO 10
      X(1) = B(1)/W(1,1)
      RETURN
   10 IP = IPIVOT(1)
      X(1) = B(IP)
      DO 15 K = 2,N
      IP = IPIVOT(K)
      KM1 = K - 1
      SUM = 0.D0
      DO 14 J = 1,KM1
   14 SUM = W(IP,J)*X(J) + SUM
   15 X(K) = B(IP) - SUM
      X(N) = X(N)/W(IP,N)
      K = N
      DO 20 NP1MK = 2,N
      KP1 = K
      K = K - 1
      IP = IPIVOT(K)
      SUM = 0.D0
      DO 19 J = KP1,N
   19 SUM = W(IP,J)*X(J) + SUM
   20 X(K) = (X(K) - SUM)/W(IP,K)
      RETURN
      END

      SUBROUTINE FACTR(W,IPIVOT,D,N,FLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL FLAG

      DIMENSION W(N,*),IPIVOT(*),D(*)

      DATA TOL/1.D-50/
!     DATA TOL/1.E-12/
      FLAG = .TRUE.
!---- INITIALIZE W,IPIVOT,D
      DO 10 I = 1,N
      IPIVOT(I) = I
      ROWMAX = 0.D0
      DO 9 J = 1,N
    9 ROWMAX = DMAX1(ROWMAX,DABS(W(I,J)))
!   9 ROWMAX = AMAX1(ROWMAX,ABS(W(I,J)))
!     IF (ROWMAX.EQ.0.0) GO TO 999
      IF (DABS(ROWMAX).LT.TOL) GO TO 999
!     IF (ABS(ROWMAX).LT.TOL) GO TO 999
   10 D(I) = ROWMAX
!---- GAUSS ELIMINATION WITH SCALED PARTIAL PIVOTING
      NM1 = N - 1
      IF (NM1.EQ.0) RETURN
      DO 20 K = 1,NM1
      J = K
      KP1 = K + 1
      IP = IPIVOT(K)
      COLMAX = DABS(W(IP,K))/D(IP)
!     COLMAX = ABS(W(IP,K))/D(IP)
      DO 11 I = KP1,N
      IP = IPIVOT(I)
      AWIKOV = DABS(W(IP,K))/D(IP)
!     AWIKOV = ABS(W(IP,K))/D(IP)
      IF (AWIKOV.LE.COLMAX) GO TO 11
      COLMAX = AWIKOV
      J = I
   11 CONTINUE
!     IF (COLMAX.EQ.0.0) GO TO 999
      IF (DABS(COLMAX).LT.TOL) GO TO 999
!     IF (ABS(COLMAX).LT.TOL) GO TO 999
      IPK = IPIVOT(J)
      IPIVOT(J) = IPIVOT(K)
      IPIVOT(K) = IPK
      DO 20 I = KP1,N
      IP = IPIVOT(I)
      W(IP,K) = W(IP,K)/W(IPK,K)
      RATIO = - W(IP,K)
      DO 20 J = KP1,N
   20 W(IP,J) = RATIO*W(IPK,J) + W(IP,J)
!     IF (W(IP,N).EQ.0.0) GO TO 999
      IF (DABS(W(IP,N)).LT.TOL) GO TO 999
!     IF (ABS(W(IP,N)).LT.TOL) GO TO 999
      RETURN
  999 FLAG = .FALSE.
      RETURN
      END
 
      SUBROUTINE SOLVE(A,B,N,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      LOGICAL FLAG
 
      DIMENSION A(N,*),B(N,*)
      DIMENSION D(50),IPIVOT(50),X(50)

!---- SOLVE A*X = B AND RETURN X IN B

      IF ((N.LE.0).OR.(N.GT.50)) THEN
      WRITE(6,'(///'' **ERROR** DETECTED BY SUBROUTINE SOLVE''
     1 /'' SYSTEM DIMENSION '',I5,'' OUT OF RANGE'')') N
      STOP
      END IF

!---- FACTORIZE MATRIX
      CALL FACTR(A,IPIVOT,D,N,FLAG)
      IF (.NOT.FLAG) THEN
      WRITE(6,'(///'' **ERROR** DETECTED BY SUBROUTINE SOLVE''
     1 /'' SINGULAR SYSTEM'')')
      STOP
      END IF

!---- BACK-SUBSTITUTE
      DO 20 IB = 1,NB
      CALL SUBST(A,B(1,IB),X,IPIVOT,N)
      DO 10 I = 1,N
   10 B(I,IB) = X(I)
   20 CONTINUE
      RETURN
      END

     SUBROUTINE MA01(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
    1 STR,EPS,Q,B,DR,VELG,ACCELG,VEL,ACCEL,CT,NDF,NDM,NEN1,
    2 NST,NSTR,NQ,NEND,FLG)
     IMPLICIT DOUBLE PRECISION (A-H,O-Z)

     LOGICAL FLG

     RETURN
     END

      SUBROUTINE MA02(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
     1 STR,EPS,Q,B,DR,VELG,ACCELG,VEL,ACCEL,CT,NDF,NDM,NEN1,
     2 NST,NSTR,NQ,NEND,FLG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      LOGICAL FLG

      RETURN
      END

      SUBROUTINE MA03(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
     1 STR,EPS,Q,B,DR,VELG,ACCELG,VEL,ACCEL,CT,NDF,NDM,NEN1,
     2 NST,NSTR,NQ,NEND,FLG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FLG

      RETURN
      END

      SUBROUTINE SUBSP(A,B,V,T,G,H,D,DP,DTOL,P,Z,JDIAG,NF,NV,NEQ,IMAS
     1                 ,TOL,SHIFT,PRT,ITS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION A(*),B(*),V(NEQ,*),T(*),G(*),H(*),D(*),DP(*),DTOL(*),
     1          P(NV,*),Z(NEQ,*),JDIAG(*)

      RETURN
      END

      SUBROUTINE PLOT1(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,STR,EPS,Q,
     1 B,DR,VELG,ACCELG,VEL,ACCEL,NDF,NDM,NEN1,NST,NSTR,NQ,NEND,KEY,IQ,
     2 PRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JDIAG(*),UL(*),XL(*),TL(*),LD(*),P(*),S(*),IE(*),D(*),
     2 IX(NEN1,*),F(*),T(*),B(*),DR(*),ID(*),X(NDM,*),IX1(8)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(*),ACCELG(*),VEL(*),ACCEL(*)
      DIMENSION SIG(6)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,NSDM,NQDM,NQUAD
!---- WRITE INPUT/OUTPUT FILES FOR G-SAFE PRE-POST PROCESSOR
      IF(KEY.EQ.0) THEN
!     OPEN(UNIT=51,FILE='CDATA.IN1',FORM='FORMATTED',ACCESS=
!    1 'SEQUENTIAL')
!     OPEN(UNIT=52,FILE='CBOUND.IN1',FORM='FORMATTED',ACCESS=
!    1 'SEQUENTIAL')
      OPEN(UNIT=53,FILE='CGEOM.IN1',FORM='FORMATTED',ACCESS=
     1 'SEQUENTIAL')
!     OPEN(UNIT=54,FILE='CMATS.IN1',FORM='FORMATTED',ACCESS=
!    1 'SEQUENTIAL')
!     OPEN(UNIT=55,FILE='ELODIN.IN3',FORM='FORMATTED',ACCESS=
!    1 'SEQUENTIAL')
!     OPEN(UNIT=56,FILE='POST.INF',FORM='FORMATTED',ACCESS=
!    1 'SEQUENTIAL')
!     OPEN(UNIT=57,FILE='STRIN.IN9',FORM='FORMATTED',ACCESS=
!    1 'SEQUENTIAL')
      IDUM = 0
      WRITE(53,5001) NUMEL,NUMNP,IDUM
      DO 500 N=1,NUMEL
      MTP=601
      IF(NEN.EQ.8) MTP=602
      IMAT=1
      IF(NEN.EQ.4) THEN
      DO 490 I = 1,4 
      IX1(I) = IX(I,N)
490   CONTINUE
      ELSE
      II = 1
      DO 491 I=1,4
      IX1(II) = IX(I,N)
      II = II + 2
491   CONTINUE
      II = 2
      DO 492 I=5,8
      IX1(II) = IX(I,N)
      II = II + 2
492   CONTINUE
      END IF
      WRITE(53,5001) N,MTP,IMAT,NEN,(IX1(I),I=1,NEN)
500   CONTINUE
      DO 501 N=1,NUMNP
      WRITE(53,5002) N,X(1,N),X(2,N)
501   CONTINUE
      END IF

      OPEN(UNIT=61,FILE='FDISP.OUT9',FORM='FORMATTED',ACCESS=
     1 'SEQUENTIAL')
      OPEN(UNIT=62,FILE='FREACT.OUT9',FORM='FORMATTED',ACCESS=
     1 'SEQUENTIAL')
      OPEN(UNIT=63,FILE='FSTRS.OUT9',FORM='FORMATTED',ACCESS=
     1 'SEQUENTIAL')
!---- WRITE DISPLACEMENT DATA
      II = 1
      DO 600 I=1,NUMNP
      WRITE(61,6001) I,B(II),B(II+1)
      II = II + 2
600   CONTINUE
      II = 1
      DO 601 I=1,NUMNP
      ID1 = ID(II)
      ID2 = ID(II+1)
      IF((ID1.NE.0).OR.(ID2.NE.0)) THEN
      DUM1 = 1.D0
      DUM2 = 1.D0
      WRITE(62,6002) I,DUM1,DUM2
      END IF
      II = II + 2
601   CONTINUE
      DUM1 = 0.D0
      DO 602 N=1,NUMEL
      WRITE(63,6001) N
      NQP1 = 1
      NQP2 = 4
      IF(NEN.GT.4) THEN
      NQP1 = 10
      NQP2 = 13
      END IF 
      DO 603 IQUAD =NQP1,NQP2    
      II =(IQUAD-1)*NSDM 
      IJ = (IQUAD-1)*NQDM
!     JJ = IQUAD*NQDM
!     WRITE(63,6003)IQUAD,Q(JJ-1,N),Q(JJ,N),STR(II+2,N),STR(II+3,N),
!---- PLOT STRESS CONTOURS
      SIG(5)=0.D0
      SIG(6)=0.D0
!---- EFFECTIVE STRESS
      DO 463 I=1,NSDM
463   SIG(I)=STR(II+I,N)
      SEFFEC = SEF(SIG)
      SIGH = (SIG(1)+SIG(2)+SIG(3))/3.D0
      IF(IQ.EQ.0) THEN
!---- PLOT CONTOURS OF SIGH, SYY, SIGM
      WRITE(63,6003)IQUAD,EPS(II+1,N),EPS(II+2,N),SIGH,
     1 STR(II+3,N),Q(IJ+2,N),DUM1,DUM1,DUM1
      END IF
      IF(IQ.EQ.1) THEN 
!---- PLOT CONTOURS OF SEFF, EMP,F
!---- PLOT CONTOURS OF INTERNAL VARIABLES
      WRITE(63,6003)IQUAD,EPS(II+1,N),EPS(II+2,N),SEFFEC,
     1 Q(IJ+3,N),Q(IJ+4,N),DUM1,DUM1,DUM1
      END IF
      IF(IQ.EQ.2) THEN 
!---- PLOT CONTOURS OF FN, F, SIGH
      WRITE(63,6003)IQUAD,EPS(II+1,N),EPS(II+2,N),Q(IJ+9,N),
     1 Q(IJ+4,N),SIGH,DUM1,DUM1,DUM1
      END IF
      IF(IQ.EQ.3) THEN
!---- PLOT CONTOURS OF SXX, SYY, SXY
      WRITE(63,6003)IQUAD,EPS(II+1,N),EPS(II+2,N),STR(II+2,N)
     1 ,STR(II+3,N),STR(II+4,N),DUM1,DUM1,DUM1
      END IF
603   CONTINUE
602   CONTINUE
      RETURN
6001  FORMAT(I10,3E14.6)
6002  FORMAT(I10,6E14.4)
6003  FORMAT(I10,8E14.6)
5001  FORMAT(15I5)
5002  FORMAT(I10,7E14.7)
      END
      SUBROUTINE PLOT2(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,STR,EPS,Q,
     1 B,DR,VELG,ACCELG,VEL,ACCEL,NDF,NDM,NEN1,NST,NSTR,NQ,NEND,KEY,IQ,
     2 PRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION JDIAG(*),UL(*),XL(*),TL(*),LD(*),P(*),S(*),IE(*),D(*),
     2 IX(NEN1,*),F(*),T(*),B(*),DR(*),ID(*),X(NDM,*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(*),ACCELG(*),VEL(*),ACCEL(*)
      DIMENSION SDEV(4)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,NSDM,NQDM,NQUAD
      DATA LPLOT/43/
      REWIND(LPLOT)
      WRITE(LPLOT) NEQ,NUMEL,NUMNP,NEN
!---- WRITE NODAL COORDS AND NODAL DISPLACEMENTS
      WRITE(LPLOT) (X(1,J),J=1,NUMNP),(X(2,J),J=1,NUMNP),
     1 (B(I),I=1,NEQ)
!---- WRITE NODAL CONNECTIVITY ARRAYS FOR ALL ELEMENTS
      WRITE(LPLOT) ((IX(I,J),I=1,NEN),J=1,NUMEL)
!---- INITALIZE
      DO 10 I=1,NEQ
      DR(I) = 0.D0
      ACCELG(I) = 0.D0
10    CONTINUE
      NQD1 = (NQUAD-1)*NSDM
!---- WRITE NODAL STRESSES, MISES EQV STRESS AND HYDROSTATIC STRESS
!---- NODAL STRESSES ARE OBTAINED BY AVERAGING OVER ALL ELEMENTS MEETING
!---- AT A NODE
      NSD2 = NSDM + 2
      DO 20 I=1,NSD2
      DO 30 N=1,NUMEL
      IF(I.LE.NSDM) THEN
      SIG = STR(NQD1+I,N)
      END IF
      IF(I.EQ.(NSDM+1)) THEN
      SIG = 0.D0
      PRES = (STR(NQD1+1,N)+STR(NQD1+2,N)+STR(NQD1+3,N))/3.D0
      DO 32 K=1,3
      SDEV(K) = STR(NQD1+K,N)-PRES
32    CONTINUE
      SDEV(4) = STR(NQD1+4,N)
      DO 33 K=1,3
      SIG = SIG + SDEV(K)*SDEV(K)
33    CONTINUE
      SIG = SIG + SDEV(4)*SDEV(4)*2.D0
      SIG = DSQRT(3.D0*SIG/2.D0)
!     SIG = SQRT(3.0*SIG/2.0)
      END IF
      IF(I.EQ.NSD2) THEN
      SIG = (STR(NQD1+1,N)+STR(NQD1+2,N)+STR(NQD1+3,N))/3.D0
      END IF
      DO 35 J=1,NEN
      II = IX(J,N)
!---- II IS GLOBAL NODE # CORR TO LOCAL NODE J
!---- ADD CONTRIBUTION TO GLOBAL NODE II
      DR(II) = DR(II) + SIG
!---- INCREMENT COUNTER FOR # CONTRIBUTIONS FROM ELEMS SURROUNDING
!---- GLOBAL NODE II
      ACCELG(II) = ACCELG(II) + 1.D0
35    CONTINUE
30    CONTINUE
      DO 40 II=1,NUMNP
      IF(ACCELG(II).NE.0.D0) THEN
      DR(II) = DR(II)/ACCELG(II)
      ELSE
      DR(II) = 0.D0
      END IF
40    CONTINUE
      WRITE(LPLOT) (DR(II),II=1,NUMNP)
      DO 45 II=1,NUMNP
      DR(II) = 0.D0
      ACCELG(II) = 0.D0
45    CONTINUE
20    CONTINUE
!---- WRITE NODAL STRAINS. NODAL STRAINS ARE OBTAINED BY ELEM AVERAGING
      DO 120 I=1,NSDM
      DO 130 N=1,NUMEL
      SIG=EPS(NQD1+I,N)
      DO 135 J=1,NEN
      II = IX(J,N)
      DR(II) = DR(II) + SIG
      ACCELG(II) = ACCELG(II) + 1.D0
135   CONTINUE
130   CONTINUE
      DO 140 II=1,NUMNP
      IF(ACCELG(II).NE.0.D0) THEN
      DR(II) = DR(II)/ACCELG(II)
      ELSE
      DR(II) = 0.D0
      END IF
140   CONTINUE
      WRITE(LPLOT) (DR(II),II=1,NUMNP)
      DO 145 II=1,NUMNP
      DR(II) = 0.D0
      ACCELG(II) = 0.D0
145   CONTINUE
120   CONTINUE
!---- WRITE INTERNAL VARIABLES 1 TO 3
      NQD2 = (NQUAD-1)*NQDM
      DO 220 I=1,7
      DO 230 N=1,NUMEL
      SIG = Q(NQD2+I,N)
      DO 235 J=1,NEN
      II = IX(J,N)
      DR(II) = DR(II) + SIG
      ACCELG(II) = ACCELG(II) + 1.D0
235   CONTINUE
230   CONTINUE
      DO 240 II=1,NUMNP
      IF(ACCELG(II).NE.0.D0) THEN
      DR(II) = DR(II)/ACCELG(II)
      ELSE
      DR(II) = 0.D0
      END IF
240   CONTINUE
      WRITE(LPLOT) (DR(II),II=1,NUMNP)
      DO 245 II=1,NUMNP
      DR(II) = 0.D0
      ACCELG(II) = 0.D0
245   CONTINUE
220   CONTINUE
      RETURN
      END

     SUBROUTINE MP101(IDL,IE,D,ID,X,IX,F,T,DR,STR,EPS,Q,VEL,ACCEL,
    1 NDF,NDM,NEN1,NSTR,NQ,III,PRTN)
     IMPLICIT DOUBLE PRECISION (A-H,O-Z)

     LOGICAL PRTN

     RETURN
     END
      SUBROUTINE MP02(IDL,IE,D,ID,X,IX,F,T,DR,STR,EPS,Q,VEL,ACCEL,
     1 NDF,NDM,NEN1,NSTR,NQ,III,PRTN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL PRTN

      RETURN
      END
      SUBROUTINE DJ(D,UL,XL,STRL,EPSL,QL,IX,TL,VEL,ACCEL,
     1 NDF,NDM,NST,NS,NQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      RETURN
      END


      SUBROUTINE JINTT(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,
     1 JDIAG,STR,EPS,Q,B,DR,VELG,ACCELG,VEL,ACCEL,CT,NDF,NDM,NEN1,
     2 NST,NSTR,NQ,NEND,FLG)
!---- DYNAMIC ENERGY RELEASE RATE (DOMAIN INTEGRAL STATIONARY CRACK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL FLG,DUMPJ
      LOGICAL AFR,BFR,CFR,AFL,BFL,CFL,DFL,EFL,GFL,TFL,IFL
     1 ,BLKFL,SRCLIN
      CHARACTER*8 HEAD
      CHARACTER*1 O
!     COMMON M(1)
      COMMON /MAINN/ M(25000000)
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,
     1 IPR,NSDM,NQDM,NQUAD
!     COMMON /PRLOD/ PROP,PROPOL
!     COMMON /PTABLE/ NPLD,KPLD,T1,P1,T2,P2
      COMMON /SUBDT/ MD,MV,MF
      COMMON /TDATA/ TIME,DT,BETA,GAMMA
      COMMON /POINTR/ NE,NN,NA,NC,NV,NM,NA2,NA3
      COMMON /FLAGS/ AFR,BFR,CFR,AFL,BFL,CFL,DFL,EFL,GFL,TFL,IFL
     1 ,BLKFL,SRCLIN
      COMMON /JINT/ NDOM(3),XC(8,4),DOM(4,10),
     1 DYNJ(8,7,10),DYNJM(8,7,10),TSWE,TKEE,TKEE1
      DIMENSION CT(3),JDIAG(*),UL(*),XL(*),TL(*),LD(*),P(*),
     2 S(*),IE(7,*),D(*),ID(*),X(NDM,*),IX(NEN1,*),F(*),T(*),
     3 B(*),DR(*),PROP(7),IDPROP(*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(*),ACCELG(*),VEL(*),ACCEL(*)
!     DIMENSION XRM(1)
!     DIMENSION XRM(10000)
!     EQUIVALENCE (M(1),XRM(1))
      DATA JDUMP/10/

      XKEY0 =4H    
      XKEY1 =4HAVRG
      DO 410 NDI = 1,8
      DO 410 I = 1,7
      DO 410 J = 1,10
      DYNJM(NDI,I,J) = 0.D0
  410 DYNJ(NDI,I,J) = 0.D0
      TSWE = 0.D0
      TKEE = 0.D0
      TKEE1 = 0.D0
      CALL PFORM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,JDIAG,
     1 STR,EPS,Q,DR,DR,DR,VELG,ACCELG,VEL,ACCEL,
     2 NDF,NDM,NEN1,NST,NSTR,NQ,11,B,
     3 M(NV),.FALSE.,.FALSE.,.FALSE.,.FALSE.)
!---- OUTPUT
!     INQUIRE(JDUMP,OPENED=DUMPJ)
!     IF(.NOT.DUMPJ) THEN
!      OPEN(JDUMP,FILE='FPJINT',STATUS='NEW',FORM='FORMATTED')
!     END IF

      DO 418 NDI = 1,NDOM(3)
      DO 418 K = 1,NDOM(1)
!---- DOUBLE FOR SYMMETRY
      IF(NDOM(2).NE.0) THEN
       DO 414 I = 1,7
!      DYNJM(NDI,I,K) = DYNJM(NDI,I,K)*2.0
 414   DYNJ(NDI,I,K) = DYNJ(NDI,I,K)*2.D0
      ENDIF
!---- DIVIDE BY VCE AREA (FOR THREE-DIMENSIONAL CRACK)
      DO 416 I = 1,7
!     DYNJM(NDI,I,K) = DYNJM(NDI,I,K)/XC(NDI,4)
  416  DYNJ(NDI,I,K) = DYNJ(NDI,I,K)/XC(NDI,4)
      IF(K.NE.1) DYNJ(NDI,6,1) = DYNJ(NDI,6,1) + DYNJ(NDI,6,K)
      IF(K.NE.1) DYNJ(NDI,7,1) = DYNJ(NDI,7,1) + DYNJ(NDI,7,K)
!     IF(K.NE.1) DYNJM(NDI,4,1) = DYNJM(NDI,4,1) + DYNJM(NDI,4,K)
      IF(CT(1).NE.XKEY1)
     1  WRITE(JDUMP,'(2I3,7(1PG12.4))') NDI,K,(DYNJ(NDI,J,K),J=1,7)
!     IF(CT(1).NE.XKEY1)
!    1  WRITE(JDUMP,'(2I3,5(1PG12.4))') NDI,K,(DYNJM(NDI,J,K),J=1,4)
 418  CONTINUE
!     WRITE(JDUMP,'(7(1PE12.4)/7(1PE11.3))') TIME,
!    1  (DYNJ(NDI,6,1)/FLOAT(NDOM(1)),NDI=1,NDOM(3)),
!    2  (DYNJ(NDI,7,1)/FLOAT(NDOM(1)),NDI=1,NDOM(3)),TSWE,TKEE,TKEE1
      WRITE(JDUMP,'(7(1PE12.4)/7(1PE11.3))') TIME,
     1  (DYNJ(NDI,6,1)/DFLOAT(NDOM(1)),NDI=1,NDOM(3)),
     2  (DYNJ(NDI,7,1)/DFLOAT(NDOM(1)),NDI=1,NDOM(3)),TSWE,TKEE,TKEE1
!     WRITE(JDUMP,'(F9.3,6(1PE11.3)/7(1PE11.3))') TIME,
!    1  (DYNJM(NDI,4,1)/FLOAT(NDOM(1)),NDI=1,NDOM(3)),TSWE,TKEE
      RETURN
      END

      SUBROUTINE QWT(XL,QW,NDM,NEL,N,LDOM)

!---- 2 - D ENERGY RELEASE RATE FOR MODE I STATIONARY CRACK

!---- COMPUTE WEIGHTS QW

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
      COMMON/JINT/ NDOM(3),XC(8,4),DOM(4,10),DYNJ(8,7,10),
     1 DYNJM(8,7,10),TSWE,TKEE,TKEE1

      DIMENSION XL(NDM,*),QW(*)

      XIFN = DOM(4,N)
      IF(XIFN.GT.1.D0) GOTO 2
!---- BILINEAR FUNCTION (RECTANGULAR DOMAIN)
      XX = DOM(1,N)
      YY = DOM(2,N)
      IF(DOM(3,N).EQ.0.D0) THEN
       XI = 0.D0
       YI = 0.D0
      ELSE
       XI = 1.D0 - DOM(3,N)/XX
       YI = 1.D0 - DOM(3,N)/YY
      ENDIF
      DO 10 I = 1,NEL
      QW(I) = 1.D0
      X = XL(1,I)
      Y = XL(2,I)
      A  = (X - XC(1,1))/XX
      B  = (Y - XC(1,2))/YY
      A = DABS(A)
!     A = ABS(A)
      B = DABS(B)
!     B = ABS(B)
      IF(A.LT.1.D0.AND.B.LT.1.D0) THEN
       IF(A.GT.XI) QW(I) = (1.D0-A)/(1.D0-XI)
       IF(B.GT.YI) QW(I) = QW(I)*(1.D0-B)/(1.D0-YI)
       LDOM = 1
      ELSE
       QW(I) = 0.D0
      END IF
 10   CONTINUE
      RETURN
!---- CONICAL FUNCTION (CIRCULAR DOMAIN)
    2 CONTINUE
      RR = DOM(1,N)
      IF (DOM(3,N).EQ.0.D0) THEN
       RI = 0.D0
      ELSE
       RI = 1.D0 - DOM(3,N)/RR
      ENDIF
      DO 20 I=1,NEL
      QW(I) = 1.D0
      X = XL(1,I)
      Y = XL(2,I)
      A  = DSQRT((X - XC(1,1))**2 + (Y - XC(1,2))**2)/RR
!     A  = SQRT((X - XC(1,1))**2 + (Y - XC(1,2))**2)/RR
      IF (A.LT.1.D0) THEN
       IF(A.GT.RI) QW(I) = (1.D0 - A)/(1.D0 - RI)
       LDOM = 1
      ELSE
       QW(I) = 0.D0
      END IF
 20   CONTINUE
      RETURN
      END

      SUBROUTINE DJFN(UL,SIG,QW,SWD,ACCEL,VEL,SHP,RHO,NEL,NDF,
     1 IDOM,DJ,DJM)
!---- 2 - D AREA INTEGRAND FOR J  STATIONARY CRACK
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /JAUX/ IAUX,STRAUX(4),EPSAUX(2),XCT(2),RAD,MODE
      DIMENSION UL(NDF,*),SIG(*),SHP(3,*),ACCEL(*),DJ(7,*),DJM(7,*),
     1 QW(*),VEL(*)
      DIMENSION EPSX(2),EPSY(2),ACC(2),DELQ(2)
      DIMENSION EPSDX(2),EPSDY(2),VELOC(2),DISP(2)
!---- COMPUTE GRADIENTS
          QQ = 0.D0
           DO 5 I=1,NEL
           QQ=QQ+SHP(3,I)*QW(I)
    5      CONTINUE
      DO 20 I = 1,NDF
        DISP(I) = 0.D0
      EPSX(I) = 0.D0
      EPSY(I) = 0.D0
      EPSDX(I) = 0.D0
      EPSDY(I) = 0.D0
      VELOC(I) = 0.D0
      ACC(I)  = 0.D0
      DELQ(I) = 0.D0
      DO 10 J = 1,NEL
      JI = (J-1)*NDF+I
           DISP(I)=DISP(I)+SHP(3,J)*UL(I,J)
      EPSX(I) = EPSX(I) + SHP(1,J)*UL(I,J)
      EPSY(I) = EPSY(I) + SHP(2,J)*UL(I,J)
      EPSDX(I) = EPSDX(I)+SHP(1,J)*VEL(JI)
      EPSDY(I) = EPSDY(I)+SHP(2,J)*VEL(JI)
      DELQ(I) = DELQ(I) + SHP(I,J)*QW(J)
   10 CONTINUE
   20 CONTINUE
!---- 'STATIC'
      IF(IAUX.EQ.0) THEN
      DJ(1,IDOM) = SIG(1)*EPSX(1)*DELQ(1) + SIG(2)*EPSX(2)*DELQ(2)
!    1           + SIG(4)*(EPSX(1)*DELQ(2) + EPSX(2)*DELQ(1))
     1           + SIG(4)*EPSX(1)*DELQ(2) + SIG(5)*EPSX(2)*DELQ(1)
!---- MODIFICATION-BEGIN ( HOOP STRESS CORRECTION REF: DELORENZI
!---- pp. 129-143, Vol.21, 1985 EFM)
           IF(MODE.EQ.2) THEN
           DJ(1,IDOM)=DJ(1,IDOM)+QQ*(DISP(1)/RAD*SIG(3)-SWD)/RAD
           ENDIF
!---- MODIFICATION ENDS
      ELSE
      DJ(1,IDOM) = SIG(1)*EPSAUX(1)*DELQ(1) + 
     1             SIG(2)*EPSAUX(2)*DELQ(2)
     1           + SIG(4)*(EPSAUX(1)*DELQ(2) + EPSAUX(2)*DELQ(1))
      DJ(1,IDOM) = DJ(1,IDOM) + STRAUX(1)*EPSX(1)*DELQ(1) + 
     1             STRAUX(2)*EPSX(2)*DELQ(2)
     1           + STRAUX(4)*(EPSX(1)*DELQ(2) + EPSX(2)*DELQ(1))
      END IF
      DJM(1,IDOM) = SIG(1)*EPSY(1)*DELQ(1) + SIG(2)*EPSY(2)*DELQ(2)
!    1           + SIG(4)*(EPSY(1)*DELQ(2) + EPSY(2)*DELQ(1))
     1           + SIG(4)*EPSY(1)*DELQ(2) + SIG(5)*EPSY(2)*DELQ(1)
      DJ(2,IDOM) = SWD*DELQ(1)
      DJM(2,IDOM) = SWD*DELQ(2)
!---- 'DYNAMIC'
      K = 1
      DO 30 I = 1,NEL
      ACC(1)  = ACC(1) + ACCEL(K)*SHP(3,I)
      ACC(2)  = ACC(2) + ACCEL(K+1)*SHP(3,I)
      VELOC(1) = VELOC(1)+VEL(K)*SHP(3,I)
      VELOC(2) = VELOC(2)+VEL(K+1)*SHP(3,I)
!     QQ = QQ + SHP(3,I)*QW(I)
   30 K = K+2
      DJ(3,IDOM) = RHO * QQ * (ACC(1)*EPSX(1) + ACC(2)*EPSX(2))
      DJM(3,IDOM) = RHO * QQ * (ACC(1)*EPSY(1) + ACC(2)*EPSY(2))
      EKE = 0.5*RHO*(VELOC(1)*VELOC(1) + VELOC(2)*VELOC(2))
      DJ(4,IDOM) = EKE*DELQ(1)
      DJM(4,IDOM) = EKE*DELQ(2)
      DJ(5,IDOM) = RHO*QQ*(VELOC(1)*EPSDX(1)+VELOC(2)*EPSDX(2))
      DJM(5,IDOM) =RHO*QQ*(VELOC(1)*EPSDY(1)+VELOC(2)*EPSDY(2))
!---- TOTAL CONTRIBUTION
      DJ(6,IDOM) = DJ(1,IDOM) - DJ(2,IDOM) + DJ(3,IDOM)
      DJM(6,IDOM) =DJM(1,IDOM) - DJM(2,IDOM) + DJM(3,IDOM)
      DJ(7,IDOM) = DJ(6,IDOM) - DJ(4,IDOM) - DJ(5,IDOM)
      DJM(7,IDOM) = DJM(6,IDOM)-DJM(4,IDOM)-DJM(5,IDOM)
      RETURN
      END



      SUBROUTINE QWT3(XL,QW,NDM,NEL,N,LDOM,NDI)

!---- 3 - D ENERGY RELEASE RATE FOR MODE I STATIONARY CRACK

!---- COMPUTE WEIGHTS QW

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /JINT/ NDOM(3),XC(8,4),DOM(4,10),DYNJ(8,7,10),
     1 DYNJM(8,7,10),TSWE,TKEE,TKEE1
      DIMENSION XL(NDM,*),QW(*)

      XIFN = DOM(4,N)
      IF (XIFN.GT.1.D0) GOTO 2
!---- BILINEAR FUNCTION (RECTANGULAR DOMAIN)
      XX = DOM(1,N)
      YY = DOM(2,N)
      IF (DOM(3,N).EQ.0.D0) THEN
       XI = 0.D0
       YI = 0.D0
      ELSE
       XI = 1.D0 - DOM(3,N)/XX
       YI = 1.D0 - DOM(3,N)/YY
      ENDIF
      DO 10 I = 1,NEL
      QW(I) = 1.D0
      IF (XC(NDI,3).GE.0.D0) THEN
       IF (XL(3,I).NE.XC(NDI,3)) QW(I) = 0.D0
       IF (XL(3,I).NE.XC(NDI,3)) GOTO 10
      ENDIF
      X = XL(1,I)
      Y = XL(2,I)
      A  = (X - XC(NDI,1))/XX
      B  = (Y - XC(NDI,2))/YY
      A = DABS(A)
!     A = ABS(A)
      B = DABS(B)
!     B = ABS(B)
      IF (A.LT.1.D0.AND.B.LT.1.D0) THEN
       IF (A.GT.XI) QW(I) = (1.D0 - A)/(1.D0 - XI)
       IF (B.GT.YI) QW(I) = QW(I)*(1.D0 - B)/(1.D0 - YI)
       LDOM = 1
      ELSE
       QW(I) = 0.D0
      END IF
 10   CONTINUE
      RETURN
!---- CONICAL FUNCTION (CIRCULAR DOMAIN)
    2 CONTINUE
      RR = DOM(1,N)
      IF(DOM(3,N).EQ.0.D0) THEN
       RI = 0.D0
      ELSE
       RI = 1.D0 - DOM(3,N)/RR
      ENDIF
      DO 20 I=1,NEL
      IF (XL(3,I).NE.XC(NDI,3)) QW(I) = 0.D0
      IF (XL(3,I).NE.XC(NDI,3)) GOTO 20
      QW(I) = 1.D0
      X = XL(1,I)
      Y = XL(2,I)
      A  = DSQRT((X - XC(NDI,1))**2 + (Y - XC(NDI,2))**2)/RR
!     A  = SQRT((X - XC(NDI,1))**2 + (Y - XC(NDI,2))**2)/RR
      IF(A.LT.1.D0) THEN
       IF(A.GT.RI) QW(I) = (1.D0 - A)/(1.D0 - RI)
       LDOM = 1
      ELSE
       QW(I) = 0.D0
      END IF
 20   CONTINUE
      RETURN
      END

      SUBROUTINE DJFN3(UL,SIG,QW,SWD,ACCEL,SHP,RHO,NEL,NDF,IDOM,DJ)
!---- AREA INTEGRAND FOR DYNAMIC J, STATIONARY CRACK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION UL(NDF,*),SIG(*),SHP(4,*),ACCEL(*),DJ(7,*),QW(*)
      DIMENSION EPSX(3),ACC(3),DELQ(3)
!---- COMPUTE GRADIENTS
      QQ = 0.D0
      DO 20 I = 1,NDF
      EPSX(I) = 0.D0
      ACC(I)  = 0.D0
      DELQ(I) = 0.D0
      DO 10 J = 1,NEL
      EPSX(I) = EPSX(I) + SHP(1,J)*UL(I,J)
      DELQ(I) = DELQ(I) + SHP(I,J)*QW(J)
   10 CONTINUE
   20 CONTINUE
!---- 'STATIC'
      DJ(1,IDOM) = SIG(1)*EPSX(1)*DELQ(1) + SIG(2)*EPSX(2)*DELQ(2)
     1           + SIG(3)*EPSX(3)*DELQ(3)
!    2           + SIG(4)*(EPSX(1)*DELQ(2) + EPSX(2)*DELQ(1))
!    3           + SIG(5)*(EPSX(2)*DELQ(3) + EPSX(3)*DELQ(2))
!    4           + SIG(6)*(EPSX(3)*DELQ(1) + EPSX(1)*DELQ(3))
     2           + SIG(4)*EPSX(1)*DELQ(2) + sig(7)*EPSX(2)*DELQ(1)
     4           + SIG(5)*EPSX(1)*DELQ(3) + sig(8)*EPSX(3)*DELQ(1)
     3           + SIG(6)*EPSX(2)*DELQ(3) + sig(9)*EPSX(3)*DELQ(2)
      DJ(2,IDOM) = SWD*DELQ(1)
      DJ(2,IDOM) = SWD*DELQ(1)
!---- 'DYNAMIC'
      K = 1
      DO 30 I = 1,NEL
      ACC(1) = ACC(1) + ACCEL(K  )*SHP(4,I)
      ACC(2) = ACC(2) + ACCEL(K+1)*SHP(4,I)
      ACC(3) = ACC(3) + ACCEL(K+2)*SHP(4,I)
      QQ = QQ + SHP(4,I)*QW(I)
   30 K = K + 3
      DJ(3,IDOM) = RHO*QQ*(ACC(1)*EPSX(1)+ACC(2)*EPSX(2)+ACC(3)*EPSX(3))
!---- TOTAL CONTRIBUTION
!     DJ(4,IDOM) = DJ(1,IDOM) - DJ(2,IDOM) + DJ(3,IDOM)
      DJ(6,IDOM) = DJ(1,IDOM) - DJ(2,IDOM) + DJ(3,IDOM)
      RETURN
      END


      SUBROUTINE JLINT(NP,TIME)
      implicit double precision (a-h,o-z)

!     J-INTEGRAL CALCULATION BY LINE INTEGRATION FOR CHARPY V NOTCH SPECIMEN.
!     READS THE VALUE OF EACH VARIABLE AT O/P QUAD PTS & USES THE AVERAGE VALUE
!     BETWEEN TWO QUAD PTS TO CALC THE VALUE OF THE INTEGRAL OVER THAT SEGMENT
!     AND SUM IT UP
!     INCORPORATED BY MANISH JHA 'OCTOBER 1993'
      COMMON /XJLINT/NNN,S11(50),S22(50),S12(50),X1(50),X2(50),
     1                DU12(50),DU22(50),WSE1(50)
      OPEN(UNIT=99,FILE='JINT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     # FORM='FORMATTED')
!     N0. OF SEGMENTS=NS
      NS = NP-1
      XJ = 0.D0
      DO 100 I=1,NS
!     CALCULATE SEGMENT LENGTH
      SD=((X1(I)-X1(I+1))**2) + ((X2(I)-X2(I+1))**2)
      D=DSQRT(SD)
!     CALCULATE SEGMNT NORMAL DIRN COSINES
      XN1=(X2(I+1)-X2(I))/D
      XN2=-(X1(I+1)-X1(I))/D
      DEN=DSQRT((XN1**2)+(XN2**2))
      XN1=XN1/DEN
      XN2=XN2/DEN
!     CALCULATE AVERAGES OVER THE SEGMENT
      AW=(WSE1(I)+WSE1(I+1))/2.D0
      ADU12=(DU12(I)+DU12(I+1))/2.D0
      ADU22=(DU22(I)+DU22(I+1))/2.D0
      AS11=(S11(I)+S11(I+1))/2.D0
      AS12=(S12(I)+S12(I+1))/2.D0
      AS22=(S22(I)+S22(I+1))/2.D0
      T1=AW * XN2
      T2=AS11 * XN1 * ADU12
      T3=AS12 * XN2 * ADU12
      T4=AS22 * XN2 * ADU22
      T5=AS12 * XN1 * ADU22
      XJ=XJ+(T1-(T2+T3+T4+T5))*D
!     FIRST TERM BETWEEN BOUNDARY AND 1ST QUAD PT.,ASSUMING QUAD PT.
!     VALUES OVER THIS SEGMENT
      IF(I.EQ.1) THEN
!     D=DIST BETWWEN THE FIRST QUAD PT AND THE INCLINED NOTCH FACE
      D=(1/(2.613))*((2.414*X1(1))+X2(1)-5.0790)
      D=DABS(D)
      FT=(WSE1(1)*XN2-(S11(1)*XN1*DU12(1)+S12(1)*XN2*DU12(1)+
     #     S12(1)*XN1*DU22(1)+S22(1)*XN2*DU22(1)))*D
      XJ=XJ+FT
!     LAST TERM BETWWEN BOUNDARY AND LAST QUAD PT.
      ELSE IF(I.EQ.NS) THEN
      TL=(WSE1(NP)-(S12(NP)*DU12(NP)+S22(NP)*DU22(NP)))*X1(NP)
      XJ=XJ+TL
      ENDIF
100   CONTINUE
!     BECAUSE OF SYMMETRICITY OF THE PROBLEM
      XJ=XJ*2
      WRITE(99,12) TIME, XJ
12    FORMAT('TIME',1X,E20.8,2X,'J-INT ',1X,E20.8)
      RETURN
      END

      subroutine shape2tc(s,xel,shp,xsj,flag)
      implicit double precision (a-h,o-z)

!  6-node composite triangular element (piecewise linear)

      logical flag
      dimension s(*),xel(2,*),shp(3,*),x(2,6)
      dimension dnds(3,6),dxds(2,3),xjac(3,2)
      dimension nsub(3,4),xc(4),yc(4),a(4),b(4),c(4)
      data nsub/1,4,6, 4,2,5, 6,5,3, 4,5,6/
      data tvol/0.d0/

!     write(6,*)'*** composite shape fn derivatives'

      shp(3,1) = s(1)*(s(1) - s(2) - s(3))
      shp(3,2) = s(2)*(s(2) - s(1) - s(3))
      shp(3,3) = s(3)*(s(3) - s(1) - s(2))
      shp(3,4) = 4.d0*s(1)*s(2)
      shp(3,5) = 4.d0*s(2)*s(3)
      shp(3,6) = 4.d0*s(3)*s(1)

      if (.not. flag) return



      xx = (xel(1,1)+xel(1,2)+xel(1,3))/3.d0
      yy = (xel(2,1)+xel(2,2)+xel(2,3))/3.d0
      do i = 1,6
         x(1,i) = xel(1,i) - xx
         x(2,i) = xel(2,i) - yy
      end do

      c11 = 0.d0
      c12 = 0.d0
      c13 = 0.d0
      c22 = 0.d0
      c23 = 0.d0
      c33 = 0.d0
      do isub = 1,4
         i1 = nsub(1,isub)
         i2 = nsub(2,isub)
         i3 = nsub(3,isub)
         asub =( x(1,i1)*x(2,i2) - x(1,i2)*x(2,i1)
     1         + x(1,i2)*x(2,i3) - x(1,i3)*x(2,i2)
     1         + x(1,i3)*x(2,i1) - x(1,i1)*x(2,i3) ) / 2.d0
         xc(isub) = (x(1,i1) + x(1,i2) + x(1,i3))/3.d0
         yc(isub) = (x(2,i1) + x(2,i2) + x(2,i3))/3.d0
         c11 = c11 + asub
         c12 = c12 + asub*xc(isub)
         c13 = c13 + asub*yc(isub)
         c22 = c22 + asub*( x(1,i1)*x(1,i1) + x(1,i2)*x(1,i2) + 
     1                      x(1,i3)*x(1,i3) + x(1,i1)*x(1,i2) +
     1                      x(1,i2)*x(1,i3) + x(1,i3)*x(1,i1) )/12.d0
         c23 = c23 + asub*( (x(1,i1)*x(2,i1) + x(1,i2)*x(2,i2) +
     1                       x(1,i3)*x(2,i3)) / 12.d0 +
     1                      (x(1,i1)*x(2,i2) + x(1,i2)*x(2,i1) +
     1                       x(1,i2)*x(2,i3) + x(1,i3)*x(2,i2) +
     1                       x(1,i3)*x(2,i1) + x(1,i1)*x(2,i3))/24.d0 )
         c33 = c33 + asub*( x(2,i1)*x(2,i1) + x(2,i2)*x(2,i2) + 
     1                      x(2,i3)*x(2,i3) + x(2,i1)*x(2,i2) +
     1                      x(2,i2)*x(2,i3) + x(2,i3)*x(2,i1) )/12.d0
      end do
      det = c11*c22*c33 + 2.d0*c12*c13*c23 - c11*c23**2 - c22*c13**2
     1    - c33*c12**2

!     write(6,1001)c11,c22,c33,c12,c13,c23,det
!1001 format('cij=',1p6g12.4,/,'det=',1pg12.4)      

      d1 = c22*c33 - c23*c23
      d2 = c13*c23 - c12*c33
      d3 = c12*c23 - c13*c22
      d4 = c13*c23 - c12*c33
      d5 = c11*c33 - c13*c13
      d6 = c12*c13 - c11*c23
      d7 = c12*c23 - c13*c22
      d8 = c12*c13 - c11*c23
      d9 = c11*c22 - c12*c12
      do isub = 1,4
         a(isub) = ( d1 + d2*xc(isub) + d3*yc(isub) ) /det
         b(isub) = ( d4 + d5*xc(isub) + d6*yc(isub) ) /det
         c(isub) = ( d7 + d8*xc(isub) + d9*yc(isub) ) /det
!        write(6,1002)isub,a(isub),b(isub),c(isub)
!1002    format('isub,a,b,c=',i4,1p3g12.4)
      end do

      xx = x(1,1)*s(1) + x(1,2)*s(2) + x(1,3)*s(3) 
      yy = x(2,1)*s(1) + x(2,2)*s(2) + x(2,3)*s(3) 
!     write(6,*)'xx,yy=',xx,yy

      shp(1,1) = ( x(2,4)-x(2,6) ) * ( a(1)+b(1)*xx+c(1)*yy )
      shp(1,2) = ( x(2,5)-x(2,4) ) * ( a(2)+b(2)*xx+c(2)*yy )
      shp(1,3) = ( x(2,6)-x(2,5) ) * ( a(3)+b(3)*xx+c(3)*yy )
      shp(1,4) = ( x(2,6)-x(2,1) ) * ( a(1)+b(1)*xx+c(1)*yy ) +
     1           ( x(2,2)-x(2,5) ) * ( a(2)+b(2)*xx+c(2)*yy ) +
     1           ( x(2,5)-x(2,6) ) * ( a(4)+b(4)*xx+c(4)*yy ) 
      shp(1,5) = ( x(2,4)-x(2,2) ) * ( a(2)+b(2)*xx+c(2)*yy ) +
     1           ( x(2,3)-x(2,6) ) * ( a(3)+b(3)*xx+c(3)*yy ) +
     1           ( x(2,6)-x(2,4) ) * ( a(4)+b(4)*xx+c(4)*yy ) 
      shp(1,6) = ( x(2,1)-x(2,4) ) * ( a(1)+b(1)*xx+c(1)*yy ) +
     1           ( x(2,5)-x(2,3) ) * ( a(3)+b(3)*xx+c(3)*yy ) +
     1           ( x(2,4)-x(2,5) ) * ( a(4)+b(4)*xx+c(4)*yy ) 

      shp(2,1) = ( x(1,6)-x(1,4) ) * ( a(1)+b(1)*xx+c(1)*yy )
      shp(2,2) = ( x(1,4)-x(1,5) ) * ( a(2)+b(2)*xx+c(2)*yy )
      shp(2,3) = ( x(1,5)-x(1,6) ) * ( a(3)+b(3)*xx+c(3)*yy )
      shp(2,4) = ( x(1,1)-x(1,6) ) * ( a(1)+b(1)*xx+c(1)*yy ) +
     1           ( x(1,5)-x(1,2) ) * ( a(2)+b(2)*xx+c(2)*yy ) +
     1           ( x(1,6)-x(1,5) ) * ( a(4)+b(4)*xx+c(4)*yy ) 
      shp(2,5) = ( x(1,2)-x(1,4) ) * ( a(2)+b(2)*xx+c(2)*yy ) +
     1           ( x(1,6)-x(1,3) ) * ( a(3)+b(3)*xx+c(3)*yy ) +
     1           ( x(1,4)-x(1,6) ) * ( a(4)+b(4)*xx+c(4)*yy ) 
      shp(2,6) = ( x(1,4)-x(1,1) ) * ( a(1)+b(1)*xx+c(1)*yy ) +
     1           ( x(1,3)-x(1,5) ) * ( a(3)+b(3)*xx+c(3)*yy ) +
     1           ( x(1,5)-x(1,4) ) * ( a(4)+b(4)*xx+c(4)*yy ) 

      do j = 1,6
         do i = 1,2
            shp(i,j) = shp(i,j)/2.d0
         end do
      end do

!     write(6,1003)((shp(i,j),j=1,6),i=1,2)
!1003 format('shp()=',1p6g12.4)

      xsj = c11/3.d0 

      tvol = tvol + xsj
!     write(6,*)'xsj,tvol=',xsj,tvol

      return
      end


      subroutine shape2t(s,x,shp,xsj,flag)
      implicit double precision (a-h,o-z)

!---- shape fns and derivatives for 6-node triangle (quadratic)

      logical flag
      dimension s(*),x(2,*),shp(3,*)
      dimension dnds(3,6),dxds(2,3),xjac(3,2)

      shp(3,1) = s(1)*(s(1) - s(2) - s(3))
      shp(3,2) = s(2)*(s(2) - s(1) - s(3))
      shp(3,3) = s(3)*(s(3) - s(1) - s(2))
      shp(3,4) = 4.d0*s(1)*s(2)
      shp(3,5) = 4.d0*s(2)*s(3)
      shp(3,6) = 4.d0*s(3)*s(1)

      if (.not. flag) return

      dnds(1,1) = 2.d0*s(1) - s(2) - s(3)
      dnds(2,1) = - s(1)
      dnds(3,1) = - s(1)

      dnds(1,2) = - s(2)
      dnds(2,2) = 2.d0*s(2) - s(1) - s(3)
      dnds(3,2) = - s(2)

      dnds(1,3) = - s(3)
      dnds(2,3) = - s(3)
      dnds(3,3) = 2.d0*s(3) - s(1) - s(2)

      dnds(1,4) = 4.d0*s(2)
      dnds(2,4) = 4.d0*s(1)
      dnds(3,4) = 0.

      dnds(1,5) = 0.
      dnds(2,5) = 4.d0*s(3)
      dnds(3,5) = 4.d0*s(2)

      dnds(1,6) = 4.d0*s(3)
      dnds(2,6) = 0.
      dnds(3,6) = 4.d0*s(1)

      do 10 j = 1,3
      do 10 i = 1,2
      dxds(i,j) = 0.
      do 10 k = 1,6
   10 dxds(i,j) = dxds(i,j) + x(i,k)*dnds(j,k)

      xsj = dxds(1,2)*dxds(2,3) - dxds(1,3)*dxds(2,2)
     1      + dxds(1,3)*dxds(2,1) - dxds(1,1)*dxds(2,3)
     2      + dxds(1,1)*dxds(2,2) - dxds(1,2)*dxds(2,1)

      xjac(1,1) = (dxds(2,2) - dxds(2,3))/xsj
      xjac(2,1) = (dxds(2,3) - dxds(2,1))/xsj
      xjac(3,1) = (dxds(2,1) - dxds(2,2))/xsj
      xjac(1,2) = (dxds(1,3) - dxds(1,2))/xsj
      xjac(2,2) = (dxds(1,1) - dxds(1,3))/xsj
      xjac(3,2) = (dxds(1,2) - dxds(1,1))/xsj

      xsj = xsj/6.d0

      do 20 j = 1,6
      do 20 i = 1,2
      shp(i,j) = 0.
      do 20 k = 1,3
   20 shp(i,j) = shp(i,j) + dnds(k,j)*xjac(k,i)
      return
      end

      subroutine shape3t(s,x,shp,xsj,flag)
      implicit double precision (a-h,o-z)

!---- shape fns and derivatives for 10-node tetrahedran (quadratic)

      logical flag
      dimension s(*),x(3,*),shp(4,*)
      dimension dnds(3,10),dxds(3,3),xjac(3,3)

!  Natural coordinates :
!     s(1) = u = 1-r-s-t
!     s(2) = r
!     s(3) = s
!     s(4) = t
!     Ni
      shp(4,1) = s(1)*(2.d0*s(1)-1.d0)
      shp(4,2) = s(2)*(2.d0*s(2)-1.d0)
      shp(4,3) = s(3)*(2.d0*s(3)-1.d0)
      shp(4,4) = s(4)*(2.d0*s(4)-1.d0)
      shp(4,5) = 4.d0*s(1)*s(2)
      shp(4,6) = 4.d0*s(2)*s(3)
      shp(4,7) = 4.d0*s(3)*s(1)
      shp(4,8) = 4.d0*s(1)*s(4)
      shp(4,9) = 4.d0*s(2)*s(4)
      shp(4,10)= 4.d0*s(3)*s(4)

      if (.not. flag) return
! dNi/dr
      dnds(1,1) = 1.d0-4.d0*s(1)
      dnds(1,2) = -1.d0+4.d0*s(2)
      dnds(1,3) = 0.d0
      dnds(1,4) = 0.d0
      dnds(1,5) = 4.d0*(s(1)-s(2))
      dnds(1,6) = 4.d0*s(3)
      dnds(1,7) = -4.d0*s(3)
      dnds(1,8) = -4.d0*s(4)
      dnds(1,9) = 4.d0*s(4)
      dnds(1,10)= 0.d0
! dNi/ds
      dnds(2,1) = 1.d0-4.d0*s(1)
      dnds(2,2) = 0.d0
      dnds(2,3) = -1.d0+4.d0*s(3)
      dnds(2,4) = 0.d0
      dnds(2,5) = -4.d0*s(2)
      dnds(2,6) = 4.d0*s(2)
      dnds(2,7) = 4.d0*(s(1)-s(3))
      dnds(2,8) = -4.d0*s(4)
      dnds(2,9) = 0.d0
      dnds(2,10)= 4.d0*s(4)
! dNi/dt
      dnds(3,1) = 1.d0-4.d0*s(1)
      dnds(3,2) = 0.d0
      dnds(3,3) = 0.d0
      dnds(3,4) = -1.d0+4.d0*s(4)
      dnds(3,5) = -4.d0*s(2)
      dnds(3,6) = 0.d0
      dnds(3,7) = -4.d0*s(3)
      dnds(3,8) = 4.d0*(s(1)-s(4))
      dnds(3,9) = 4.d0*s(2)
      dnds(3,10)= 4.d0*s(3)

      do 10 j = 1,3
         do 10 i = 1,3
            dxds(i,j) = 0.d0
            do 10 k = 1,10
               dxds(i,j) = dxds(i,j) + x(i,k)*dnds(j,k)
  10  continue

      xsj =   dxds(1,1)*(dxds(2,2)*dxds(3,3)-dxds(3,2)*dxds(2,3))
     1      - dxds(1,2)*(dxds(2,1)*dxds(3,3)-dxds(3,1)*dxds(2,3))
     2      + dxds(1,3)*(dxds(2,1)*dxds(3,2)-dxds(3,1)*dxds(2,2))

      xjac(1,1) =  (dxds(2,2)*dxds(3,3)-dxds(3,2)*dxds(2,3))/xsj
      xjac(2,1) = -(dxds(2,1)*dxds(3,3)-dxds(3,1)*dxds(2,3))/xsj
      xjac(3,1) =  (dxds(2,1)*dxds(3,2)-dxds(3,1)*dxds(2,2))/xsj
      xjac(1,2) = -(dxds(1,2)*dxds(3,3)-dxds(3,2)*dxds(1,3))/xsj
      xjac(2,2) =  (dxds(1,1)*dxds(3,3)-dxds(3,1)*dxds(1,3))/xsj
      xjac(3,2) = -(dxds(1,1)*dxds(3,2)-dxds(3,1)*dxds(1,2))/xsj
      xjac(1,3) =  (dxds(1,2)*dxds(2,3)-dxds(2,2)*dxds(1,3))/xsj
      xjac(2,3) = -(dxds(1,1)*dxds(2,3)-dxds(2,1)*dxds(1,3))/xsj
      xjac(3,3) =  (dxds(1,1)*dxds(2,2)-dxds(2,1)*dxds(1,2))/xsj

      xsj = xsj/6.d0

      do 20 j = 1,10
         do 20 i = 1,3
            shp(i,j) = 0.d0
            do 20 k = 1,3
               shp(i,j) = shp(i,j) + dnds(k,j)*xjac(k,i)
  20  continue
      return
      end

