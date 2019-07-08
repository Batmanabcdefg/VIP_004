      SUBROUTINE RNSTCN(UL,XL,TL,P,S,D,X,F,STR,EPS,Q,VEL,ACCEL,
     1 HB1,FB1,HB2,FB2,NDF,NDM,NST,NSTR,NQ,NELGAM,NGAMA,N2GAMA,N4GAMA,          
     2 NRING,NBAND1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---- PERFORMS RING-BY-RING STATIC CONDENSATION OF OUTER MESH USED
C---- IN SIMULATION OF MIXED-MODE SSY PROBLEM WITH DISPLACEMENTS (K1-K2 FIELD)
C---- PRESCRIBED ON OUTERMOST RING OF NODES
c
c---- modified on 26/2/98 to treat arbitrary numbering of outer
c---- ring of nodes. Note NSET(2,I+1), I=1,ngama contains the
c---- the global node numbers of the outer ring of nodes in clockwise
c---- sense starting from upper crack flank and ending in lower crack
c---- flank.
C
      DIMENSION UL(NDF,*),XL(NDM,*),TL(*),P(*),S(NST,*),D(40,*),
     1 X(NDM,*),F(*),STR(NSTR,*),EPS(NSTR,*),Q(NQ,*),VEL(*),ACCEL(*)
      DIMENSION HB1(N4GAMA,*),HB2(N4GAMA,*),FB1(*),FB2(*)
      DIMENSION LM(4),IDL(2,4)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL
      COMMON /ISTCND/ IFLSTC
      COMMON /PRTSET/ NSET(9,201)
      COMMON /MAINN/ M(25000000) 
      COMMON /MDATA/ NN,N0,N1,N2,N3,N301,N302,N4,N5,N6,N7,N8,N9,N91,
     1 N10,N11,N12,N13,N14,N15,N16,N161,N162,N163,N164,N165,N166,
     1 N17,N18,N181,N182,N183,N184,N185,N186,N19,N20
      DATA PI/3.141592653589793/
      MA = 1
      MODE = D(3,MA)
      NEN1 = NEN +1
c--------- modified  on 26/2/98 ------------------------
c---- set the profile of the stiffness matrix   
      CALL PROFL1(M(N12),M(N7),M(N9),NDF,NEN1,NAD,2)
c--------- 26/2/98  ------------------------------------
C---- FOLL ARE FOR RATE INDEP PROG: PL2D.F
c     emod = d(6,ma)
c     POIS = D(7,MA)
      emod = 2500.0
      pois = 0.25
c     EG = D(8,MA)
      eg = emod/(2.d0*(1.d0+pois))
c     XSI = D(13,MA)
      xsi = 0.d0
C---- FOR FPL2DNEW.F   XSI = D(21,MA)
C---- FOLL ARE FOR RATE DEP PROG: FPL2DSUM.F
C     POIS = D(21,MA)
C     EG = D(8,MA)
C     XSI = D(22,MA)
C--------------------------------------------------------
C     XK2 = SIN(XSI*PI/180.E0)
C     XK1 = COS(XSI*PI/180.E0)
      XK2 = DSIN(XSI*PI/180.D0)
      XK1 = DCOS(XSI*PI/180.D0)
C---- modified on 26/2/98 ---------------------------
C     R0 = ABS(X(1,NUMNP))
C     RNTCH = ABS(X(2,NUMNP))
      NLAST = NSET(2, NGAMA+1)
      R0 = DABS(X(1,NLAST))
c     RNTCH = DABS(X(2,NLAST))
      RNTCH = DABS(X(2,NSET(2,2)))
c----------  26/2/98  ---------------------------------------
C     R0 = DABS(X(1,NUMNP))
C     RNTCH = DABS(X(2,NUMNP))
C     DTHT = 2.E0*PI/FLOAT(NELGAM)
c---- following is for full ring (mixed-mode case).
      DTHT = 2.D0*PI/DFLOAT(NELGAM)
c---- following is for mode I symm case (1/2 ring)
c     DTHT = PI/DFLOAT(NELGAM)
c
C     VAL1 = FLOAT(NRING)*LOG(1.E0+DTHT)
      VAL1 = DFLOAT(NRING)*DLOG(1.D0+DTHT)
C     TOTRAD = EXP(VAL1)      
      TOTRAD = DEXP(VAL1)      
      RAD2=TOTRAD*R0
C     CONST=SQRT(RAD2/(2.0*PI))/(2.0*EG)
      CONST=DSQRT(RAD2/(2.D0*PI))/(2.D0*EG)
      const1 = (1.d0 - pois*pois)*rad2/emod
      const2 = -pois*(1.d0+pois)*rad2/emod
      NRING1=NRING-1
      EKAP = 3.D0 -4.D0*POIS
      IF(MODE.EQ.3)
     1 EKAP=(3.D0-POIS)/(1.D0+POIS)
      NEL = 4
C---- FOLL STATEMENT IS FOR RATE INDEP CASE: FPL2DNEW
      IEL = 1
C---- FOLL IS FOR RATE DEPEND CASE: FPL2DSUM
C     IEL = 2
      RAD1=RAD2/(1.D0+DTHT)
C***************************
      WRITE(6,9001) NELGAM, NRING,NGAMA,N2GAMA,N4GAMA,NBAND1
9001  FORMAT(2X,'NELGAM NRING NGAMA N2GAMA N4GAMA NBAND1',1X,6I5)
c------------  modified on 26/2/98 --------------
      WRITE(6,9002) R0, RNTCH, RAD1,RAD2,DTHT
9002  FORMAT(2X,'R0 RNTCH RAD1 RAD2 DTHT',1X,5E13.5)
c--------------  26/2/98  ---------------------------
      WRITE(6,9003) POIS,EG,CONST
9003  FORMAT(2X,'POIS EG CONST',1X,3E13.5)
C***************************************
      DO 101 KKK=1,2
      RAD2=TOTRAD*R0
      RAD1=RAD2/(1.D0+DTHT)
c      KKK = 1
C----  KKK=1 : ELASTIC MIXED-MODE K-FIELD IS IMPOSED (RESULTS ADDED TO F)
C----  KKK=2:  T-STRESS TERM IMPOSED  (RESULTS RETAINED IN FB1)
      IRING=1
      IFLSTC = 1
      CALL PZERO(HB1,N4GAMA*NBAND1)
      CALL PZERO(FB1,N4GAMA)
C *** LOOP OVER ELEMENTS (4-NODED QUADS).
      DO 100 N=1,NELGAM
      NELMT = N 
C *** FORM LM ARRAY
      LM(1)=N+1
      LM(2)=N
      LM(3)=NGAMA+N
      LM(4)=NGAMA+N+1
C *** FORM IDL ARRAY ***
      DO 110 I=1,2
      IDL(1,I)=0
      IDL(2,I)=0
110   CONTINUE
      II=3
      DO 120 I=3,4
      II=II-1
      IDL(1,I)=2*LM(II)-1
      IDL(2,I)=2*LM(II)
120   CONTINUE
c---- following is for mode I symm case (1/2 ring); comment out for full ring
c---- for last element of the ring (i.e., adjacant to theta = 0 or symm line
c---- restrain y dof
c     if(n.eq.nelgam) then
c     idl(2,1) = -idl(2,1)
c     idl(2,4) = -idl(2,4)
c---- following is for anti-symm mode II case (1/2 ring)
c     idl(1,1) = -idl(1,1)
c     idl(1,4) = -idl(1,4)
c     end if
c
C *** FORM XL ***
      CALL CORELN(XL,RAD1,RAD2,THT1,THT2,RNTCH,DTHT,N,NDM,NELGAM)
      DO 121 I=1,4
      DO 121 J=1,2
121   UL(J,I) = 0.D0
C---- FORM THE ELEMENT ELASTIC STIFFNESS MATRIX S(NST,NST)
      CALL PZERO(VEL,NST)
      CALL PZERO(ACCEL,NST)
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,NUMEL),EPS(1,NUMEL),Q(1,NUMEL),
     1 LM,TL,S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,3)
C
      IF(KKK.EQ.1) THEN
C *** FORM UL FOR OUTERMOST NODE RING USING (K1-K2) FIELD
C---- FIRST THE MODE-I PART
      UL(1,1)=((EKAP-0.5)*DCOS(0.5*THT1)-0.5*DCOS(1.5*THT1))*CONST*XK1
      UL(2,1)=((EKAP+0.5)*DSIN(0.5*THT1)-0.5*DSIN(1.5*THT1))*CONST*XK1
      UL(1,2)=((EKAP-0.5)*DCOS(0.5*THT2)-0.5*DCOS(1.5*THT2))*CONST*XK1
      UL(2,2)=((EKAP+0.5)*DSIN(0.5*THT2)-0.5*DSIN(1.5*THT2))*CONST*XK1
C---- NOW FOR MODE-II PART
      UL(1,1)= UL(1,1) +
     1 ((EKAP+1.5)*DSIN(0.5*THT1)+0.5*DSIN(1.5*THT1))*CONST*XK2
      UL(2,1)= UL(2,1) -
     1 ((EKAP-1.5)*DCOS(0.5*THT1)+0.5*DCOS(1.5*THT1))*CONST*XK2
      UL(1,2)= UL(1,2) +
     1 ((EKAP+1.5)*DSIN(0.5*THT2)+0.5*DSIN(1.5*THT2))*CONST*XK2
      UL(2,2)= UL(2,2) -
     1 ((EKAP-1.5)*DCOS(0.5*THT2)+0.5*DCOS(1.5*THT2))*CONST*XK2
      ELSEIF(KKK.EQ.2) THEN
C *** FORM UL FOR OUTERMOST RING USING T-STRESS TERM (WITH T=1)
      ul(1,1) = const1*dcos(tht1)
      ul(2,1) = const2*dsin(tht1)
      ul(1,2) = const1*dcos(tht2)
      ul(2,2) = const2*dsin(tht2)
      ENDIF
C---- FORM P DUE TO SPECIFIED DISPLACEMENT
      CALL PZERO(P,NST)
      DO 122 I=5,8
      P(I) = -S(I,1)*UL(1,1) - S(I,2)*UL(2,1) - S(I,3)*UL(1,2)
     1 - S(I,4)*UL(2,2)      
122   CONTINUE
      CALL ASSMBL(HB1,FB1,S,P,IDL,N4GAMA,NBAND1,8,1,1)
100   CONTINUE
C *** LOOP OVER RINGS.IRING=2,NRING-1
      DO 200 IRING=2,NRING1
      CALL PZERO(HB2,N4GAMA*NBAND1)
      CALL PZERO(FB2,N4GAMA)
      RAD2=RAD1
      RAD1=RAD2/(1.D0+DTHT)
C *** LOOP OVER ELEMENTS (4-NODED QUADS) IN RING. 
      DO 210 N=1,NELGAM
      NELMT = N
C *** FORM LM ARRAY ***
      LM(1)=N+1
      LM(2)=N
      LM(3)=NGAMA+N
      LM(4)=NGAMA+N+1
C *** FORM IDL ARRAY ***
      DO 220 I=1,4
      IDL(1,I)=2*LM(I)-1
      IDL(2,I)=2*LM(I)
220   CONTINUE
c---- following is for mode I symm case (1/2 ring); comment out for full ring
c---- for last element of the ring (i.e., adjacant to theta = 0 or symm line
c---- restrain y dof
c     if(n.eq.nelgam) then
c     idl(2,1) = -idl(2,1)
c     idl(2,4) = -idl(2,4)
c---- following is for anti-symm mode II case (1/2 ring)
c     idl(1,1) = -idl(1,1)
c     idl(1,4) = -idl(1,4)
c     end if
c
C *** FORM XL ***
      CALL CORELN(XL,RAD1,RAD2,THT1,THT2,RNTCH,DTHT,N,NDM,NELGAM)
      DO 221 I=1,4
      DO 221 J=1,2
221   UL(J,I) = 0.D0
C---- FORM THE ELEMENT ELASTIC STIFFNESS MATRIX S(NST,NST)
      CALL PZERO(VEL,NST)
      CALL PZERO(ACCEL,NST)
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,NUMEL),EPS(1,NUMEL),Q(1,NUMEL),
     1 LM,TL,S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,3)
      CALL PZERO(P,NST)
      CALL ASSMBL(HB2,FB2,S,P,IDL,N4GAMA,NBAND1,8,1,0)
210   CONTINUE
C *** TRANSFER UPPER PART OF HB1 TO UPPER PART OF HB2 ***
C *** TRANSFER U.P. OF FB1 TO U.P. OF FB2 ***
      DO 240 I=1,N2GAMA
      FB2(I)=FB2(I)+FB1(I)
      DO 250 J=1,NBAND1
      HB2(I,J)=HB2(I,J)+HB1(I,J)
250   CONTINUE
240   CONTINUE
C---- PERFORM PARTIAL GAUSS REDUCTION OF UPPER PART (UP TO ROW OF N2GAMA)
C---- OF HB2 AND FB2
      CALL BSOLV2(HB2,FB2,N4GAMA,NBAND1,N2GAMA)
      CALL PZERO(HB1,N4GAMA*NBAND1)
      CALL PZERO(FB1,N4GAMA)
C *** TRANSFER LOWER PART OF HB2 TO UPPER PART OF HB1 ***
C *** TRANSFER L.P. OF FB2 TO U.P. OF FB1 ***
      DO 260 I=1,N2GAMA
      FB1(I)=FB2(N2GAMA+I)
      DO 260 J=1,NBAND1
      HB1(I,J)=HB2(N2GAMA+I,J)
260   CONTINUE
200   CONTINUE
      IRING =NRING
      CALL PZERO(HB2,N4GAMA*NBAND1)
      CALL PZERO(FB2,N4GAMA)
      RAD2=RAD1
      RAD1=R0
      write(6,9101) rad2, rad1
9101  format(10x,'rad2, rad1 for last ring',2e13.5)
c----- modified on 26/2/98 ----------------
C     IND1=NUMNP-NGAMA
c------------------- 26/2/98 -----------------
C *** LOOP OVER ELEMENTS (4-NODED QUADS) IN RING. 
      DO 310 N=1,NELGAM
      NELMT = N
c--------  modified on 26/2/98-----------------
C     IND1=IND1+1
C     IND2 = IND1 + 1
      IND1 = NSET(2, N+1)
      IND2 = NSET(2, N+2)
c---------------  26/2/98  -------------------
C *** FORM LM ARRAY
      LM(1)=N+1
      LM(2)=N
      LM(3)=NGAMA+N
      LM(4)=NGAMA+N+1
C *** FORM IDL ARRAY
      DO 320 I=1,4
      IDL(1,I)=2*LM(I)-1
      IDL(2,I)=2*LM(I)
320   CONTINUE
C *** FORM XL
c---- following is for mode I symm case (1/2 ring); comment out for full ring.
c---- for last element of the ring (i.e., adjacant to theta = 0 or symm line
c---- restrain y dof
c     if(n.eq.nelgam) then
c     idl(2,1) = -idl(2,1)
c     idl(2,4) = -idl(2,4)
c---- following is for anti-symm mode II case (1/2 ring)
c     idl(1,1) = -idl(1,1)
c     idl(1,4) = -idl(1,4)
c     end if
c---- following two statements are for mixed-mode case (full ring)
      THT1=DFLOAT((NELGAM-N))*DTHT - PI
      THT2=DFLOAT((NELGAM-N+1))*DTHT - PI
c---- following two statements are for mode I symm case (1/2 ring)
c     tht1 = dfloat((nelgam-n))*dtht
c     tht2 = dfloat((nelgam-n+1))*dtht
c
      XL(1,1)=RAD2*DCOS(THT1)
      XL(2,1)=RAD2*DSIN(THT1)
      XL(1,2)=RAD2*DCOS(THT2)
      XL(2,2)=RAD2*DSIN(THT2)
      XL(1,3)=X(1,IND1)
      XL(2,3)=X(2,IND1)
c----- modified on 26/2/98 ----------
      XL(1,4)=X(1,IND2)
      XL(2,4)=X(2,IND2)
c-----------   26/2/98  --------------
c---- following two statements apply to notch - full ring; comment out for
c     sharp crack.
      IF(N.EQ.1) XL(2,2) = RNTCH
      IF(N.EQ.NELGAM) XL(2,1) = -RNTCH      
c
      DO 321 I=1,4
      DO 321 J=1,2
321   UL(J,I) = 0.D0
C---- FORM THE ELEMENT ELASTIC STIFFNESS MATRIX S(NST,NST)
      CALL PZERO(VEL,NST)
      CALL PZERO(ACCEL,NST)
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,NUMEL),EPS(1,NUMEL),Q(1,NUMEL),
     1 LM,TL,S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,3)
      CALL PZERO(P,NST)
      CALL ASSMBL(HB2,FB2,S,P,IDL,N4GAMA,NBAND1,8,1,0)
310   CONTINUE
C *** TRANSFER U.P. OF HB1 TO U.P.OF HB2
C *** TRANSFER U.P. OF FB1 TO U.P. OF FB2
      DO 340 I=1,N2GAMA
      FB2(I)=FB2(I)+FB1(I)
      DO 350 J=1,NBAND1
      HB2(I,J)=HB2(I,J)+HB1(I,J)
350   CONTINUE
340   CONTINUE
C---- PERFORM PARTIAL FORWARD GAUSS REDUCTION OF UPPER PART OF HB2 AND FB2
      CALL BSOLV2(HB2,FB2,N4GAMA,NBAND1,N2GAMA)
      CALL PZERO(HB1,N4GAMA*NBAND1)
      CALL PZERO(FB1,N4GAMA)
C *** TRANSFER L.P. OF HB2 TO U.P. OF HB1
C *** TRANSFER L.P. OF FB2 TO U.P. OF FB1
      DO 360 I=1,N2GAMA
      FB1(I)=FB2(N2GAMA+I)
      DO 360 J=1,NBAND1
      HB1(I,J)=HB2(N2GAMA+I,J)
360   CONTINUE
C---- TRANSFER CONDENSED FORCE VECTOR CORR TO NODES ON OUTERMOST RADIUS OF
C---- ACTIVE MESH TO F ONLY FOR KKK=1 (K-FIELD TERMS). NOTE THAT THE
C---- T-STRESS TERMS WILL BE AVAILABLE IN CONDENSED FORCE VECTOR FB1 WHICH
C---- WILL BE RETURNED TO MAIN PROGRAM AFTER EXECUTING KKK=2 LOOP. THIS
C---- HAS TO BE ACCOUNTED FOR IN FADD WHILE ASSEMBLING FORCE VECTORS.
      IF(KKK.EQ.1) THEN
c---- modified on 26/2/98 -------------
c---- comment out the following if you are inputing concentrated nodal
c---- loads through  FORC  macro
      CALL PZERO(F,NEQ)
c     IEQ=NEQ-N2GAMA
      DO 370 I=1,N2GAMA
      IN = INT(I/2) + MOD(I,2)
      NI = NSET(2, IN+1)
      IEQ = 2*NI - MOD(I,2) 
c     IEQ=IEQ+1
      F(IEQ)=FB1(I)
370   CONTINUE
c------     26/2/98 ----------
      END IF
C     REWIND MT3
C     WRITE(MT3) F
C     REWIND MT11
C     WRITE(MT11) ((HB1(I,J),I=1,N2GAM1),J=1,NBAND1)
101   CONTINUE
      IFLSTC = 0
      RETURN
      END
C
      SUBROUTINE CORELN(COORL,RAD1,RAD2,THT1,THT2,RNTCH,DTHT,N,NDM,
     1 NELGAM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C---- FORMS COORDS OF NODES OF ELEMENT N IN THE RING
C
      DIMENSION COORL(NDM,*)
      DATA PI/3.141592653589793/
c---- following two statements are for mixed-mode case (full ring)
      THT1=DFLOAT((NELGAM-N))*DTHT - PI
      THT2=DFLOAT((NELGAM-N+1))*DTHT - PI
c---- following two statements are for mode I symm case (1/2 ring)
c     tht1 = dfloat((nelgam-n))*dtht
c     tht2 = dfloat((nelgam-n+1))*dtht
c
      COORL(1,1)=RAD2*DCOS(THT1)
      COORL(2,1)=RAD2*DSIN(THT1)
      COORL(1,2)=RAD2*DCOS(THT2)
      COORL(2,2)=RAD2*DSIN(THT2)
      COORL(1,3)=RAD1*DCOS(THT2)
      COORL(2,3)=RAD1*DSIN(THT2)
      COORL(1,4)=RAD1*DCOS(THT1)
      COORL(2,4)=RAD1*DSIN(THT1)
c---- following statements are for notch - full ring; comment out for sharp crack.
      IF(N.EQ.1) THEN
      COORL(2,2) = RNTCH
      COORL(2,3) = RNTCH
      END IF
      IF(N.EQ.NELGAM) THEN
      COORL(2,1) = -RNTCH
      COORL(2,4) = -RNTCH
      END IF
c
      RETURN
      END
C
      SUBROUTINE BSOLV1(A,B,N,NBD,N1)
C     PERFORMS PARTIAL FORWARD REDN OF AX=B WHERE A IS A SYMMETRIC
C     MATRIX IN BANDED FORM OF SIZE N*NBD AND B IS A VECTOR OF SIZE N*1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,*),B(*)
C *** PARTIAL FORWARD REDUCTION OF A ***
      AVOD=0.D0
      DO 10 I=1,N
C     AVOD=AVOD+ABS(A(I,1))
      AVOD=AVOD+DABS(A(I,1))
10    CONTINUE
C     AVOD=AVOD/FLOAT(N)
      AVOD=AVOD/DFLOAT(N)
C     TOL=AVOD*1.0E-06
      TOL=AVOD*1.0D-06
      DMIN=1.0D+30
      DO 100 I=1,N1
      IROW=I
      IF(A(I,1).LT.DMIN) DMIN=A(I,1)
      IF(DMIN.LT.TOL)GO TO 149
      MIN=MIN0(NBD,N-I+1)
      DO 105 J=2,MIN
      ST=A(I,J)/A(I,1)
      IJ1=I+J-1
      DO 110 K=J,MIN
      KJ1=K-J+1
  110 A(IJ1,KJ1)=A(IJ1,KJ1)-ST*A(I,K)
  105 A(I,J)=ST
  100 CONTINUE
C *** PARTIAL FORWARD REDUCTION OF B
      DO 120 I=1,N1
      MIN=MIN0(NBD,N-I+1)
      DO 125 J=2,MIN
      JJ=I+J-1
  125 B(JJ)=B(JJ)-B(I)*A(I,J)
  120 CONTINUE
      RETURN
  149 WRITE(6,1002)
 1002 FORMAT(//,1X,'AVOD,DMIN,TOL,IROW (DMIN.LT.TOL IN PARTIAL FORW
     1 ARD REDUCTION)')
      WRITE(6,1003) AVOD,DMIN,TOL,IROW
 1003 FORMAT(1X,3E13.5,I4)
      STOP
      END
C
      SUBROUTINE BSOLV2(A,B,N,NBD,N1)
C     PERFORMS PARTIAL FORWARD REDN OF AX=B WHERE A IS A SYMMETRIC
C     MATRIX IN BANDED FORM OF SIZE N*NBD AND B IS A VECTOR OF SIZE N*1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,*),B(*)
C *** PARTIAL FORWARD REDUCTION OF A ***
      AVOD=0.D0
      DO 10 I=1,N
C     AVOD=AVOD+ABS(A(I,1))
      AVOD=AVOD+DABS(A(I,1))
10    CONTINUE
C     AVOD=AVOD/FLOAT(N)
      AVOD=AVOD/DFLOAT(N)
      TOL=AVOD*1.0D-06
      DMIN=1.0D+30
      DO 100 I=1,N1
      IROW=I
      IF(A(I,1).LT.DMIN) DMIN=A(I,1)
      IF(DMIN.LT.TOL)GO TO 149
      MIN=MIN0(NBD,N-I+1)
      DO 105 J=2,MIN
      ST=A(I,J)/A(I,1)
      IJ1=I+J-1
      DO 110 K=J,MIN
      KJ1=K-J+1
  110 A(IJ1,KJ1)=A(IJ1,KJ1)-ST*A(I,K)
  105 A(I,J)=ST
  100 CONTINUE
      WRITE(6,1000)
1000  FORMAT(//,1X,'AVOD,DMIN')
      WRITE(6,1001) AVOD,DMIN
1001  FORMAT(1X,2E13.5)
C *** PARTIAL FORWARD REDUCTION OF B
      DO 120 I=1,N1
      MIN=MIN0(NBD,N-I+1)
      DO 125 J=2,MIN
      JJ=I+J-1
  125 B(JJ)=B(JJ)-B(I)*A(I,J)
  120 CONTINUE
      RETURN
  149 WRITE(6,1002)
 1002 FORMAT(//,1X,'AVOD,DMIN,TOL,IROW (DMIN.LT.TOL IN PARTIAL FORW
     1 ARD REDUCTION)')
      WRITE(6,1003) AVOD,DMIN,TOL,IROW
 1003 FORMAT(1X,3E13.5,I4)
      STOP
      END
C
C
      SUBROUTINE ASSMBL(H,F,HE,FE,IDL,NUMEQ,NBAND,NED,M1,M2)
C---- ASSEMBLES ELEM STIFF MATRIX HE INTO GLOBAL (BANDED, SYMM) MATRIX
C---- H IF M1.NE.0 AND FORCE VECTOR FE INTO GLOBAL FORCE VECTOR F IF M2.NE.0.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(NUMEQ,*),F(*),HE(NED,*),FE(*),IDL(*)
      IF(M1.EQ.0) GO TO 2
C---- ASSEMBLE HE INTO H.
      DO 10 I=1,NED
      II=IDL(I)
      IF(II.EQ.0) go to 10
      if(ii.lt.0) then
c---- specified displacement dof ; to avoid zero pivot set entry in 1st col as 1
      H(-II,1) = 1.d0
      GO TO 10
      end if
c
      DO 20 J=1,NED
      IJ=IDL(J)-II+1
      IF(IJ.GT.0) H(II,IJ)=H(II,IJ)+HE(I,J)
   20 CONTINUE
   10 CONTINUE
    2 IF(M2.EQ.0) RETURN
C---- ASSEMBLE FE INTO F.
      DO 30 I=1,NED
      II=IDL(I)
      IF(II.EQ.0) GO TO 30
      F(II)=F(II)+FE(I)
   30 CONTINUE
      RETURN
      END
c
      SUBROUTINE HADD(A,C,HB1,JDIAG,N2GAMA,N4GAMA,NEQ,AFL,CFL)
c
c---- modified on 26/2/98 to treat arbitrary numbering of outer
c---- ring of nodes. Note NSET(2,I+1), I=1,ngama contains the
c---- the global node numbers of the outer ring of nodes in clockwise
c---- sense starting from upper crack flank and ending in lower crack
c---- flank.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---- ADDS STIFFNESS FROM FAR-FIELD CONDENSED ELASTIC REGION TO NODES
C---- ALONG OUTER CURVE OF ACTIVE REGION (GAMA CURVE).
      DIMENSION A(*),C(*),HB1(N4GAMA,*),JDIAG(*)
      COMMON /PRTSET/ NSET(9,201)
      LOGICAL AFL,CFL
      DO 100 I=1,N2GAMA
      IN = INT(I/2) + MOD(I,2)
      NI = NSET(2, IN+1)
      NID = 2*NI - MOD(I,2)
      NJE = N2GAMA-I+1
      DO 100 J=1,NJE
      J1 = I+J-1
      JN = INT(J1/2) + MOD(J1,2)
      NJ  = NSET(2, JN+1)
      NJD = 2*NJ - MOD(J1,2)
      IF(NJD.GE.NID) THEN
      M = JDIAG(NJD) - NJD + NID
      ELSE
      M = JDIAG(NID) - NID + NJD
      END IF
      IF(AFL) A(M) = A(M) + HB1(I,J)
      IF(CFL) C(M) = C(M) + HB1(I,J)
100   CONTINUE
      RETURN
      END
C
      SUBROUTINE FADD(HB1,FB1,PROP,DR,UTOT,N2GAMA,N4GAMA,NUMEQ)
c---- modified on 26/2/98
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION HB1(N4GAMA,*),DR(*),UTOT(*),FB1(*),PROP(7)
      COMMON /PRTSET/ NSET(9,201)
C *** SUBTRACT CONTRIBUTION DUE TO ADDED STIFFNESS IN RING OF ELEMENTS
C     ALONG GAMA CURVE.
      DO 200 I=1,N2GAMA
      IN = INT(I/2) + MOD(I,2)
      NI = NSET(2,IN+1)
      II = 2*NI - MOD(I,2)
      I1=I-1
      DO 210 J=I,N2GAMA
      JN = INT(J/2) + MOD(J,2)
      NJ = NSET(2, JN+1)
      JJ = 2*NJ - MOD(J,2)
      JI1=J-I1
      DR(II)=DR(II)-HB1(I,JI1)*UTOT(JJ)
C *** HB1UNBAND(I,J)=HB1BAND(I,J-I+1)
210   CONTINUE
      IF(I.EQ.1) GO TO 200
      DO 220 J=1,I1
      JN = INT(J/2) + MOD(J,2)
      NJ = NSET(2, JN+1)
      JJ = 2*NJ - MOD(J,2)
      IJ1=I-J+1
      DR(II)=DR(II)-HB1(J,IJ1)*UTOT(JJ)
C *** HB1UNBAND(I,J)=HB1UNBAND(J,I)=HB1BAND(J,I-J+1)
220   CONTINUE
200   CONTINUE
C----  following are for T-stress terms (comment out if not used)
c----  prop(2) = scaling factor for T-stress.
c----  prop(1) = scaling factor for K-field terms
c----  FB1 will contain the condensed forced vector of outer elastic
c----  region due to prescription of T-terms on remote boundary
c     IEQ=NUMEQ-N2GAMA
c     DO 230 I=1,N2GAMA
c     IEQ=IEQ+1
c     DR(IEQ)=DR(IEQ)+PROP(2)*FB1(I)
c230   CONTINUE
      DO 230 I=1,N2GAMA
      IN = INT(I/2) + MOD(I,2)
      NI = NSET(2, IN+1)
      IEQ = 2*NI - MOD(I,2)
      DR(IEQ)= DR(IEQ)+PROP(2)*FB1(I)
230   CONTINUE
c------------------------------------------
      RETURN
      END
C
      SUBROUTINE HADD1(A,C,HB1,JDIAG,N2GAMA,N4GAMA,NEQ,AFL,CFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---- ADDS STIFFNESS FROM FAR-FIELD CONDENSED ELASTIC REGION TO NODES
C---- ALONG OUTER CURVE OF ACTIVE REGION (GAMA CURVE).
      DIMENSION A(*),C(*),HB1(N4GAMA,*),JDIAG(*)
      LOGICAL AFL,CFL
      NMARK = NEQ-N2GAMA
      DO 100 I=1,N2GAMA
      NMI = NMARK+I
      NJ = N2GAMA-I+1
      DO 100 J=1,NJ
      NMJ = NMI+J-1
      M=JDIAG(NMJ)-NMJ+NMI
      IF(AFL) A(M) = A(M) + HB1(I,J)
      IF(CFL) C(M) = C(M) + HB1(I,J)
100   CONTINUE
      RETURN
      END
C
      SUBROUTINE FADD1(HB1,FB1,PROP,DR,UTOT,N2GAMA,N4GAMA,NUMEQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION HB1(N4GAMA,*),DR(*),UTOT(*),FB1(*),PROP(7)
C *** SUBTRACT CONTRIBUTION DUE TO ADDED STIFFNESS IN RING OF ELEMENTS
C     ALONG GAMA CURVE.
      II1=NUMEQ-N2GAMA
      JJ1=NUMEQ-N2GAMA
      DO 200 I=1,N2GAMA
      II=II1+I
      I1=I-1
      DO 210 J=I,N2GAMA
      JJ=J+JJ1
      JI1=J-I1
      DR(II)=DR(II)-HB1(I,JI1)*UTOT(JJ)
C *** HB1UNBAND(I,J)=HB1BAND(I,J-I+1)
210   CONTINUE
      IF(I.EQ.1) GO TO 200
      DO 220 J=1,I1
      JJ=JJ1+J
      IJ1=I-J+1
      DR(II)=DR(II)-HB1(J,IJ1)*UTOT(JJ)
C *** HB1UNBAND(I,J)=HB1UNBAND(J,I)=HB1BAND(J,I-J+1)
220   CONTINUE
200   CONTINUE
      IEQ=NUMEQ-N2GAMA
      DO 230 I=1,N2GAMA
      IEQ=IEQ+1
      DR(IEQ)=DR(IEQ)+PROP(2)*FB1(I)
230   CONTINUE
      RETURN
      END
C
      SUBROUTINE PROFL1(JDIAG,ID,IX,NDF,NEN1,NAD,IOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- COMPUTE PROFILE OF GLOBAL ARRAYS TAKING INTO ACCOUNT LAST
C---- RING OF NODES WHICH MAY BE RANDOMLY NUMBERED
C
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /PRTSET/ NSET(9,201)
      DIMENSION JDIAG(*),ID(NDF,*),IX(NEN1,*)
      DATA WNOD/4HNODE/,WDIA/4HDIAG/,WDOF/4H DOF/,WEQN/4H EQN/
      GO TO (1,2), IOP
C---- SET UP SYSTEM SIZE
1     NEQ = NUMNP*NDF
      RETURN
C---- COMPUTE COLUMN HEIGHTS
2     CONTINUE
c---- first initialize jdiag array
      do 50 i=1,neq
      jdiag(i) = 0
50    continue
c
      DO 500 N = 1,NUMEL
      DO 400 I = 1,NEN
      II = IX(I,N)
      IF(II.EQ.0) GO TO 400
      DO 300 K = 1,NDF
      KK = (II-1)*NDF + K
      DO 200 J = I,NEN
      JJ = IX(J,N)
      IF(JJ.EQ.0) GO TO 200
      DO 100 L = 1,NDF
      LL = (JJ-1)*NDF + L
      M = MAX0(KK,LL)
      JDIAG(M) = MAX0(JDIAG(M),IABS(KK-LL))
100   CONTINUE
200   CONTINUE
300   CONTINUE
400   CONTINUE
500   CONTINUE
c---- scan thro' set of node nos corresp to outer boundary and identify
c---- the minimum one
      minnd = 500000
      ngama = nset(2,1)
      if(ngama.le.0) then
      write(6,2006) 
      stop
      end if
      do 510 i=1,ngama
      if(nset(2,i+1).lt.minnd) minnd = nset(2,i+1)
510   continue
      write(6,2005) minnd
      mindof = (minnd-1)*ndf + 1
c---- now raise all columns corresp to dofs on outer boundary at least
c---- up to mindof
      do 515 i=1,ngama
      j = nset(2,i+1)
      do 520 k=1,ndf
      jj = (j-1)*ndf + k
      jjm = jj - mindof
      if(jj.gt.mindof) jdiag(jj) = max0(jdiag(jj),jjm)
520   continue
515   continue
c
C---- COMPUTE DIAGONAL POINTERS FOR PROFILE
      NAD = 1
      JDIAG(1) = 1
      IF(NEQ.EQ.1) RETURN
      DO 600 N = 2,NEQ
600   JDIAG(N) = JDIAG(N) + JDIAG(N-1) + 1
      NAD = JDIAG(NEQ)
      WRITE(6,2004) NAD
      RETURN
2000  FORMAT(10(4X,A4))
2001  FORMAT(10I8)
2002  FORMAT(5X,6(1X,A4,1X,A4))
2003  FORMAT((5X,6(1X,I3,1H=,I5)))
2004  FORMAT(//5X,'SIZE OF REVISED STIFFNESS MATRIX IN PROFILE FORM='
     1 ,I8)
2005  FORMAT(//5X,'MIN NODE NUMBER ON OUTER BDRY =',I8)
2006  FORMAT(//5X,'***ERROR *** NUMBER OF NODES ON OUTER BDRY INPUT',
     1 ' AS LESS THAN EQ 0 IN SET 2')
      END
