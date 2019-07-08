      SUBROUTINE PFORM_SWAPNIL(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,
     1 IDPROP,PROP,F,T, JDIAG,STR,EPS,Q,B,A,C,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,ISW,U,UD,AFL,BFL,CFL,DFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!---- COMPUTE ELEMENT ARRAYS AND ASSEMBLE GLOBAL ARRAYS

      LOGICAL AFL,BFL,CFL,DFL,FLAG
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
!---- COMMON/PRLOD/ PROP
      COMMON /CNVG1/ICONV,IDIVG,IPATH
      COMMON /STEADY/ NELROW,IELROW,NELOUT
      DIMENSION XL(NDM,NEN),LD(NDF,NEN),P(NST),S(NST,*),IE(7,*),D(40,*),
     1 ID(NDF,*),X(NDM,*),IX(NEN1,*),F(NDF,*),JDIAG(*),B(*),A(*),C(*),
     2 UL(NDF,2*NEN),TL(NEN),T(*),U(NDF,*),UD(NDF,*)
      DIMENSION PROP(7),IDPROP(NDF,*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(NDF,*),ACCELG(NDF,*),VEL(NDF,NEN),ACCEL(NDF,NEN)

!---- LOOP ON ELEMENTS
      IEL = 0
      IELROW = NELOUT + 1

c*********************************************************************************************************************************
!$ACC REGION

      DO 110 N = 1,NEQ
!---- SET UP LOCAL ARRAYS
      MA = IX(NEN1,N)
      FLAG = .FALSE.
      DO 108 I = 1,NEN
      II = IX(I,N)
      IF(II.NE.0) GO TO 105
      TL(I) = 0.D0
      DO 103 J = 1,NDM
103   XL(J,I) = 0.D0
      DO 104 J = 1,NDF
      VEL(J,I) = 0.D0
      ACCEL(J,I) = 0.D0
      UL(J,I) = 0.D0
      UL(J,I+NEN)=0.D0
104   LD(J,I) = 0
      GO TO 108
105   CONTINUE
      NEL = I
      TL(I) = T(II)
      DO 106 J = 1,NDM
106   XL(J,I) = X(J,II)
      DO 107 J = 1,NDF
      JJ = IE(J,MA)
      IF (ID(JJ,II).NE.0) FLAG = .TRUE.
      VEL(J,I) = VELG(JJ,II)
      ACCEL(J,I) = ACCELG(JJ,II)
      IF(JJ.LE.0) GO TO 107
      UL(J,I) = U(JJ,II)
      UL(J,I+NEN)=UD(JJ,II)
      LD(J,I) = (II-1)*NDF + JJ
107   CONTINUE
108   CONTINUE
C     IF (ISW.EQ.5) FLAG = .FALSE.
!---- FORM ELEMENT ARRAY
      IF(IE(7,MA).NE.IEL) MCT = 0
      IEL = IE(7,MA)
!------------------------------------------------------------------------
!---- INTEGRATION OF CONSTITUTIVE EQUATIONS FOR STEADY STATE CRACK GROWTH
      IF(ISW.EQ.9) THEN
      IF(NELROW.GT.0) THEN
      IF((N.LE.NELOUT).OR.(N.EQ.IELROW)) THEN
      DO 111 I=1,NSTR
      STR(I,N) = 0.D0
      EPS(I,N) = 0.D0
111   CONTINUE
!---- DO 112 I=1,NQ
!112- Q(I,N) = 0.D0
      END IF
      IF((N.GT.NELOUT).AND.(N.NE.IELROW)) THEN
      DO 114 I =1,NSTR
      STR(I,N) = STR(I,N-1)
      EPS(I,N) = EPS(I,N-1)
114   CONTINUE
!---- IQUAD = 1
!---- IQDM = (IQUAD-1)*NQDM
      DO 115 I = 1,NQ
!
!---- QL(10,IQUAD) AND QL(11,IQUAD) WHICH CONTAIN THE VALUES OF DELTA X
!---- AND DELTA WP/DELTA X1 AT GAUSS POINTS IQUAD MUST NOT BE ALTERED.
!---- IF(I.EQ.(IQDM+10)) GO TO 115
!---- IF(I.EQ.(IQDM+11)) THEN
!---- IQUAD = IQUAD + 1
!---- IQDM = (IQUAD-1)*NQDM
!---- GO TO 115
!---- END IF

      Q(I,N) = Q(I,N-1)
115   CONTINUE
      END IF
      END IF
      END IF
!--------------------------------------------------------------------------
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,ISW)

      IF(ISW.EQ.9) THEN
      IF(NELROW.GT.0) THEN
      IF(N.EQ.IELROW) IELROW = IELROW+NELROW
      END IF
      END IF

110   CONTINUE

!$ACC END REGION
!********************************************************************************************************************************
      RETURN
      END


      SUBROUTINE PFORM22(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,
     1   JDIAG,STR,EPS,Q,B,A,C,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,ISW,U,UD,
     3   AFL,BFL,CFL,DFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTE ELEMENT ARRAYS AND ASSEMBLE GLOBAL ARRAYS
!---- (OLD VERSION WITHOUT CONSISTENT IMPLEMENTATION OF
!---- DISPLACEMENT BOUNDARY CONDITIONS)

      LOGICAL AFL,BFL,CFL,DFL
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
      DIMENSION XL(NDM,NEN),LD(NDF,NEN),P(NST),S(NST,NST),IE(7,*),
     1  ID(NDF,*),X(NDM,*),IX(NEN1,*),F(NDF,*),JDIAG(*),B(*),A(*),C(*),
     2  UL(NDF,2*NEN),TL(NEN),T(*),U(*),UD(NDF,*),D(40,*)
      DIMENSION IDPROP(NDF,*),PROP(7)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(NDF,*),ACCELG(NDF,*),VEL(NDF,NEN),ACCEL(NDF,NEN)
!---- LOOP ON ELEMENTS
      IEL = 0

!*********************************************************************************************************************************
!$ACC REGION

      DO 110 N = 1,NEQ

!---- SET UP LOCAL ARRAYS
      MA = IX(NEN1,N)
      DO 108 I = 1,NEN
      II = IX(I,N)
      IF(II.NE.0) GO TO 105
      TL(I) = 0.d0
      DO 103 J = 1,NDM
103   XL(J,I) = 0.d0
      DO 104 J = 1,NDF
      VEL(J,I) = 0.D0
      ACCEL(J,I) = 0.D0
      UL(J,I) = 0.D0
      UL(J,I+NEN)=0.D0
104   LD(J,I) = 0
      GO TO 108
105   IID = II*NDF - NDF
      NEL = I
      TL(I) = T(II)
      DO 106 J = 1,NDM
106   XL(J,I) = X(J,II)
      DO 107 J = 1,NDF
      JJ = IE(J,MA)
      VEL(J,I) = VELG(JJ,II)
      ACCEL(J,I) = ACCELG(JJ,II)
      IF(JJ.LE.0) GO TO 107
!     K = IABS(ID(JJ,II))
      K = (II-1)*NDF + JJ
      IND=IDPROP(JJ,II)
      UL(J,I) = F(JJ,II)*PROP(IND)
      UL(J,I+NEN)=UD(JJ,II)
      IF(K.GT.0) UL(J,I) = U(K)
      IF(DFL) K = IID + JJ
107   LD(J,I) = K
108   CONTINUE
!---- FORM ELEMENT ARRAY
      IF(IE(7,MA).NE.IEL) MCT = 0
      IEL = IE(7,MA)
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,ISW)
!---- ADD TO TOTAL ARRAY
      IF(AFL.OR.BFL.OR.CFL) then 
      CALL ADDSTF_SWAPNIL(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      end if
110   CONTINUE

!$ACC END REGION
!********************************************************************************************************************************
      RETURN
      END


      SUBROUTINE ADDSTF_SWAPNIL(A,B,C,S,P,JDIAG,LD,NST,NEL,AFL,BFL,CFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ASSEMBLE GLOBAL ARRAYS

      LOGICAL AFL,BFL,CFL
      DIMENSION A(*),B(*),JDIAG(*),P(*),S(NST,*),LD(*) ,C(*)
      DO 200 J = 1,NEL
      K = LD(J)
      IF(K.EQ.0) GO TO 200
      IF(BFL) then	 
      B(K) = B(K) + P(J)
	  end if 
      IF(.NOT.AFL.AND..NOT.CFL) GO TO 200
      L = JDIAG(K) - K
      DO 100 I = 1,NEL
      M = LD(I)
      IF(M.GT.K.OR.M.EQ.0) GO TO 100
      M = L + M
      IF(AFL) then	 
      A(M) = A(M) + S(I,J)
	  end if 
      IF(CFL) then	 
      C(M) = C(M) + S(J,I)
	  end if
100   CONTINUE
200   CONTINUE
      RETURN
      END

!********************************************************************************************************************************
!********************************************************************************************************************************
!********************************************************************************************************************************


      SUBROUTINE PFORM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,
     1   JDIAG,STR,EPS,Q,B,A,C,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,ISW,U,UD,AFL,BFL,CFL,DFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTE ELEMENT ARRAYS AND ASSEMBLE GLOBAL ARRAYS

      LOGICAL AFL,BFL,CFL,DFL,FLAG
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
!---- COMMON/PRLOD/ PROP
      COMMON /CNVG1/ICONV,IDIVG,IPATH
      COMMON /STEADY/ NELROW,IELROW,NELOUT
      DIMENSION XL(NDM,*),LD(NDF,*),P(*),S(NST,*),IE(7,*),D(40,*),
     1 ID(NDF,*),X(NDM,*),IX(NEN1,*),F(NDF,*),JDIAG(*),B(*),A(*),C(*),
     2 UL(NDF,*),TL(*),T(*),U(NDF,*),UD(NDF,*)
      DIMENSION PROP(7),IDPROP(NDF,*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(NDF,*),ACCELG(NDF,*),VEL(NDF,*),ACCEL(NDF,*)
!---- LOOP ON ELEMENTS
      IEL = 0
      IELROW = NELOUT + 1
      DO 110 N = 1,NUMEL
!---- SET UP LOCAL ARRAYS
      MA = IX(NEN1,N)
      FLAG = .FALSE.
      DO 108 I = 1,NEN
      II = IX(I,N)
      IF(II.NE.0) GO TO 105
      TL(I) = 0.D0
      DO 103 J = 1,NDM
103   XL(J,I) = 0.D0
      DO 104 J = 1,NDF
      VEL(J,I) = 0.D0
      ACCEL(J,I) = 0.D0
      UL(J,I) = 0.D0
      UL(J,I+NEN)=0.D0
104   LD(J,I) = 0
      GO TO 108
105   CONTINUE
      NEL = I
      TL(I) = T(II)
      DO 106 J = 1,NDM
106   XL(J,I) = X(J,II)
      DO 107 J = 1,NDF
      JJ = IE(J,MA)
      IF (ID(JJ,II).NE.0) FLAG = .TRUE.
      VEL(J,I) = VELG(JJ,II)
      ACCEL(J,I) = ACCELG(JJ,II)
      IF(JJ.LE.0) GO TO 107
      UL(J,I) = U(JJ,II)
      UL(J,I+NEN)=UD(JJ,II)
      LD(J,I) = (II-1)*NDF + JJ
107   CONTINUE
108   CONTINUE
!---- IF (ISW.EQ.5) FLAG = .FALSE.
!---- FORM ELEMENT ARRAY
      IF(IE(7,MA).NE.IEL) MCT = 0
      IEL = IE(7,MA)
!------------------------------------------------------------------------
!---- INTEGRATION OF CONSTITUTIVE EQUATIONS FOR STEADY STATE CRACK GROWTH
      IF(ISW.EQ.9) THEN
      IF(NELROW.GT.0) THEN
      IF((N.LE.NELOUT).OR.(N.EQ.IELROW)) THEN
      DO 111 I=1,NSTR
      STR(I,N) = 0.D0
      EPS(I,N) = 0.D0
111   CONTINUE
!---- DO 112 I=1,NQ
!112- Q(I,N) = 0.D0
      END IF
      IF((N.GT.NELOUT).AND.(N.NE.IELROW)) THEN
      DO 114 I =1,NSTR
      STR(I,N) = STR(I,N-1)
      EPS(I,N) = EPS(I,N-1)
114   CONTINUE
!     IQUAD = 1
!     IQDM = (IQUAD-1)*NQDM
      DO 115 I = 1,NQ

!---- QL(10,IQUAD) AND QL(11,IQUAD) WHICH CONTAIN THE VALUES OF DELTA X
!---- AND DELTA WP/DELTA X1 AT GAUSS POINTS IQUAD MUST NOT BE ALTERED.
!     IF(I.EQ.(IQDM+10)) GO TO 115
!     IF(I.EQ.(IQDM+11)) THEN
!     IQUAD = IQUAD + 1
!     IQDM = (IQUAD-1)*NQDM
!     GO TO 115
!     END IF

      Q(I,N) = Q(I,N-1)
115   CONTINUE
      END IF
      END IF
      END IF
!--------------------------------------------------------------------------
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,ISW)
!
      IF(ISW.EQ.9) THEN
      IF(NELROW.GT.0) THEN
      IF(N.EQ.IELROW) IELROW = IELROW+NELROW
      END IF
      END IF

!---- CHECK FOR PLASTIC DISSIPATION AND TOO LARGE PLASTIC STRAIN INCREMENTS
      IF(IDIVG.EQ.999) RETURN

      IF ((.NOT.FLAG).AND.(AFL.OR.BFL.OR.CFL))
     1 CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      IF (FLAG.AND.(AFL.OR.CFL)) THEN
      CALL MODIFY(S,P,IDPROP,PROP,F,U,ID,IX,N,NST,NEN,NDF,NEN1,AFL,
     1 BFL,CFL)
      CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      END IF
      IF (FLAG.AND.BFL.AND.(.NOT.AFL).AND.(.NOT.CFL)) THEN
      CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      IF (ISW.NE.5) THEN
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,3)
      DO 120 I = 1,NST
120   P(I) = 0.D0
      CALL MODIFY(S,P,IDPROP,PROP,F,U,ID,IX,N,NST,NEN,NDF,NEN1,AFL,
     1 BFL,CFL)
      CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      END IF
      END IF
110   CONTINUE
      IF (ISW.EQ.5) RETURN
      IF(AFL.OR.BFL.OR.CFL)
     1 CALL DIAG(A,B,C,IDPROP,PROP,F,U,JDIAG,ID,NDF,NUMNP,AFL,BFL,CFL)
      RETURN
      END


      SUBROUTINE MODIFY(S,P,IDPROP,PROP,F,U,ID,IX,N,NST,NEN,NDF,NEN1,
     1 AFL,BFL,CFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- MODIFY LOCAL STIFFNESS AND FORCE ARRAYS TO ACCOUNT FOR BC

      LOGICAL AFL,BFL,CFL
!     COMMON/PRLOD/ PROP
      DIMENSION S(NST,*),P(*),ID(NDF,*),IX(NEN1,*),F(NDF,*),U(NDF,*)
      DIMENSION IDPROP(NDF,*),PROP(7)

      IF (AFL.OR.CFL) THEN
      DO 10 J = 1,NEN
      NODE = IX(J,N)
      IF (NODE.EQ.0) GO TO 10
      DO 20 I = 1,NDF
      IS = (J-1)*NDF + I
      IF (ID(I,NODE).NE.0) THEN
      DO 30 JS = 1,NST
      S(JS,IS) = 0.D0
   30 S(IS,JS) = 0.D0
      END IF
   20 CONTINUE
   10 CONTINUE
      END IF
      IF (BFL) THEN
      DO 35 K = 1,NEN
      KL = IX(K,N)
      IF (KL.NE.0) THEN
      DO 40 L = 1,NDF
      IF (ID(L,KL).EQ.0) THEN
      IS = (K-1)*NDF + L
      DO 50 J = 1,NEN
      IJ = IX(J,N)
      IF (IJ.NE.0) THEN
      DO 60 I = 1,NDF
      IF (ID(I,IJ).NE.0) THEN
      JS = (J-1)*NDF + I
      IND=IDPROP(I,IJ)
      P(IS) = P(IS) - S(IS,JS)*(PROP(IND)*F(I,IJ) - U(I,IJ))
      END IF
   60 CONTINUE
      END IF
   50 CONTINUE
      END IF
   40 CONTINUE
      END IF
   35 CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE PFORM2(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,
     1   JDIAG,STR,EPS,Q,B,A,C,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,ISW,U,UD,
     3   AFL,BFL,CFL,DFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTE ELEMENT ARRAYS AND ASSEMBLE GLOBAL ARRAYS
!---- (OLD VERSION WITHOUT CONSISTENT IMPLEMENTATION OF
!---- DISPLACEMENT BOUNDARY CONDITIONS)

      LOGICAL AFL,BFL,CFL,DFL
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
!     COMMON/PRLOD/ PROP
      DIMENSION XL(NDM,*),LD(NDF,*),P(*),S(NST,*),IE(7,*),D(40,*),
     1  ID(NDF,*),X(NDM,*),IX(NEN1,*),F(NDF,*),JDIAG(*),B(*),A(*),C(*),
     2   UL(NDF,*),TL(*),T(*),U(*),UD(NDF,*)
      DIMENSION IDPROP(NDF,*),PROP(7)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(NDF,*),ACCELG(NDF,*),VEL(NDF,*),ACCEL(NDF,*)
!---- LOOP ON ELEMENTS
      IEL = 0
      DO 110 N = 1,NUMEL
!---- SET UP LOCAL ARRAYS
      MA = IX(NEN1,N)
      DO 108 I = 1,NEN
      II = IX(I,N)
      IF(II.NE.0) GO TO 105
      TL(I) = 0.
      DO 103 J = 1,NDM
103   XL(J,I) = 0.
      DO 104 J = 1,NDF
      VEL(J,I) = 0.D0
      ACCEL(J,I) = 0.D0
      UL(J,I) = 0.D0
      UL(J,I+NEN)=0.D0
104   LD(J,I) = 0
      GO TO 108
105   IID = II*NDF - NDF
      NEL = I
      TL(I) = T(II)
      DO 106 J = 1,NDM
106   XL(J,I) = X(J,II)
      DO 107 J = 1,NDF
      JJ = IE(J,MA)
      VEL(J,I) = VELG(JJ,II)
      ACCEL(J,I) = ACCELG(JJ,II)
      IF(JJ.LE.0) GO TO 107
!     K = IABS(ID(JJ,II))
      K = (II-1)*NDF + JJ
      IND=IDPROP(JJ,II)
      UL(J,I) = F(JJ,II)*PROP(IND)
      UL(J,I+NEN)=UD(JJ,II)
      IF(K.GT.0) UL(J,I) = U(K)
      IF(DFL) K = IID + JJ
107   LD(J,I) = K
108   CONTINUE
!---- FORM ELEMENT ARRAY
      IF(IE(7,MA).NE.IEL) MCT = 0
      IEL = IE(7,MA)
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,ISW)
!---- ADD TO TOTAL ARRAY
      IF(AFL.OR.BFL.OR.CFL) CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,
     1   AFL,BFL,CFL)
110   CONTINUE
      RETURN
      END

      SUBROUTINE DIAG(A,B,C,IDPROP,PROP,F,U,JDIAG,ID,NDF,NUMNP,AFL,
     1 BFL,CFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ADJUST DIAGONAL OF STIFFNESS MATRIX AND RHS TO ACCOUNT FOR BC

      LOGICAL AFL,BFL,CFL
!     COMMON/PRLOD/ PROP
      DIMENSION A(*),B(*),C(*),JDIAG(*),ID(NDF,*),F(NDF,*),U(NDF,*)
      DIMENSION IDPROP(NDF,*),PROP(7)

      IDF = 0
      DO 10 J = 1,NUMNP
      DO 10 I = 1,NDF
      IDF = IDF + 1
      IF (ID(I,J).EQ.0) GO TO 10
      IF (AFL) A(JDIAG(IDF)) = 1.
      IF (CFL) C(JDIAG(IDF)) = 1.
      IF (BFL) THEN
      IND=IDPROP(I,J)
      B(IDF)=PROP(IND)*F(I,J)-U(I,J)
      END IF
!     IF (BFL) B(IDF) = PROP*F(I,J) - U(I,J)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE RESETQ(Q,QNEW,IQ,NUMEL,NQUAD,NQDM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- RESET IQTH ENTRY OF INTERNAL VARIABLE ARRAY TO QNEW

      DIMENSION Q(*)
      IF ((IQ.LE.0).OR.(IQ.GT.NQDM)) THEN
      WRITE(6,'(///5X,''** FATAL ERROR 07** ILLEGAL IQ ='',
     1 I5,'' IN MACRO RESET''///)') IQ
      STOP
      END IF
      II = IQ
      DO 10 I = 1,NUMEL
      DO 10 J = 1,NQUAD
      Q(II) = QNEW
   10 II = II + NQDM
      RETURN
      END
      SUBROUTINE ACTCOL(A,B,JDIAG,NEQ,AFAC,BACK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL AFAC,BACK,PRF
      COMMON/ENGYS/ AENGY
      COMMON/STURM/NEIGLO
      COMMON /PRINT/ PRF
      DIMENSION A(*),B(*),JDIAG(*)
      DATA TOL/1.D-15/

!---- ACTIVE COLUMN PROFILE SYMMETRIC EQUATION SOLVER

!---- FACTOR A TO UT*D*U, REDUCE B
      AENGY = 0.D0
      NEIGLO = 0
      JR = 0
      DO 600 J = 1,NEQ
      JD = JDIAG(J)
      JH = JD - JR
      IS = J - JH + 2
      DG = A(JD)
!     IF(JH-2) 550,300,100
      IF(JH-2) 501,300,100
100   IF(.NOT.AFAC) GO TO 500
      IE = J - 1
      K = JR + 2
      ID = JDIAG(IS - 1)
!---- REDUCE ALL EQUATIONS EXCEPT DIAGONAL
      DO 200 I = IS,IE
      IR = ID
      ID = JDIAG(I)
      IH = MIN0(ID-IR-1,I-IS+1)
      IF(IH.GT.0) A(K) = A(K) - DOT(A(K-IH),A(ID-IH),IH)
200   K = K + 1
!---- REDUCE DIAGONAL TERM
300   IF(.NOT.AFAC) GO TO 500
      IR = JR + 1
      IE = JD - 1
      K = J - JD
      DO 400 I = IR,IE
      ID = JDIAG(K+I)
      D = A(I)
      A(I) = A(I)*A(ID)
      DG = DG - D*A(I)
400   CONTINUE
501   CONTINUE
      IF(DG*A(JD).LT.0.0D0) WRITE(6,2000) J
      IF(DABS(DG).LT.TOL*DABS(A(JD))) THEN
!     IF(ABS(DG).LT.TOL*ABS(A(JD))) THEN
      WRITE(6,2001) J
      STOP
      END IF
      IF(DG.EQ.0.0D0) THEN
      WRITE(6,2002) J
      STOP
      END IF
!---- REDUCE RHS
500   IF(BACK) B(J) = B(J) - DOT(A(JR+1),B(IS-1),JH-1)
550   CONTINUE
      IF(DG.LT.0.D0) NEIGLO = NEIGLO +1
      IF(DG.NE.0.0D0.AND.AFAC) A(JD) = 1./DG
600   JR = JD
      IF(.NOT.BACK) RETURN
!---- DIVIDE BY DIAGONAL PIVOTS
      DO 700 I = 1,NEQ
      ID = JDIAG(I)
      DG = B(I)
      B(I) = B(I)*A(ID)
700   AENGY = AENGY + B(I)*DG
!---- BACKSUBSTITUTE
      J = NEQ
      JD = JDIAG(J)
800   D = B(J)
      J = J - 1
      IF(J.LE.0) GO TO 1100
      JR = JDIAG(J)
      IF(JD-JR.LE.1) GO TO 1000
      IS = J - JD + JR + 2
      K = JR - IS + 1
      DO 900 I = IS,J
900   B(I) = B(I) - A(I+K)*D
1000  JD = JR
      GO TO 800
1100  IF(PRF) WRITE(6,2003) AENGY
      RETURN
2000  FORMAT('  **WARNING** SIGN OF EQUATION',I5,' CHANGED IN ACTCOL')
2001  FORMAT(///'  **ERROR** DETECTED BY SUBROUTINE ACTCOL'
     1 /' EQUATION',I5,' LOST AT LEAST 7 SIGNIFICANT DIGITS')
2002  FORMAT(///'  **ERROR** DETECTED BY SUBROUTINE ACTCOL'
     1 /' PIVOT OF EQUATION',I5,' IS ZERO')
2003  FORMAT('  **ENERGY COMPUTED IN ACTCOL IS',E18.10)
      END
      SUBROUTINE ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL,AFL,BFL,CFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ASSEMBLE GLOBAL ARRAYS

      LOGICAL AFL,BFL,CFL
      DIMENSION A(*),B(*),JDIAG(*),P(*),S(NST,*),LD(*) ,C(*)
      DO 200 J = 1,NEL
      K = LD(J)
      IF(K.EQ.0) GO TO 200
      IF(BFL)
     1B(K) = B(K) + P(J)
      IF(.NOT.AFL.AND..NOT.CFL) GO TO 200
      L = JDIAG(K) - K
      DO 100 I = 1,NEL
      M = LD(I)
      IF(M.GT.K.OR.M.EQ.0) GO TO 100
      M = L + M
      IF(AFL)
     1A(M) = A(M) + S(I,J)
      IF(CFL)
     1C(M) = C(M) + S(J,I)
100   CONTINUE
200   CONTINUE
      RETURN
      END
      SUBROUTINE PLOAD(ID,IX,PDIST,IDIST,F,B,NEQ,IDPROP,PROP,NEN,
     1 NEN1,NDF,NST,NDIST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- INITIALIZE LOAD VECTOR

      DIMENSION ID(*),F(*),B(*)
      DIMENSION IX(NEN1,*),IDIST(2,*),PDIST(NST,*)
      DIMENSION IDPROP(*),PROP(7)

      DO 10 I = 1,NEQ
!     IF (ID(I).EQ.0) B(I) = P*F(I)
      IF (ID(I).EQ.0) THEN
!---- NO DISPLACEMENT HAS BEEN SPECIFIED.
      IND=IDPROP(I)
      B(I)=PROP(IND)*F(I)
      END IF
   10 CONTINUE

      IF(NDIST.EQ.0) RETURN
!---- ADD TRACTION LOADS
      DO 20 N=1,NDIST
      IELEM=IDIST(1,N)
      IND=IDIST(2,N)
      DO 25 I=1,NEN
      II=IX(I,IELEM)
      IF(II.EQ.0) GO TO 25
!---- II IS GLOBAL NODE# CORR TO LOCAL NODE I OF ELEM IELEM
      DO 26 J=1,NDF
      JS=(II-1)*NDF+J
      KS=(I-1)*NDF+J
!---- JS IS GLOBAL EQN# CORR TO GLOBAL NODE II AND DOF J.
!---- KS IS LOCAL EQN# CORR  TO LOCAL  NODE I  AND DOF J.
      IF(ID(JS).EQ.0) THEN
!---- NO DISPLACEMENT HAS BEEN SPECIFIED AT GLOBAL DOF JS
      B(JS)=B(JS)+PDIST(KS,N)*PROP(IND)
      END IF
26    CONTINUE
25    CONTINUE
20    CONTINUE
      RETURN
      END
      SUBROUTINE PRTDIS(ID,X,B,IDPROP,PROP,F,NDM,NDF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- OUTPUT NODAL VALUES

      LOGICAL PCOMP
!     COMMON/PRLOD/ PROP
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      COMMON /TDATA/ TIME,DT,BET1,GAM1
      DIMENSION X(NDM,*),B(NDF,*),UL(6),ID(NDF,*),F(NDF,*),CD(3),
     1 DI(2),UCLEAN(6)
      DIMENSION IDPROP(NDF,*),PROP(7)
      DATA BL/4HBLAN/,TOL/1.D-7/
      DATA CD/4H COO,4HRDIN,4HATES/,DI/4H DIS,2HPL/
      DO 102 II = 1,NUMNP,50
      WRITE(6,2000) O,HEAD,TIME,(PROP(I),I=1,7),(I,CD(1),CD(2),I=1,
     1 NDM),(I,DI(1),DI(2),I=1,NDF)
      JJ = MIN0(NUMNP,II+49)
      DO 102 N = II,JJ
      IF(PCOMP(X(1,N),BL)) GO TO 101
      DO 100 I = 1,NDF
      IF (ID(I,N).EQ.0) THEN
      UL(I) = B(I,N)
      ELSE
      IND=IDPROP(I,N)
      UL(I) = F(I,N)*PROP(IND)
!     UL(I) = F(I,N)*PROP
      END IF
100   CONTINUE
!---- FILTER OFF NOISE
      UMAX = 0.D0
      DO 200 I = 1,NDF
      UABS = DABS(UL(I))
!     UABS = ABS(UL(I))
      IF (UMAX.LT.UABS) UMAX = UABS
200   CONTINUE
      UMAX = TOL*UMAX
      DO 210 I = 1,NDF
      UCLEAN(I) = UL(I)
      IF (DABS(UCLEAN(I)).LT.UMAX) UCLEAN(I) = 0.D0
!     IF (ABS(UCLEAN(I)).LT.UMAX) UCLEAN(I) = 0.0
210   CONTINUE
      ndm1 = max0(ndm,3)
      WRITE(6,2001) N,(X(I,N),I=1,NDM1),(UCLEAN(I),I=1,NDF)
101   CONTINUE
102   CONTINUE
      RETURN
2000  FORMAT(A1,20A4//5X,19HNODAL DISPLACEMENTS,5X,4HTIME,E13.5/
     1       2X,'PROP. LD.',7E13.5//6X,4HNODE,9(I7,A4,A2))
2001  FORMAT(I10,9E13.4)
      END
      SUBROUTINE PSETM(NA,NE,NJ,AFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PLONG/ KMAX

!---- SET POINTER FOR ARRAYS

      LOGICAL AFL
      NA = MAX0(KMAX,NE)
      NE = NE + NJ
      AFL = .FALSE.
      CALL SETMEM(NE)
      RETURN
      END
      SUBROUTINE UACTCL(A,C,B,JDIAG,NEQ,AFAC,BACK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL AFAC,BACK
      DIMENSION A(*),B(*),JDIAG(*),C(*)

!---- UNSYMMETRIC, ACTIVE COLUMN PROFILE EQUATION SOLVER

!---- FACTOR A TO UT*D*U, REDUCE B TO Y
      JR = 0
      DO 300 J = 1,NEQ
      JD = JDIAG(J)
      JH = JD - JR
      IF(JH.LE.1) GO TO 300
      IS = J + 1 - JH
      IE = J - 1
      IF(.NOT.AFAC) GO TO 250
      K = JR + 1
      ID = 0
!---- REDUCE ALL EQUATIONS EXCEPT DIAGONAL
      DO 200 I = IS,IE
      IR = ID
      ID = JDIAG(I)
      IH = MIN0(ID - IR - 1,I - IS)
      IF(IH.EQ.0) GO TO 150
      A(K) = A(K) - DOT(A(K-IH),C(ID-IH),IH)
      C(K) = C(K) - DOT(C(K-IH),A(ID-IH),IH)
150   IF(A(ID).NE.0.0D0) C(K) = C(K)/A(ID)
200   K = K + 1
!---- REDUCE DIAGONAL TERM
      A(JD) = A(JD) - DOT(A(JR+1),C(JR+1),JH-1)
!---- FORWARD REDUCE THE R.H.S.
250   IF(BACK) B(J) = B(J) - DOT(C(JR+1),B(IS),JH-1)
300   JR = JD
      IF(.NOT.BACK) RETURN
!---- BACKSUBSTITUTION
      J = NEQ
      JD = JDIAG(J)
500   IF(A(JD).NE.0.0D0) B(J) = B(J)/A(JD)
      D = B(J)
      J = J - 1
      IF(J.LE.0) RETURN
      JR = JDIAG(J)
      IF(JD - JR.LE.1) GO TO 700
      IS = J -JD + JR + 2
      K = JR - IS + 1
      DO 600 I = IS,J
600   B(I) = B(I) - A(I+K)*D
700   JD = JR
      GO TO 500
      END
      SUBROUTINE NUMASS(B,NEQ,MQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(*)
      NN = 0
      DO 10 N = 1,NEQ
      IF(B(N).NE.0.0D0) NN = NN + 1
10    CONTINUE
      IF(NN.LT.MQ) WRITE(6,2000) NN
      MQ = MIN0(MQ,NN)
      RETURN
2000  FORMAT(' SUBSPACE SIZE REDUCED TO',I4,' BY NUMBER OF NONZERO LUMPE
     1D MASS TERMS')
      END
      SUBROUTINE PRTREA(R,ID,NDF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- PRINT NODAL REACTIONS

      LOGICAL BCFLG
      DIMENSION R(NDF,*),RSUM(6),ASUM(6),ID(*)
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      DO 50 K = 1,NDF
      RSUM(K) = 0.
50    ASUM(K) = 0.
      WRITE(6,2000) O,HEAD,(K,K=1,NDF)
      KOUNT = 0
      DO 100 N = 1,NUMNP
      BCFLG = .FALSE.
      DO 75 K = 1,NDF
      KOUNT = KOUNT + 1
      R(K,N) = -R(K,N)
      RSUM(K) = RSUM(K) + R(K,N)
      ASUM(K) = ASUM(K) + DABS(R(K,N))
!     ASUM(K) = ASUM(K) + ABS(R(K,N))
      IF (ID(KOUNT).NE.0) BCFLG = .TRUE.
75    CONTINUE
      IF (BCFLG) WRITE(6,2001) N,(R(K,N),K=1,NDF)
100   CONTINUE
!---- PRINT STATICS CHECK
      WRITE(6,2002) (RSUM(K),K=1,NDF)
      WRITE(6,2003) (ASUM(K),K=1,NDF)
      RETURN
2000  FORMAT(A1,20A4//5X,15HNODAL REACTIONS//6X,4HNODE,
     1  6(I9,4H DOF))
2001  FORMAT(I10,6E13.4)
2002  FORMAT(/7X,3HSUM,6E13.4)
2003  FORMAT(/3X,7HABS SUM,6E13.4)
      END
      SUBROUTINE PRTNOD(R,TITLE,NDF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*14 TITLE

!---- PRINT NODAL ARRAYS

      DIMENSION R(NDF,*)
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      DO 100 N = 1,NUMNP,50
      J = MIN0(NUMNP,N+49)
      WRITE(6,2000) O,HEAD,TITLE,(K,K=1,NDF)
      DO 100 I = N,J
100   WRITE(6,2001) I,(R(K,I),K=1,NDF)
      RETURN
2000  FORMAT(A1,20A4//5X,'NODAL ',A12//6X,4HNODE,
     1  6(I9,4H DOF))
2001  FORMAT(I10,6E13.4)
      END
      SUBROUTINE PEIGS(A,B,IDPROP,PROP,F,X,Y,Z,ID,IX,JDIAG,
     1 NDF,NDM,NEN1,DFL,AFR,TOL,ITS,SHFTEG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- OLD VERSION: COMPUTE LOWEST EIGENVALUE BY INVERSE ITERATION
!---- NEW VERSION: COMPUTE ANY EIGENVALUE BY INVERSE VECTOR ITERATION
!---- WITH SHIFTING. IF SHFTEG=0.0, THEN LOWEST EIGENVALUE WILL BE COMPUTED.

      LOGICAL DFL,AFR
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      COMMON /ENGYS/ AENGY
      COMMON /STURM/ NEIGLO
      DIMENSION A(*),B(*),F(*),X(*),Y(*),Z(*),ID(*),IX(*),JDIAG(*)
      DIMENSION IDPROP(*),PROP(7)
!     DATA ITS/100/,TOL/1.D-9/

      WRITE(6,2002) SHFTEG
!---- COMPUTE SHIFTED STIFFNESS MATRIX
      NSTMX = JDIAG(NEQ)
      IF(DFL) THEN
!---- DIAGONAL (LUMPED)MASS MATRIX
      DO 10 I=1,NEQ
      JK = JDIAG(I)
      JM = I
      A(JK) = A(JK) - SHFTEG*B(JM)
10    CONTINUE
      ELSE
!---- CONSISTENT MASS MATRIX
      DO 15 I=1,NSTMX
      A(I) = A(I) - SHFTEG*B(I)
15    CONTINUE
      END IF
!---- GET START VECTOR FROM DIAGONAL OF MASS MATRIX
!     DO 100 I = 1,NEQ
!     J = JDIAG(I)
!     IF(DFL) J = I
!100   Y(I) = B(J)
!---- IN THE PRESENT VERSION, THE START VECTOR 'Y' IS USER-INPUT THRO
!---- MACRO 'INIT' IN PMESH
!---- RESET DOFS CORR TO FIXED SUPPORTS TO ZERO IN START VECTOR 'Y'
      DO 102 I=1,NEQ
      IF(ID(I).NE.0) Y(I) = 0.0
102   CONTINUE

      EIGP = 0.
      CALL ACTCOL(A,Z,JDIAG,NEQ,AFR,.FALSE.)
      WRITE(6,2003) NEIGLO
      DO 200 I = 1,ITS
      CALL PZERO(Z,NEQ)
!---- {Z} = [M]{Y}
      IF(.NOT.DFL)
     1 CALL PROMUL(B,Y,Z,JDIAG,NEQ)
      IF(DFL) THEN
!---- B IS A LUMPED MASS MATRIX
      DO 110 II=1,NEQ
      Z(II) = B(II)*Y(II)
110   CONTINUE
      END IF
!---- RAYLEIGH QUOTIENT
!---- EIG = <Y>[K]{Y}/<Y>[M]<Y>
      EIG = AENGY/DOT(Y,Z,NEQ)
      IF(DABS(EIG-EIGP).LT.TOL*DABS(EIG)) GO TO 300
!     IF(ABS(EIG-EIGP).LT.TOL*ABS(EIG)) GO TO 300
!---- NORMALIZE: {Y} = [M]{Y}/<Y>[M]{Y}
      CALL NORM(Y,Z,NEQ,.TRUE.)
      EIGP = EIG
!---- INVERSE ITERATION
!---- SOLVE FOR {YNEW} FROM: [K]{YNEW} = {Y}
200   CALL ACTCOL(A,Y,JDIAG,NEQ,.FALSE.,.TRUE.)
      WRITE(6,2001) ITS
      RETURN
300   CONTINUE
!---- COMPUTE ACTUAL EIGENVALUE BY ADDING SHIFT
      EIG = EIG + SHFTEG
      WRITE(6,2000) O,HEAD,EIG,I
!---- NORMALIZE EIGEN VECTOR: {Y} = {Y}/<Y>[M]{Y}
      CALL NORM(Z,Y,NEQ,.TRUE.)
      CALL PRTDIS(ID,X,Z,IDPROP,PROP,F,NDM,NDF)
      RETURN
2000  FORMAT(A1,20A4//5X,14HEIGENVALUE =  ,E15.7/5X,14HITERATIONS =  ,
     1   I9/)
2001  FORMAT(5X,57H**FATAL ERROR 09** NO CONVERGENCE IN EIGENVALUES, ITS
     1 =  ,I5)
2002  FORMAT(//,5X,'SHIFT USED IN INVERSE VECTOR ITERATION=',E15.7)
2003  FORMAT(5X,'THERE ARE',I5,2X,'EIGENVALUES LESS THAN THE SHIFT')
      END
      DOUBLE PRECISION FUNCTION PRPLD1(T)
!     REAL FUNCTION PRPLD1(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTE PROPORTIONAL LOAD FACTOR

      COMMON /PTABLE/ NPLD,KPLD,T1,P1,T2,P2
!---- CHECK IF TIME IS PAST ARRIVAL TIME
      IF (T.GT.T1) GO TO 10
!     PROPLD = P1
      PRPLD1=P1
      RETURN
   10 CONTINUE
!---- CHECK WHETHER TIME IS BEYOND CURRENT INTERVAL
      IF (T.LE.T2) GO TO 20
!---- CHECK WHETHER PROP TABLE IS EXHAUSTED
      IF (KPLD.LT.NPLD) GO TO 15
!     PROP = P2
      PRPLD1=P2
      RETURN
   15 CONTINUE
!---- SHIFT ONE TIME INTERVAL FORWARD
      T1 = T2
      P1 = P2
      READ (5,'(2F10.0)') T2,P2
      KPLD = KPLD + 1
      GO TO 10
   20 CONTINUE
!---- CHECK FOR NEGATIVE INTERVALS
      DT = T2 -T1
      IF (DT.GE.0.D0) GO TO 25
      WRITE(5,'(///5X,''**FATAL ERROR 17** NEGATIVE INTERVAL IN '',
     1 ''PROP TABLE'')')
      STOP
   25 CONTINUE
!---- CHECK FOR TIME AT END OF CURRENT INTERVAL
      IF (T.NE.T2) GO TO 30
!     PROPLD = P2
      PRPLD1 = P2
      RETURN
   30 CONTINUE
!---- LINEARLY INTERPOLATE PROP TABLE
      RR = (T - T1)/DT
!     PROPLD = (1. - RR)*P1 + RR*P2
      PRPLD1 = (1. - RR)*P1 + RR*P2
      RETURN
      END

      SUBROUTINE PROPLD(TTABLE,PTABLE,PROP,TIME,NPLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTES PROPORTIONAL LOAD FACTORS
!---- MAXIMUM OF 7 PROP LOAD FACTORS ARE PERMITTED.
!---- NPLD IS MAX # OF PROP/TIME CARDS WHICH ARE INPUT.
!---- (TTABLE(I),I=1,NPLD) CONTAINS TIME VALUES IN
!---- INCREASING ORDER. ((PTABLE(I,J),I=1,7),J=1,NPLD)
!---- CONTAINS CORR VALUES OF THE 7 PROP FACTORS.

      DIMENSION TTABLE(*),PTABLE(7,*),PROP(7)

!---- CHECK IN WHICH INTERVAL THE CURRENT TIME FALLS:
      IF(TIME.GT.TTABLE(1)) GO TO 10
      DO 20 I=1,7
      PROP(I)=PTABLE(I,1)
20    CONTINUE
      RETURN
10    CONTINUE
      IF(TIME.LT.TTABLE(NPLD)) GO TO 30
      DO 40 I=1,7
      PROP(I)=PTABLE(I,NPLD)
40    CONTINUE
      RETURN
30    CONTINUE
      IND=0
50    IND=IND+1
      IF(TTABLE(IND+1).LT.TTABLE(IND)) THEN
      WRITE(6,1000) IND,TTABLE(IND),TTABLE(IND+1)
      STOP
      END IF
      IF(TIME.LE.TTABLE(IND+1)) GO TO 60
      GO TO 50
60    CONTINUE
!---- TTABLE(IND).LT.TIME.LE.TTABLE(IND+1)  WHERE   1.LE.IND.LE.(NPLD-1)
!---- LINEARLY INTERPOLATE FOR OBTAINING PROP LOAD FACTORS
      RR=(TIME-TTABLE(IND))/(TTABLE(IND+1)-TTABLE(IND))
      DO 70 I=1,7
      PROP(I)=(1.D0-RR)*PTABLE(I,IND)+RR*PTABLE(I,IND+1)
70    CONTINUE
      RETURN
1000  FORMAT(///5X,'*** FATAL ERROR 17 *** NEGATIVE INTERVAL IN',
     1 ' PROP TABLE'/5X,'TIME CARD #=',I5,5X,'TIME1=',E13.5,5X,
     2 'TIME2=',E13.5)
      END


!---------------------------------------------------------------------
!---- THE FOLLOWING SEGMENT WAS ADDED BY R.NARASIMHAN ON FEB 11, 1987.
!---- AN OUT-OF-CORE EQUATION SOLVER FOR SYMMETRIC MATRICES, STORED
!---- IN BLOCKED,PROFILE FORM.
!---- ADAPTED FROM: D.P.MONDKAR & G.H.POWELL, COMP & STRUCT, VOL.4,1974,P.699

      SUBROUTINE OPTBLOK(A,B,PIVOT,NCOL,AVDIAG,NEQ,KEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /BLOCK/ MAXAV,LBLOK,NBLOK,MAXC,NELCOR
      COMMON /TAPE/ NMAT,NRED,NPVT,NELSTF,NBFGS,LDUMP,LREST
     1 ,LSEFF,LDISP,LPRTS,LRESF,LSMTH
      DIMENSION A(LBLOK),B(LBLOK),PIVOT(MAXC),NCOL(NBLOK)
      DATA TOL/1.D-09/
!-----------------------------------------------------------------------
!     A LARGE CAPACITY SUBROUTINE TO SOLVE A SYMMETRIC SET OF LINEAR
!     ALGEBRIAC EQUATIONS A.X=B USING BLOCK PARTITIONING.
!----                   ------------
!     THE VARIABLES IN THE ARGUMENT LIST HAVE THE FOLLOWING MEANING:

!        NEQ = NUMBER OF EQUATIONS
!      NBLOK = NUMBER OF MATRIX BLOCKS
!      LBLOK = STORAGE SPACE AVAILABLE FOR A MATRIX BLOCK.
!      LBLOK2= LBLOK/2
!       MAXC = MAXIMUM NUMBER OF COLUMNS IN MATRIX BLOCKS.

!          A = VECTOR OF DIMENSION LBLOK.
!          B = VECTOR OF DIMENSION LBLOK.
!       PIVOT= VECTOR OF DIMENSION MAXC  (SEE NOTE 1).
!       NCOL = VECTOR OF DIMENSION NBLOK, WHICH UPON ENTRY CONTAINS
!              NUMBER OF COLUMNS IN EACH MATRIX BLOCK.

!       NMAT = TAPE NUMBER, WHICH CONTAINS SEQUENTIALLY THE UNREDUCED
!              MATRIX BLOCKS.
!       NRED = TAPE NUMBER, TO STORE SEQUENTIALLY THE REDUCED MATRIX
!              BLOCKS.
!       NPVT = TAPE NUMBER, TO STORE PIVOT BLOCKS.
!        KEX = EXECUTION PARAMETER. IF UPON ENTRY,
!    (1) KEX =1 REDUCTION OF COEFF MATRIX. REQUIRES TAPES
!               NMAT, NRED AND NPVT ONLY.
!    (2) KEX =2 REDUCTION OF LOAD VECTOR (PARTITION TYPE-A)
!               AND BACK SUBSTITUTION. REQUIRES TAPES NRED AND NPVT ONLY.
! -----                     ----------------
!                                NOTES
!    (1) THE ENTRY PIVOT(1) HAS THE SAME ADDRESS AS THE ENTRY B(1).
!    (2) FOR KEX=2, UPON ENTRY VECTOR A CONTAINS THE LOAD VECTOR
!        (PARTITION TYPE- A : A LOAD BLOCK CONTAINS A COMPLETE LOAD
!         VECTOR. (I.E.) PARTITIONING OF LOAD VECTOR IS NOT PERMITTED.)
!        UPON EXIT, VECTOR A CONTAINS THE SOLUTION VECTOR.
! -----                      -------------------
!     TOTAL CORE STORAGE = 2*LBLOK + 1*NBLOK
!     ****************************************************************
!
      LVEC=1
      SOLTOL=AVDIAG*TOL
      DIAMIN=1.D+30
      GO TO (20,280), KEX
!     ****************************************************************
!     REDUCE COEFFICIENT MATRIX
!     ****************************************************************
   20 JF=1
      JL=0
      JCP=0
      REWIND NMAT
      REWIND NRED
      REWIND NPVT

!     MINA(1)=1
!----                 ----------
!----  LOOP OVER UNREDUCED STIFFNESS BLOCKS
      DO 270 NPB=1,NBLOK
      READ(NMAT) A
      NCA=NCOL(NPB)
      JL=JL+NCA
      JNA=LBLOK-NCA
      IF(NPB.EQ.1) GO TO 180
      NAJP=1
      MIN=JF
      JK=JNA
      DO 30 J=JF,JL
      JK=JK+1
      NAJ=IDINT(A(JK))
!     NAJ=IFIX(A(JK))
      IF=J-NAJ+NAJP
      IF (MIN.GT.IF) MIN=IF
30    NAJP=NAJ+1
!     MINA(NPB)=MIN
      NSBL=NPB-1
!----             -------------
!---- LOOP OVER PREVIOUSLY REDUCED BLOCKS
      REWIND NRED
      JCP=0
      DO 100 NSB=1,NSBL
      NCB=NCOL(NSB)
      JC=JCP+NCB
      IF (MIN.GE.JC) GO TO 90
      READ (NRED) B
      JNB=LBLOK-NCB
      JK=JNB-JCP
      JJ=JNA
      NAJP=1
      JCP1=JCP+1
      DO 80 J=JF,JL
      JJ=JJ+1
      NAJ=IDINT(A(JJ))
!     NAJ=IFIX(A(JJ))
      NAJJ=NAJ-J
      IF=NAJP-NAJJ
      IF (IF.GE.JC) GO TO 80
      IS=MAX0(IF+1,JCP1)
      II=IS+JK
      NBIP=1
      IF (II.EQ.JNB+1) GO TO 40
      NBIP=IDINT(B(II-1))+1
!     NBIP=IFIX(B(II-1))+1
40    DO 70 I=IS,JC
      NBI=IDINT(B(II))
!     NBI=IFIX(B(II))
      NBII=NBI-I
      KF=NBIP-NBII
      IF (KF.GE.I) GO TO 60
      KF=NAJJ+MAX0(IF,KF)
      JIA=NAJJ+I
      KL=JIA-1
      JM=NBII-NAJJ
      AA=0.D0
      DO 50 K=KF,KL
50    AA=AA+B(JM+K)*A(K)
      A(JIA)=A(JIA)-AA
60    II=II+1
70    NBIP=NBI+1
80    NAJP=NAJ+1
      GO TO 100
90    READ (NRED)
100   JCP=JC
! ----                 -----------------------
      REWIND NPVT
      JCP=0
      DO 170 NSB=1,NSBL
      NCB=NCOL(NSB)
      JC=JCP+NCB
      IF(MIN.GT.JC) GO TO 160
      READ (NPVT) PIVOT
      JJ=JNA
      JCP1=JCP+1
      NAJP=1
      DO 150 J=JF,JL
      JJ=JJ+1
      NAJ=IDINT(A(JJ))
!     NAJ=IFIX(A(JJ))
      NAJJ=NAJ-J
      IF=NAJP-NAJJ
      IF (IF.GT.JC) GO TO 150
      IL=J-1
      KL=NAJJ+JC
      IF (JF.GT.IL) GO TO 130
      II=JNA
      NAIP=1
      DO 120 I=JF,IL
      II=II+1
      NAI=IDINT(A(II))
!     NAI=IFIX(A(II))
      NAII=NAI-I
      KF=NAIP-NAII
      IF (KF.GT.JC) GO TO 120
      KF=NAJJ+MAX0(IF,KF,JCP1)
      JM=NAII-NAJJ
      AA=0.D0
      DO 110 K=KF,KL
110   AA=AA+A(JM+K)*A(K)
      JIA=NAJJ+I
      A(JIA)=A(JIA)-AA
120   NAIP=NAI+1

130   KF=NAJJ+MAX0(IF,JCP1)
      JM=-(NAJJ+JCP)
      AA=0.D0
      DO 140 K=KF,KL
      DIV=A(K)/PIVOT(JM+K)
      AA=AA+DIV*A(K)
140   A(K)=DIV
      A(NAJ)=A(NAJ)-AA
150   NAJP=NAJ+1
      GO TO 170
160   READ (NPVT)
170   JCP=JC
! -----             ------------------------
180   JF1=JF+1
      IF (JF1.GT.JL) GO TO 250
      JJ=JNA+1
      NAJP=IDINT(A(JJ))
!     NAJP=IFIX(A(JJ))
      JK=JNA-JCP
      DO 240 J=JF1,JL
      JJ=JJ+1
      NAJ=IDINT(A(JJ))
!     NAJ=IFIX(A(JJ))
      NAJJ=NAJ-J
      IF=NAJP-NAJJ
      IF (IF.GE.J) GO TO 240
      IS=MAX0(IF+1,JF1)
      IL=J-1
      IF (IS.GT.IL) GO TO 220
      II=IS+JK
      NAIP=IDINT(A(II-1))
!     NAIP=IFIX(A(II-1))
      DO 210 I=IS,IL
      NAI=IDINT(A(II))
!     NAI=IFIX(A(II))
      NAII=NAI-I
      KF=NAIP-NAII
      IF(KF.GE.I) GO TO 200
      KF=NAJJ+MAX0(IF,KF,JF)
      JIA=NAJJ+I
      KL=JIA-1
      JM=NAII-NAJJ
      AA=0.D0
      DO 190 K=KF,KL
190   AA=AA+A(JM+K)*A(K)
      A(JIA)=A(JIA)-AA
200   II=II+1
210   NAIP=NAI+1

220   KF=MAX0(IF,JF)
      II=KF+JK
      KF=KF+NAJJ
      KL=NAJ-1
      AA=0.D0
      DO 230 K=KF,KL
      NAI=IDINT(A(II))
!     NAI=IFIX(A(II))
      IF(A(NAI).LT.DIAMIN) DIAMIN=A(NAI)
      IF(DIAMIN.LT.SOLTOL) THEN
      WRITE(6,2002)
      KCOLM=II-LBLOK+NCA
      KCOLM=KCOLM+JL-NCA
      WRITE(6,2003) AVDIAG,DIAMIN,SOLTOL,KCOLM
      STOP
      END IF
      DIV=A(K)/A(NAI)
      AA=AA+DIV*A(K)
      A(K)=DIV
230   II=II+1
      A(NAJ)=A(NAJ)-AA
240   NAJP=NAJ+1
!-----             ------------------------
250   WRITE (NRED) A
      DO 260 J=1,NCA
      NAJ=IDINT(A(JNA+J))
!     NAJ=IFIX(A(JNA+J))
      PIVOT(J)=A(NAJ)
      IF(PIVOT(J).LT.DIAMIN) DIAMIN=PIVOT(J)
      IF(DIAMIN.LT.SOLTOL) THEN
      WRITE(6,2002)
      KCOLM=J+JL-NCA
      WRITE(6,2003) AVDIAG,DIAMIN,SOLTOL,KCOLM
      STOP
      END IF
260   CONTINUE
      WRITE (NPVT) PIVOT
!----              --------------------------
270   JF=JL+1
      WRITE(6,2001) AVDIAG,DIAMIN
      RETURN
!     *****************************************************************
!     REDUCE LOAD VECTOR (PARTITION TYPE - A) AND BACK SUBSTITUTE
!     *****************************************************************
280   MIN=NEQ
      JJ=0
      DO 310 J=1,LVEC
      DO 290 I=1,NEQ
      IF(A(JJ+I).NE.0.D0) GO TO 300
290   CONTINUE
      I=NEQ
300   NLV=I
      IF(MIN.GT.I) MIN=I
310   JJ=JJ+NEQ
! -----             --------------------------
      REWIND NRED
      JCP=0
      DO 380 NSB=1,NBLOK
      NCB=NCOL(NSB)
      JC=JCP+NCB
      IF(MIN.GE.JC) GO TO 370
      READ (NRED) B
      JNB=LBLOK-NCB
      JK=JNB-JCP
      JCP1=JCP+1
      JJ=0
      DO 360 J=1,LVEC
      IF=NLV
      IF (IF.GE.JC) GO TO 360
      IS=MAX0(IF+1,JCP1)
      II=IS+JK
      NBIP=1
      IF (II.EQ.JNB+1) GO TO 320
      NBIP=IDINT(B(II-1))+1
!     NBIP=IFIX(B(II-1))+1
320   DO 350 I=IS,JC
      NBI=IDINT(B(II))
!     NBI=IFIX(B(II))
      NBII=NBI-I
      KF=NBIP-NBII
      IF(KF.GE.I) GO TO 340
      KF=JJ+MAX0(IF,KF)
      JIA=JJ+I
      KL=JIA-1
      JM=NBII-JJ
      AA=0.D0
      DO 330 K=KF,KL
330   AA=AA+B(JM+K)*A(K)
      A(JIA)=A(JIA)-AA
340   II=II+1
350   NBIP=NBI+1
360   JJ=JJ+NEQ
      GO TO 380
370   READ (NRED)
380   JCP=JC
! -----         ----------------------
      REWIND NPVT
      JCP=0
      DO 420 NSB=1,NBLOK
      NCB=NCOL(NSB)
      JC=JCP+NCB
      IF(MIN.GT.JC) GO TO 410
      READ (NPVT) PIVOT
      JCP1=JCP+1
      JJ=JCP
      DO 400 J=1,LVEC
      IF=NLV
      IF (IF.GT.JC) GO TO 400
      IS=MAX0(IF,JCP1)-JCP
      DO 390 I=IS,NCB
      JIA=JJ+I
390   A(JIA)=A(JIA)/PIVOT(I)
400   JJ=JJ+NEQ
      GO TO 420
410   READ (NPVT)
420   JCP=JC
! -----           ----------------------------
      N=NEQ
      NPB=NBLOK
      DO 510 NSB=1,NBLOK
      BACKSPACE NRED
      READ (NRED) B
      BACKSPACE NRED
      NCB=NCOL(NPB)
      II=LBLOK
      NBI=IDINT(B(II))
!     NBI=IFIX(B(II))
      IF(NCB.EQ.1) GO TO 470
      DO 460 I=2,NCB
      NBIP=IDINT(B(II-1))
!     NBIP=IFIX(B(II-1))
      NBII=NBI-N
      KF=NBIP-NBII+1
      IF(KF.GE.N) GO TO 450
      JM=NBII
      JJ=0
      DO 440 J=1,LVEC
      JKA=JJ+KF
      JIA=JJ+N
      IKA=JIA-1
      AA=A(JIA)
      DO 430 K=JKA,IKA
430   A(K)=A(K)-B(JM+K)*AA
      JM=JM-NEQ
440   JJ=JJ+NEQ
450   II=II-1
      N=N-1
460   NBI=NBIP
470   IF (NPB.EQ.1) GO TO 510
      NBII=NBI-N
      KF=-NBII+1
      IF (KF.GE.N) GO TO 500
      JM=NBII
      JJ=0
      DO 490 J=1,LVEC
      JKA=JJ+KF
      JIA=JJ+N
      IKA=JIA-1
      AA=A(JIA)
      DO 480 K=JKA,IKA
480   A(K)=A(K)-B(JM+K)*AA
      JM=JM-NEQ
490   JJ=JJ+NEQ
500   N=N-1
      NPB=NPB-1
510   CONTINUE
      RETURN
2001  FORMAT(//5X,'AVDIAG=',E15.7,5X,'DIAMIN=',E15.7)
2002  FORMAT(//5X,'*** MIN DIAG.LT.TOLERANCE IN ROUTINE OPTBLOK',
     1 5X,'EXECUTION TERMINATED'/,3X,'AVDIAG',7X,'DIAMIN',7X,'SOLTOL',
     2 6X,'COLM#')
2003  FORMAT(3E13.5,2X,I5)
      END

      SUBROUTINE ADBLKST(A,B,IX,ID,NCOL,JDIAG,AVDIAG,NST,NEN,NEN1,
     1 NUMEL,NEQ,NDF,ISW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ASSEMBLES ELEM STIFF MATRICES INTO BLOCKED, PROFILE GLOBAL ARRAYS.
!---- ONE BLOCK IS ASSEMBLED AT A TIME. THE ELEM STIFF MATRICES ARE
!---- GROUPED INTO MATRICES B(NST,NST,NELCOR) WHICH HAVE BEEN PREVIOUSLY
!---- WRITTEN (IN ROUTINE PFORM1) INTO TAPE NELSTF. IN THIS ROUTINE,
!---- THE ASSEMBLED GLOBAL ARRAYS, A(LBLOK), WILL BE WRITTEN SEQUENTIALLY
!---- INTO TAPE NMAT.

      COMMON /TAPE/ NMAT,NRED,NPVT,NELSTF,NBFGS,LDUMP,LREST
     1 ,LSEFF,LDISP,LPRTS,LRESF,LSMTH
      COMMON /BLOCK/ MAXAV,LBLOK,NBLOK,MAXC,NELCOR
      DIMENSION A(LBLOK), B(NST,NST,NELCOR), NCOL(NBLOK), IX(NEN1,*),
     1 JDIAG(*),ID(*)

      NCAEND=0
      REWIND NMAT
      AVDIAG=0.D0
!---- LOOP OVER BLOCKS
      DO 100 IBLOK=1,NBLOK
      DO 110 I=1,LBLOK
      A(I)=0.D0
110   CONTINUE
!---- DETERMINE FIRST COLM # NCABEG, LAST COLM # NCAEND OF CURRENT BLOCK.
!---- WRITE DIAGONAL ADDRESSES OF CURRENT BLOCK RIGHT JUSTIFIED IN A(LBLOK).
      NCL=LBLOK-NCOL(IBLOK)
      NCA=NCAEND
      NCABEG=NCAEND+1
      LBLOK1=0
      LDIAG=0
115   NCA=NCA+1
      IF(NCA.GT.NEQ) GO TO 120
      NCOLHT=JDIAG(NCA)
      IF(NCA.GT.1) NCOLHT=NCOLHT-JDIAG(NCA-1)
      LBLOK1=LBLOK1+NCOLHT+1
      IF(LBLOK1.GT.LBLOK) GO TO 120
      LDIAG=LDIAG+NCOLHT
      NCL=NCL+1
      A(NCL)=DFLOAT(LDIAG)
!     A(NCL)=FLOAT(LDIAG)
      GO TO 115
120   CONTINUE
!---- BLOCK IBLOK HAS BEEN FULLY PACKED WITH COLUMNS
      NCA=NCA-1
      NCAEND=NCA
      REWIND NELSTF
!---- LOOP OVER ELEM STIFF GROUPINGS
200   CONTINUE
      READ(NELSTF) NELEM1,NELEM2,B
      KELEM=0
      DO 300 IELEM=NELEM1,NELEM2
      KELEM=KELEM+1
      DO 400 INODE=1,NEN
      II=IX(INODE,IELEM)
      IF(II.EQ.0) GO TO 400
      DO 450 IDF=1,NDF
      IEQ=(II-1)*NDF+IDF
      IF((IEQ.LT.NCABEG).OR.(IEQ.GT.NCAEND)) GO TO 450
      IBLKEQ=IEQ-NCABEG+1
      NCL=LBLOK-NCOL(IBLOK)+IBLKEQ
      IG=IDINT(A(NCL))-IEQ
!     IG=IFIX(A(NCL))-IEQ
      IGPREV=0
      IF(IBLKEQ.GT.1) IGPREV=IDINT(A(NCL-1))
!     IF(IBLKEQ.GT.1) IGPREV=IFIX(A(NCL-1))
      IL=(INODE-1)*NDF+IDF
      DO 500 JNODE=1,NEN
      JJ=IX(JNODE,IELEM)
      IF(JJ.EQ.0) GO TO 500
      DO 550 JDF=1,NDF
      JEQ=(JJ-1)*NDF+JDF
      IF(JEQ.GT.IEQ) GO TO 550
!---- ASSEMBLING CONTRIBUTION TO ENTRY IN ROW JEQ, COLUMN IEQ
!---- (ABOVE DIAGONAL ONLY) IN GLOBAL STIFF MATRIX.
      JL=(JNODE-1)*NDF+JDF
      IG1=IG+JEQ
      IF(IG1.LE.IGPREV) THEN
      WRITE(6,2001) IEQ,JEQ,IG1,IBLOK
      STOP
      END IF
!---- THE LOCATION IG1 IN THE BLOCK IBLOK OF THE PROFILE STIFF ARRAY
!---- CORRESPONDS TO ROW JEQ, COLOUMN IEQ, IN GLOBAL STIFF ARRAY.
      A(IG1)=A(IG1)+B(JL,IL,KELEM)
550   CONTINUE
500   CONTINUE
450   CONTINUE
400   CONTINUE
300   CONTINUE
!---- END OF LOOP OVER ELEM STIFF ARRAYS IN THE CURRENT ELEM GROUPING.
!---- CHECK TO SEE IF THERE ARE ANY MORE ELEMS LEFT TO BE ASSEMLED.
      IF(NELEM2.LT.NUMEL) GO TO 200
!---- ASSEMBLY OF BLOCK IBLOK OF THE PROFILE STIFF ARRAY HAS NOW BEEN
!---- COMPLETED. ADJUST DIAGONAL OF STIFF MATRIX TO ACCOUNT FOR B.C.
!---- ADD UP STIFF VALUES ALONG DIAGONAL TO AVDIAG.
      IF(ISW.EQ.5) GO TO 135
      NCL=LBLOK-NCOL(IBLOK)
      DO 130 NCB=NCABEG,NCAEND
      NCL=NCL+1
      LDIAG=IDINT(A(NCL))
!     LDIAG=IFIX(A(NCL))
      IF(ID(NCB).EQ.0) GO TO 131
      A(LDIAG)=1.D0
131   AVDIAG=AVDIAG+A(LDIAG)
130   CONTINUE
135   CONTINUE
      WRITE(NMAT) A
100   CONTINUE
!---- END OF LOOP OVER STIFF BLOCKS
      AVDIAG=AVDIAG/NEQ
      RETURN
2001  FORMAT(//5X,'*** ERROR IN ASSEMBLY ROUTINE ADBLKST'/5X,
     1 'ATTEMPT TO ASSEMBLE ENTRY IN COL',I6,' ROW',I6,
     2 ' IN THE POSITION',I7,' OF STIFF BLOCK',I6)
      END

      SUBROUTINE BLKCMP(JDIAG,NA,NE,NEQ,NST,IPR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTES NUMBER OF BLOCKS OF PROFILE STIFF ARRAY AND THE
!---- AVAILABLE SPACE FOR EACH BLOCK BASED ON REM CORE SPACE.
!---- TWO BLOCKS SHOULD BE ACCOMODATED IN CORE. THE SIZE OF
!---- A BLOCK MUST BE AT LEAST EQUAL TO THE NUMBER OF EQUATIONS.
!---- SPACE IS PROVIDED AT THE END OF EACH BLOCK TO STORE THE
!---- THE DIAGONAL ADDRESSES OF THE COLS IN THE BLOCK.

      COMMON /PSIZE/ MAX
      COMMON /PLONG/ KMAX
      COMMON /BLOCK/ MAXAV,LBLOK,NBLOK,MAXC,NELCOR
      DIMENSION JDIAG(*)
      NAD=JDIAG(NEQ)
      NA=NE
      NE=NA+(NAD+NEQ)*IPR*2+1
      NREMSP=MAX-NE
      IF(NREMSP.GE.0) THEN
      NBLOK=1
      LBLOK=NAD+NEQ
      IF(LBLOK.LT.NST*NST) THEN
      LBLOK=NST*NST
      NE=NA+LBLOK*IPR*2+1
      IF(NE.GT.MAX) THEN
      WRITE(6,2007) NE,MAXAV
      STOP
      END IF
      END IF
      WRITE(6,2006) NBLOK,LBLOK
      CALL SETMEM(NE)
      RETURN
      END IF
!---- DETERMINE NUMBER OF BLOCKS
!---- FIRST ESTIMATE OF NUMBER OF BLOCKS IS NBLOK1
      NBLOK1=NAD/NEQ+1
      LBLOK=(MAXAV-NA-NBLOK1)/(2*IPR)
      IF(LBLOK.LT.NEQ) THEN
      NE=NA+2*IPR*NEQ+NBLOK1
      WRITE(6,2005) LBLOK,NEQ,NE
      STOP
      END IF
      IF(LBLOK.LT.NST*NST) THEN
      LBLOK=NST*NST
      NE=NA+2*IPR*LBLOK+NBLOK1
      IF(NE.GT.MAX) THEN
      WRITE (6,2007) NE,MAXAV
      STOP
      END IF
      END IF
      NBLOK=0
      NCA=0
610   CONTINUE
      NBLOK=NBLOK+1
      LBLOK1=0
615   NCA=NCA+1
      IF(NCA.GT.NEQ) GO TO 625
      NCOLHT=JDIAG(NCA)
      IF(NCA.GT.1) NCOLHT=NCOLHT-JDIAG(NCA-1)
      LBLOK1=LBLOK1+NCOLHT+1
      IF(LBLOK1.GT.LBLOK) GO TO 620
      GO TO 615
620   CONTINUE
!---- BLOCK NBLOK HAS BEEN COMPLETED. BUT ALL COLS IN STIFF ARRAY HAVE
!---- NOT BEEN ACCOMODATED. SO, GO TO 610 AND FILL NEX BLOCK.
      NCA=NCA-1
      GO TO 610
625   CONTINUE
!---- ALL COLS IN STIFF ARRAY HAVE NOW BEEN ACCOUNTED FOR.
      NE=NA+2*LBLOK*IPR+NBLOK
      WRITE(6,2006) NBLOK,LBLOK
      IF(NE.GT.MAX) THEN
      WRITE(6,2007) NE,MAXAV
      STOP
      END IF
      CALL SETMEM(NE)
      RETURN
2005  FORMAT(//5X,'**SIZE OF STIFF BLOCK',I7,' IS LESS THAN NEQ=',I7/
     1 5X,'INCREASE CORE SIZE MDIM TO',I8)
2006  FORMAT(//5X,'BASED ON AVAILABLE CORE SPACE, THE STIFF MATRIX',
     1 ' WILL BE BLOCKED INTO',I7/,5X,'BLOCKS OF SIZE',I7,' EACH.')
2007  FORMAT(//5X,'*** INCREASE CORE SIZE MDIM TO',I8,
     1 ' KEEP MDIMAV AT',I7)
      END

      SUBROUTINE BLKCOL(JDIAG,NCOL,NEQ,NST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- DETERMINES NUMBER OF COLS IN EACH STIFF BLOCK, NCOL(NBLOK).
!---- MAXC IS THE MAXIMUM NUMBER OF COLS IN ANY STIFF BLOCK.
!---- ALSO DETERMINES THE NUMBER OF ELEM STIFF ARRAYS WHICH WILL
!---- FIT IN CORE (NELCOR).

      COMMON /BLOCK/ MAXAV,LBLOK,NBLOK,MAXC,NELCOR
      DIMENSION JDIAG(NEQ),NCOL(NBLOK)
      MAXC=1
      NCA=0
!---- DETERMINE NUMBER OF COLS FOR EACH BLOCK
      DO 100 IBLOK=1,NBLOK
      LBLOK1=0
      NCOL(IBLOK)=0
115   NCA=NCA+1
      IF(NCA.GT.NEQ) GO TO 120
      NCOLHT=JDIAG(NCA)
      IF(NCA.GT.1) NCOLHT=NCOLHT-JDIAG(NCA-1)
      LBLOK1=LBLOK1+NCOLHT+1
      IF(LBLOK1.GT.LBLOK) GO TO 120
      NCOL(IBLOK)=NCOL(IBLOK)+1
      GO TO 115
120   CONTINUE
      MAXC=MAX0(MAXC,NCOL(IBLOK))
      NCA=NCA-1
100   CONTINUE
      NELCOR=LBLOK/(NST*NST)
      NELCOR=MAX0(NELCOR,1)
      WRITE(6,2001) MAXC,NELCOR
      RETURN
2001  FORMAT(//5X,'MAX NUMBER OF COLS IN ANY STIFF BLOCK (MAXC)=',I6/
     1     1X,'MAX NUMBER OF ELEM STIFF ARRAYS IN CORE (NELCOR)=',I6)
      END

      SUBROUTINE PRBLKML(A,B,C,NCOL,NEQ,NTAPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ROUTINE TO FORM C = C + A*B, WHERE A IS A SYMM SQUARE MATRIX
!---- STORED IN BLOCKS, IN PROFILE FORM, IN TAPE NTAPE AND B AND C
!---- ARE VECTORS. THE DIAGONAL ADDRESSES OF EACH STIFF BLOCK IS
!---- STORED AT THE END OF THE BLOCK.

      COMMON /BLOCK/ MAXAV,LBLOK,NBLOK,MAXC,NELCOR
      DIMENSION A(LBLOK),B(*),C(*),NCOL(NBLOK)
      REWIND NTAPE
      NCABEG=1
      DO 300 IBLOK=1,NBLOK
      READ (NTAPE) A
      NCAEND=NCABEG+NCOL(IBLOK)-1
      NCL=LBLOK-NCOL(IBLOK)
      JS=1
      DO 200 J=NCABEG,NCAEND
      NCL=NCL+1
      JD=IDINT(A(NCL))
!     JD=IFIX(A(NCL))
      IF(JS.GT.JD) GO TO 200
      BJ=B(J)
      AB=A(JD)*BJ
      IF(JS.EQ.JD) GO TO 150
      JB=J-JD
      JE=JD-1
      DO 100 JJ=JS,JE
      AB=AB+A(JJ)*B(JJ+JB)
100   C(JJ+JB)=C(JJ+JB)+A(JJ)*BJ
150   C(J)=C(J)+AB
200   JS=JD+1
      NCABEG=NCAEND+1
300   CONTINUE
      RETURN
      END

      SUBROUTINE PSWITCH(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- SWITCHES ARRAYS A AND B.

      DIMENSION A(*),B(*)
      DO 100 I=1,N
      AA=A(I)
      A(I)=B(I)
      B(I)=AA
100   CONTINUE
      RETURN
      END

      SUBROUTINE PFORM1(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,
     1   JDIAG,STR,EPS,Q,B,A,A2,NCOL,VELG,ACCELG,VEL,ACCEL,AVDIAG,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,ISW,U,UD,AFL,BFL,CFL,DFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTE ELEMENT ARRAYS AND ASSEMBLE GLOBAL ARRAYS
!---- NELCOR ELEM STIFF MATRICES ARE GROUPED IN CORE. THE GLOBAL
!---- STIFF MATRICES ARE ASSEMBLED IN BLOCKED, PROFILE FORM USING
!---- ROUTINE ADBLKST. THE FORCE AND DISPLACEMENT ARRAYS ARE STORED
!---- ENTIRELY IN-CORE. GLOBAL FORCE ARRAY IS ASSEMBLED IN ROUTINE
!---- ADDFRC. SPECIFIED DISPL ARE ACCOUNTED IN THE FORCE ARRAY IN
!---- ROUTINE DIAG1.

      LOGICAL AFL,BFL,CFL,DFL,FLAG
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
!     COMMON/PRLOD/ PROP
      COMMON /CNVG1/ICONV,IDIVG,IPATH
      COMMON /TAPE/ NMAT,NRED,NPVT,NELSTF,NBFGS,LDUMP,LREST
     1 ,LSEFF,LDISP,LPRTS,LRESF,LSMTH
      COMMON /BLOCK/ MAXAV,LBLOK,NBLOK,MAXC,NELCOR
      DIMENSION XL(NDM,*),LD(NDF,*),P(*),S(NST,*),IE(7,*),D(40,*),
     1 ID(NDF,*),X(NDM,*),IX(NEN1,*),F(NDF,*),JDIAG(*),B(*),
     2 UL(NDF,*),TL(*),T(*),U(NDF,*),UD(NDF,*)
      DIMENSION A(LBLOK),A2(NST,NST,NELCOR),NCOL(NBLOK)
      DIMENSION PROP(7),IDPROP(NDF,*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(NDF,*),ACCELG(NDF,*),VEL(NDF,*),ACCEL(NDF,*)
!---- LOOP ON ELEMENTS
      IEL = 0
      REWIND NELSTF
      KELEM=0
      NELEM1=1
      NELEM2=1
      DO 110 N = 1,NUMEL
!---- SET UP LOCAL ARRAYS
      KELEM=KELEM+1
      NELEM2=N
      MA = IX(NEN1,N)
      FLAG = .FALSE.
      DO 108 I = 1,NEN
      II = IX(I,N)
      IF(II.NE.0) GO TO 105
      TL(I) = 0.
      DO 103 J = 1,NDM
103   XL(J,I) = 0.
      DO 104 J = 1,NDF
      VEL(J,I) = 0.D0
      ACCEL(J,I) = 0.D0
      UL(J,I) = 0.D0
      UL(J,I+NEN)=0.D0
104   LD(J,I) = 0
      GO TO 108
105   CONTINUE
      NEL = I
      TL(I) = T(II)
      DO 106 J = 1,NDM
106   XL(J,I) = X(J,II)
      DO 107 J = 1,NDF
      JJ = IE(J,MA)
      IF (ID(JJ,II).NE.0) FLAG = .TRUE.
      VEL(J,I) = VELG(JJ,II)
      ACCEL(J,I) = ACCELG(JJ,II)
      IF(JJ.LE.0) GO TO 107
      UL(J,I+NEN)=UD(JJ,II)
      UL(J,I) = U(JJ,II)
107   LD(J,I) = (II-1)*NDF + JJ
108   CONTINUE
C     IF (ISW.EQ.5) FLAG = .FALSE.
!---- FORM ELEMENT ARRAY
      IF(IE(7,MA).NE.IEL) MCT = 0
      IEL = IE(7,MA)
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,ISW)
      IF ((.NOT.FLAG).AND.(BFL))
     1 CALL ADDFRC(B,P,LD,NST,NEL*NDF)
!---- IF AT LEAST 1 DOF IS RESTRAINED AND ELEM STIFF MX HAS BEEN COMPUTED:
      IF (FLAG.AND.(AFL.OR.CFL)) THEN
!---- MODIFY ELEM STIFF MX FOR SPECIFIED DISPL. IF (BFL), THE ELEM
!---- INTERNAL FORCE VECTOR IS ALSO MODIFIED FOR SPECIFIED DISPL.
      CALL MODIFY(S,P,IDPROP,PROP,F,U,ID,IX,N,NST,NEN,NDF,NEN1,AFL,
     1 BFL,CFL)
!---- IF (BFL), ASSEMBLE ELEM INTERNAL FORCE VECTOR INTO GLOBAL ARRAY:
      IF(BFL)
     1 CALL ADDFRC(B,P,LD,NST,NEL*NDF)
      END IF
!---- IF AT LEAST 1 DOF IS RESTRAINED AND ELEM STIFF MX HAS NOT BEEN COMPUTED:
      IF (FLAG.AND.BFL.AND.(.NOT.AFL).AND.(.NOT.CFL)) THEN
!---- FIRST ASSEMBLE ELEMENT INTERNAL FORCE VECTOR INTO GLOBAL ARRAY:
      CALL ADDFRC(B,P,LD,NST,NEL*NDF)
      IF (ISW.NE.5) THEN
!---- COMPUTE ELEM STIFF MX:
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,3)
      DO 120 I = 1,NST
120   P(I) = 0.D0
!---- COMPUTE ELEM FORCE VECTOR CORR TO SPECIFIED DISPL
      CALL MODIFY(S,P,IDPROP,PROP,F,U,ID,IX,N,NST,NEN,NDF,NEN1,AFL,
     1 BFL,CFL)
!---- ASSEMBLE THE ABOVE ELEM FORCE VECTOR INTO GLOBAL ARRAY:
      CALL ADDFRC(B,P,LD,NST,NEL*NDF)
      END IF
      END IF

!---- STORE COMPUTED ELEM STIFF MX IN BLOCK A2
      IF(AFL.OR.CFL) THEN
      DO 121 J1=1,NST
      DO 122 J2=1,NST
      A2(J1,J2,KELEM)=S(J1,J2)
122   CONTINUE
121   CONTINUE
      IF((KELEM.EQ.NELCOR).OR.(N.EQ.NUMEL)) THEN
      WRITE(NELSTF) NELEM1,NELEM2,A2
      KELEM=0
      NELEM1=N+1
      DO 124 I=1,NELCOR
      DO 125 J1=1,NST
      DO 125 J2=1,NST
      A2(J1,J2,I)=0.D0
125   CONTINUE
124   CONTINUE
      END IF
      END IF
110   CONTINUE
      IF(AFL.OR.CFL) THEN
      CALL ADBLKST(A,A2,IX,ID,NCOL,JDIAG,AVDIAG,NST,NEN,NEN1,NUMEL,
     1 NEQ,NDF,ISW)
      END IF
      IF (ISW.EQ.5) RETURN
      IF(BFL)
     1 CALL DIAG1(B,IDPROP,PROP,F,U,ID,NDF,NUMNP)
      RETURN
      END

      SUBROUTINE DIAG1(B,IDPROP,PROP,F,U,ID,NDF,NUMNP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ADJUST RHS TO ACCOUNT FOR BC.

      DIMENSION B(*),ID(NDF,*),F(NDF,*),U(NDF,*),IDPROP(NDF,*),PROP(7)

      IDF=0
      DO 10 J=1,NUMNP
      DO 10 I=1,NDF
      IDF=IDF+1
      IF(ID(I,J).EQ.0) GO TO 10
      IND=IDPROP(I,J)
      B(IDF)=PROP(IND)*F(I,J)-U(I,J)
10    CONTINUE
      RETURN
      END

      SUBROUTINE ADDFRC(B,P,LD,NST,NEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- ASSEMBLE GLOBAL FORCE VECTOR

      DIMENSION B(*),P(*),LD(*)
      DO 200 J=1,NEL
      K=LD(J)
      IF(K.EQ.0) GO TO 200
      B(K)=B(K)+P(J)
200   CONTINUE
      RETURN
      END

      SUBROUTINE PFORMS(NPRT,UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,
     1   PROP,F,T,JDIAG,STR,EPS,Q,B,A,C,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,ISW,U,UD,AFL,BFL,CFL,DFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTE ELEMENT ARRAYS OF NSET AND ASSEMBLE GLOBAL ARRAYS

      LOGICAL AFL,BFL,CFL,DFL,FLAG
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL

!     COMMON/PRLOD/ PROP
      COMMON /CNVG1/ICONV,IDIVG,IPATH
      COMMON /PRTSET/ NSET(9,201)
      COMMON /XJLINT/ NNN,S11(50),S12(50),S22(50),X1(50),X2(50),
     1  DU12(50),DU22(50),WSE(50)
      DIMENSION XL(NDM,*),LD(NDF,*),P(*),S(NST,*),IE(7,*),D(40,*),
     1 ID(NDF,*),X(NDM,*),IX(NEN1,*),F(NDF,*),JDIAG(*),B(*),A(*),C(*),
     2 UL(NDF,*),TL(*),T(*),U(NDF,*),UD(NDF,*)
      DIMENSION PROP(7),IDPROP(NDF,*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(NDF,*),ACCELG(NDF,*),VEL(NDF,*),ACCEL(NDF,*)
!---- LOOP ON ELEMENTS
      IEL = 0
      DO 110 NNN = 1,NSET(NPRT,1)
      N=NSET(NPRT,NNN+1)
!---- SET UP LOCAL ARRAYS
      MA = IX(NEN1,N)
      FLAG = .FALSE.
      DO 108 I = 1,NEN
      II = IX(I,N)
      IF(II.NE.0) GO TO 105
      TL(I) = 0.
      DO 103 J = 1,NDM
103   XL(J,I) = 0.
      DO 104 J = 1,NDF
      VEL(J,I) = 0.D0
      ACCEL(J,I) = 0.D0
      UL(J,I) = 0.
      UL(J,I+NEN)=0.
104   LD(J,I) = 0
      GO TO 108
105   CONTINUE
      NEL = I
      TL(I) = T(II)
      DO 106 J = 1,NDM
106   XL(J,I) = X(J,II)
      DO 107 J = 1,NDF
      JJ = IE(J,MA)
      IF (ID(JJ,II).NE.0) FLAG = .TRUE.
      VEL(J,I) = VELG(JJ,II)
      ACCEL(J,I) = ACCELG(JJ,II)
      IF(JJ.LE.0) GO TO 107
      UL(J,I+NEN)=UD(JJ,II)
      UL(J,I) = U(JJ,II)
107   LD(J,I) = (II-1)*NDF + JJ
108   CONTINUE
!     IF (ISW.EQ.5) FLAG = .FALSE.
!---- FORM ELEMENT ARRAY
      IF(IE(7,MA).NE.IEL) MCT = 0
      IEL = IE(7,MA)
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,ISW)

!---- CHECK FOR PLASTIC DISSIPATION AND TOO LARGE PLASTIC STRAIN INCREMENTS
      IF(IDIVG.EQ.999) RETURN

      IF ((.NOT.FLAG).AND.(AFL.OR.BFL.OR.CFL))
     1 CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      IF (FLAG.AND.(AFL.OR.CFL)) THEN
      CALL MODIFY(S,P,IDPROP,PROP,F,U,ID,IX,N,NST,NEN,NDF,NEN1,AFL,
     1 BFL,CFL)
      CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      END IF
      IF (FLAG.AND.BFL.AND.(.NOT.AFL).AND.(.NOT.CFL)) THEN
      CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      IF (ISW.NE.5) THEN
      CALL ELMLIB(D(1,MA),UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,3)
      DO 120 I = 1,NST
120   P(I) = 0.D0
      CALL MODIFY(S,P,IDPROP,PROP,F,U,ID,IX,N,NST,NEN,NDF,NEN1,AFL,
     1 BFL,CFL)
      CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
      END IF
      END IF
110   CONTINUE
      IF (ISW.EQ.5) RETURN
      IF(AFL.OR.BFL.OR.CFL)
     1 CALL DIAG(A,B,C,IDPROP,PROP,F,U,JDIAG,ID,NDF,NUMNP,AFL,BFL,CFL)
      RETURN
      END

      DOUBLE PRECISION FUNCTION SEF(SIG)
!     REAL FUNCTION SEF(SIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!---- COMPUTE EFFECTIVE STRESS

      DIMENSION SIG(6),SIGD(6)

      P=(SIG(1)+SIG(2)+SIG(3))/3.D0
      DO 10 I=1,3
      SIGD(I)=SIG(I)-P
10    SIGD(I+3)=SIG(I+3)
      XJ2=0.D0
      DO 20 I=1,3
      XJ2=XJ2+SIGD(I)*SIGD(I)
      XJ2=XJ2+2.D0*SIGD(I+3)*SIGD(I+3)
20    CONTINUE
      XJ2=DSQRT(1.5*XJ2)
!     XJ2=SQRT(1.5*XJ2)
      SEF=XJ2
      RETURN
      END
