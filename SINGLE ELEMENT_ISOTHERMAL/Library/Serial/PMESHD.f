      SUBROUTINE PMESH(IDL,XL,PDIST,IDIST,S,IE,D,ID,X,IX,IDPROP,F,T,
     1 DR,STR,EPS,Q,VEL,ACCEL,NDF,NDM,NEN1,NSTR,NQ,NST,NDIST,III,PRTN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- DATA INPUT ROUTINE FOR MESH DESCRIPTION
C
      LOGICAL PRTN,PRT,ERR,PCOMP
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
      COMMON /PMPRTN/ PRT
      COMMON /PRTSET/ NSET(9,201)
      COMMON /WEDGE/ NSURF,ISURF(50),CPSURF(2,50),NDSURF(50),
     1    THTSRF,PENSRF
      COMMON /CNTAC/ CP(2,50),BET(50),DEL(2,50),IC(50),IM(100),
     1 IAB(2,50),MSEG(50),NDEST(50),ANG(50),DUSL(2,50),
     1  SMU,DMU,NC,NM,ALP,FRTOL
      COMMON /COHES/ NINTND,INTND(2,100),ALPAA,DLTAA,SIGMXX,
     1 XLNTND(100),ANGINT(100),CPINT(2,100),
     2 XLN(100), ISTAT(100)
      DIMENSION IE(7,*),D(40,*),ID(NDF,*),X(NDM,*),IX(NEN1,*),
     1 IDL(*),F(NDF,*),T(*),VA(2),XM(3),XV(3),WD(26),
     2 CD(3),TE(3),FD(3),U0(3)
      DIMENSION PDIST(NST,*),IDIST(2,*),S(NST,*),XL(NDM,*)
      DIMENSION IDPROP(NDF,*), PD(3)
      DIMENSION DR(NDF,*),STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VEL(NDF,*),ACCEL(NDF,*)
      DATA WD/4HCOOR,4HELEM,4HMATE,4HBOUN,4HFORC,4HTEMP,4HEND ,4HPRIN,
     1        4HNOPR,4HPAGE,4HBLOC,4HPOLA,4HINIT,4HSTRE,4HSTRA,4HIVAR,
     1        4HMP01,4HVELO,4HMASS,4HMP02,4HDIST,4HDOMA,4HSET ,4HCONT,
     1        4HINTF,4HSURF/
      DATA BL/4HBLAN/,LIST/26/
      DATA DOF/3HDOF/
      DATA VA/4H VAL,2HUE/
      DATA XM/4H MAS,4HSES ,4H    /,XV/4H VEL,4HOCIT,4HIES /,
     1     CD/4H COO,4HRDIN,4HATES/,TE/4H TEM,4HPERA,4HTURE/,
     2     FD/4H FOR,4HCE/D,4HISPL/,U0/4H DIS,4HPL. ,4H    /,
     3     PD/4H PRP,4HID  ,4H    /
      DATA PI/3.141592653/
C---- INITIALIZE ARRAYS
      ERR = .FALSE.
      IF(III.LT.0) GO TO 10
      PRT = .TRUE.
      DO 101 N = 1,NUMNP
      DO 100 I = 1,NDF
      IDPROP(I,N)=1
      ID(I,N) = 0
      DR(I,N) = 0.D0
      VEL(I,N) = 0.D0
      ACCEL(I,N) = 0.D0
100   F(I,N) = 0.D0
      IF(III.EQ.0) X(1,N) = BL
101   IF(III.EQ.0) T(N) = 0.D0
C---- READ MACRO CARDS
10    CONTINUE
      READ(5,1000) CC
      DO 20 I = 1,LIST
20    IF(PCOMP(CC,WD(I))) GO TO 30
      GO TO 10
C---- PROCESS MACROS
30    GO TO (1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,22,23,24,
     1 25,26,27,28),I
C
C---- MACRO 'COOR'
C---- NODAL COORDINATE DATA INPUT
1     CALL GENVEC(NDM,X,CD,PRT,ERR,.TRUE.)
      GO TO 10
C
C---- MACRO 'ELEM'
C---- ELEMENT DATA INPUT
2     CALL GRID(IDL,IX,NEN1,PRT,ERR)
      GO TO 10
C
C---- MACRO 'MATE'
C---- MATERIAL DATA INPUT
3     IF (PRT) WRITE(6,2004) O,HEAD
      DO 304 N = 1,NUMMAT
      READ(5,1002) MA,IEL,(IDL(I),I=1,NDF)
      IF(MA.LE.0) GO TO 10
      DO 300 I = 1,NDF
300   IE(I,MA) = IDL(I)
      DO 301 I = 1,NDF
      IF(IDL(I).NE.0) GO TO 303
301   CONTINUE
C---- RESET ALL ZERO INPUTS
      DO 302 I = 1,NDF
302   IE(I,MA) = I
303   IE(7,MA) = IEL
      IF (PRT) WRITE(6,2003)MA,IEL,(I,IE(I,MA),I=1,NDF)
304   CALL ELMLIB(D(1,MA),DUM,X,DUM,DUM,DUM,IX,DUM,DUM,DUM,DUM,DUM,
     1 NDF,NDM,1,NSDM,NQDM,1)
      GO TO 10
C
C---- MACRO 'BOUN'
C---- READ IN THE RESTRAINT CONDITIONS FOR EACH NODE
    4 III = 1
      CALL BOUND(IDL,ID,NDF,PRT,ERR)
      GO TO 10
C
C---- MACRO 'FORC'
C---- FORCE/DISPL DATA INPUT
5     CONTINUE
C     CALL GENVEC(NDF,F,FD,PRT,ERR,.FALSE.)
      CALL GENFRC(NDF,F,IDPROP,FD,PD,PRT,ERR,.FALSE.)
      GO TO 10
C
C---- MACRO 'TEMP'
C---- TEMPERATURE DATA INPUT
6     CALL GENVEC(1,T,TE,PRT,ERR,.FALSE.)
      GO TO 10
C
C---- MACRO 'END '
C---- TERMINATE MESH INPUT
7     IF(ERR) STOP
      PRTN = PRT
      RETURN
C
C---- MACRO 'PRIN'
C---- SET UP PRINT OPTION
8     PRT = .TRUE.
      GO TO 10
C
C---- MACRO 'NOPR'
C---- TURN OFF PRINT OPTION
9     PRT = .FALSE.
      GO TO 10
C---- MACRO 'PAGE'
C---- SET CARRIAGE CONTROL CHARACTER
11    CONTINUE
      READ(5,1000) O
      GO TO 10
C
C---- MACRO 'BLOC'
C---- GENERATE BLOCK OF NODES AND 4-NODE ELEMENTS
12    IF(III.LT.0) WRITE(6,3003)
      CALL BLKGEN(NDM,NDF,NEN,NEN1,X,IX,PRT)
      GO TO 10
C
C---- MACRO 'POLA'
C---- POLAR COORDINATE NODAL INPUT
13    CALL POLAR(X,NDM,PRT)
      GO TO 10
C
C---- MACRO 'INIT'
C---- READ INITIAL CONDITIONS
   14 CALL GENVEC(NDF,DR,U0,PRT,ERR,.FALSE.)
      GO TO 10
C
C---- MACRO 'STRE'
C---- READ INITIAL STRESSES
   15 CALL GENEL(NSDM,NQUAD,NSTR,'INITIAL STRESSES  ',STR,PRT,ERR)
      GO TO 10
C
C---- MACRO 'STRA'
C---- READ INTIAL STRAINS
   16 CALL GENEL(NSDM,NQUAD,NSTR,'INITIAL STRAINS   ',EPS,PRT,ERR)
      GO TO 10
C
C---- MACRO 'IVAR'
C---- READ INITIAL VALUES OF INTERNAL VARIABLES
   17 CALL GENEL(NQDM,NQUAD,NQ  ,'INTERNAL VARIABLES',Q,PRT,ERR)
      GO TO 10
C
C---- MACRO 'MP01'
C---- USER SUPPLIED MACRO
   18 CONTINUE
      CALL MP01(IDL,IE,D,ID,X,IX,F,T,DR,STR,EPS,Q,VEL,ACCEL,
     1 NDF,NDM,NEN1,NSTR,NQ,III,PRTN)
      GO TO 10
C
C---- MACRO 'VELO'
C---- READ INITIAL NODAL VELOCITIES
   19 CONTINUE
      CALL GENVEC(NDF,VEL,XV,PRT,ERR,.FALSE.)
      GO TO 10
C
C---- MACRO 'MASS'
C---- READ NODAL LUMPED MASSES (ACTIVATED THROUGH MACRO 'LMAS')
   21 CONTINUE
      CALL GENVEC(NDF,ACCEL,XM,PRT,ERR,.FALSE.)
      GO TO 10
C
C---- MACRO 'MP02'
C---- USER SUPPLIED MACRO
   22 CONTINUE
      CALL MP02(IDL,IE,D,ID,X,IX,F,T,DR,STR,EPS,Q,VEL,ACCEL,
     1 NDF,NDM,NEN1,NSTR,NQ,III,PRTN)
      GO TO 10
C
C---- MACRO 'DIST'
C---- INPUT DISTRIBUTED EDGE LOADS AND COMPUTE EQUIV ELEM NODAL LOADS.
C---- CAUTION: THIS MACRO SHOULD BE USED ONLY AFTER 'ELEM' AND
C---- 'MATE'.
C---- GENERAL INSTRUCTIONS FOR PREPARING DATA CARDS FOR 'DIST':
C---- 1) ASSOCIATE EACH OF THE NDIST DISTRIBUTED LOAD CARDS WITH
C----    AN ELEMENT# AND A PROPORTIONAL LOAD TABLE #.
C---- 1A)FORMAT(16I5) (IDIST(1,I),IDIST(2,I),I=1,NDIST)
C---- 2) IF AN ELEMENT HAS DISTRIBUTED LOADS ACTING ON SEVERAL
C----    EDGES, THEN:
C---- 2A)INPUT ONLY ONE SET OF CARDS IF ALL THE DISTRIBUTED LOADS
C----    REFER TO THE SAME PROP LOAD TABLE.
C---- 2B)OTHERWISE, INPUT AS MANY SETS OF CARDS AS THERE ARE
C----    PROP LOAD TABLE REFERENCES.
C---- 3) IF ON THE SAME EDGE OF AN ELEMENT, THE TANGENTIAL AND NORMAL
C----    TRACTIONS VARY DIFFERENTLY WITH TIME (I.E. REFER TO DIFFERENT
C----    PROP LOAD TABLES), 2 SETS OF CARDS SHOULD BE GIVEN.
C---- 4) FOR SMALL STRAIN, PLANE STRESS/PLANE STRAIN/AXISYMM ELEM:
C---- 4A)CARD1: # OF EDGES WITH DISTRIBUTED LOADS
C----    FORMAT(I5): NEDGE
C---- 4B)CARD2: 3 GLOBAL NODE #S ALONG EDGE (ORDER ANTI-CLOCKWISE WRT ELEM),
C----           TANGENTIAL TRACTIONS AT THE 3 NODES, NORMAL TRACTION AT
C----           THE THREE NODES. GIVE ZERO CORR TO MISSING MID-SIDE NODES.
C----           GIVE NEDGE NUMBER OF CARD2
C----    FORMAT(3I5,6F10.0) (NODEG(I),I=1,3),(PRESST(I),I=1,3),(PRESSN(I),
C----                        I=1,3)
C
   23 CONTINUE
      IF(NDIST.LE.0) THEN
      WRITE(6,2010) NDIST
      STOP
      END IF
C---- READ ELEM #S, PROP ID#S CORR TO NDIST DISTRIBUTED LOADS.
      READ(5,1002) (IDIST(1,I),IDIST(2,I),I=1,NDIST)
      DO 230 I=1,NDIST
      IF((IDIST(1,I).LE.0).OR.(IDIST(1,I).GT.NUMEL)) THEN
      WRITE(6,2011) IDIST(1,I)
      STOP
      END IF
      IF((IDIST(2,I).LE.0).OR.(IDIST(2,I).GT.7)) THEN
      WRITE(6,2012) IDIST(2,I)
      STOP
      END IF
230   CONTINUE
C---- LOOP OVER DISTRIBUTED LOAD CARDS
      IF(PRT) WRITE(6,2014)
      DO 231 KDIST=1,NDIST
C---- IDENTIFY ELEM# AND MATERIAL GROUP# TO WHICH THIS ELEM BELONGS.
      N=IDIST(1,KDIST)
      MA=IX(NEN1,N)
C---- IDENTIFY ELEM TYPE
      IEL=IE(7,MA)
C---- SET UP LOCAL COORDINATE ARRAY FOR ELEMENT
      DO 238 I=1,NEN
      II=IX(I,N)
      IF(II.NE.0) GO TO 235
      DO 233 J=1,NDM
      XL(J,I)=0.D0
233   CONTINUE
      GO TO 238
235   CONTINUE
      NEL=I
      DO 236 J=1,NDM
      XL(J,I)=X(J,II)
236   CONTINUE
238   CONTINUE
C---- ENTER ELEM LIBRARY WITH ISW=12 TO READ TRACTIONS ALONG EDGES OF
C---- ELEM# N AND COMPUTE EQUIVALENT NODAL LOADS.
      CALL ELMLIB(D(1,MA),DUM,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),
     1 DUM,S,PDIST(1,KDIST),DUM,DUM,NDF,NDM,NST,NSDM,NQDM,12)
C---- WRITE EQUIVALENT NODAL LOADS FOR DISTRIBUTED LOAD# KDIST
      IF(PRT)
     1 WRITE(6,2015) KDIST,IDIST(1,KDIST),IDIST(2,KDIST),MA,IEL,
     2 (DOF,KDF,KDF=1,NDF)
      DO 239 I=1,NEN
      IF(IX(I,N).EQ.0) GO TO 239
      JS=(I-1)*NDF+1
      JS1=I*NDF
      DO 232 KS=JS,JS1
      IF(PDIST(KS,KDIST).NE.0.D0) GO TO 234
232   CONTINUE
      GO TO 239
234   IF(PRT) WRITE(6,2016) IX(I,N),(PDIST(KS,KDIST),KS=JS,JS1)
239   CONTINUE
231   CONTINUE
      GO TO 10
C
C---- MACRO 'DOMA'
C---- DOMAIN J-INTEGRAL
24    CONTINUE
      CALL DOMCOR
      GO TO 10
C
C---- MACRO 'SET '
C---- READ ELEMENT OR NODE SETS FOR PRINT OUT
C---- MAXIMUM OF 9 NODE OR ELEM SETS ARE PERMITTED
C---- EACH NODE/ELEM SET MAY CONTAIN A MAXIMUM OF 200 NODES/ELEMENTS
C---- ADAPTED FROM FEAP VERSION OF TOSHIO NAKAMURA.
25    CONTINUE
      READ(5,1002) NPS,NPST,IGEN,IGENT
C
C---- NPS  = NODE/ELEM SET # (.LE.9)
C---- NPST = # OF NODES/ELEM IN THIS SET (.LE.200)
C---- IGEN = GENERATION INCREMENT USED TO GENERATE NODES/ELEMS IN THIS SET
C---- IGENT= LAST NODE/ELEM # TO BE OBTAINED BY GENERATION.
C---- NOTE 1) FOR EACH NODE/ELEM SET ONLY 1 GENERATION IS POSSIBLE
C---- NOTE 2) IF IGENT = 0, GENERATION TO MAX # (NPST) OF NODES/ELEMS WILL
C----         BE CARRIED OUT (IF IGEN.GT.0)
C----         IF IGENT.NE.0, THEN THE REMAINING (NPST-IGENT) NODES/ELEMS
C----         MUST BE INPUT WITHOUT GENERATION.
C
      IF ((NPS.GT.9).OR.(NPST.GT.200)) THEN
      WRITE (6,2017) NPS,NPST
      STOP
      END IF
      NSET(NPS,1) = NPST
      IF(IGEN.EQ.0) THEN
      READ(5,1003) (NSET(NPS,NN),NN=2,NPST+1)
      ELSE
      READ(5,1003) NSET(NPS,2)
      IF(IGENT.EQ.0) IGENT=NPST
      DO 252 I=2,IGENT
252   NSET(NPS,I+1)=NSET(NPS,I)+IGEN
      IF(IGENT.LT.NPST) READ(5,1003) (NSET(NPS,NN+1),NN=IGENT+1,NPST)
      END IF
      IF(PRT) WRITE(6,2018) NPS,(NSET(NPS,NN),NN=1,NPST+1)
      GO TO 10
C
C---- MACRO 'CONT'
26    CONTINUE
C---- READ NODE NUMBERS OF TENTATIVE SLAVE AND MASTER NODES ON SLIDELINE
      VOID = 1.50
      DO 261 J=1,50
      IC(J) = 0
      IM(2*J-1) = 0
      IM(2*J) = 0
      NDEST(J) = 0
      MSEG(J) = 0
      IAB(1,J) = 0  
      IAB(2,J) = 0
      BET(J) = VOID
      ANG(J) = 0.D0
      DO 262 K=1,NDF
      DUSL(K,J) = 0.D0
      DEL(K,J) = 0.D0
      CP(K,J)  = 0.D0
262   CONTINUE
261   CONTINUE
      READ(5,1005) NM,ALP,SMU,DMU
      WRITE(6,4000)NM,ALP,SMU,DMU
4000  FORMAT('NO. OF MASTER NODES =',I3,/,'PENALTY NUMBER=',E13.5,/,
     1 'COEFF. OF STICKING FRICTION=',E13.5,/,'COEFF. OF ', 
     1 'SLIDING  FRICTION=',E13.5)
      READ(5,1003) (IC(J),J=1,NC),(IM(J),J=1,NM)
      WRITE(6,4001)
4001  FORMAT('        NODE NO.S OF SLAVE NODES         ')
      WRITE(6,1003) (IC(J),J=1,NC)
      WRITE(6,4002)
4002  FORMAT('       NODE NO.S OF MASTER NODES         ')
      WRITE(6,1003) (IM(J),J=1,NM)
      GO TO 10
C
C---- MACRO 'INTF'
C---- INPUT DATA FOR INTERFACE MODELLING (COHESIVE ZONE WITH LUMPED
C---- STIFFNESS MATRIX)
27    CONTINUE
c     Read NINTND, (INTND(1,I),I=1,NINTND+1),ALPAA,DLTAA,SIGMXX
c     For an open loop enter junk value for INTND(1,NINTND+1) and for a 
c     closed loop enter INTND(NINTND+1)=INTND(1) in the input data file.
c     INTND(1,*) are the interface nodes of matrix material.
c     INTND(2,*) are the intf nodes of inclusion, initially both INTND(1,*) 
c     and INTND(2,*) are having same coordinates.
      READ(5,1004) NINTND,ALPAA,DLTAA,SIGMXX
      READ(5,1003) (INTND(1,I),I=1,NINTND+1)      
      READ(5,1003) (INTND(2,I),I=1,NINTND+1)      
      READ(5,1003) (ISTAT(I),I=1,NINTND)      
      WRITE(6,2019) NINTND,ALPAA,DLTAA,SIGMXX
      WRITE(6,2020) 
      WRITE(6,1003) (INTND(1,I),I=1,NINTND+1)
      WRITE(6,2129) 
      WRITE(6,1003) (INTND(2,I),I=1,NINTND+1)
      WRITE(6,*)'INTF NODE STATUS: 0 --> BONDED, NON-ZERO --> DEBONDED'
      WRITE(6,1003) (ISTAT(I),I=1,NINTND)
C____ ******** INTERFACE FORMULATIONS ********
c____ This part calculates the segment lengths (xln(*)), the segment length
c---- associated with each node (xlntnd(*)) and the angle Theta (angint(*))
c---- (by the tangent with the x-axis).
c____ Read nintnd, (intnd(1,i),i=1,nintnd+1),alpaa, dltaa, sigmxx
c____ For an open loop enter junk value for INTND(1,NINTND+1) and for a
c     closed loop. 
c____ INTND(1,NINTND+1)=INTND(1,1) in the input data file.
      NINTT = NINTND
      IF(INTND(1,NINTND+1) .NE. INTND(1,1))NINTT=NINTND-1
      DO 271 I=1,NINTT
      DLX = X(1,INTND(1,I+1)) - X(1,INTND(1,I))
      DYY = X(2,INTND(1,I+1)) - X(2,INTND(1,I))
      XLN(I) = DSQRT(DLX**2 + DYY**2)
271   CONTINUE
c____ Calculate the segment length associated with each node
      DO 272 I=2,NINTT
      XLNTND(I) = (XLN(I-1) + XLN(I))*0.50
272   CONTINUE
      IF(NINTT .EQ. NINTND) THEN
      XLNTND(1) = (XLN(1) + XLN(NINTND))*0.50
      ELSE
      XLNTND(1) = 0.50*XLN(1)
      XLNTND(NINTND) = 0.50*XLN(NINTT)
      ENDIF
c____ Calculate the angle at each node
      DO 273 I=2,NINTT
      DXX = X(1,INTND(1,I)) - X(1,INTND(1,I+1))
      DYY = X(2,INTND(1,I)) - X(2,INTND(1,I+1))
      TETAA1 = DATAN2(DYY,DXX)
      IF (TETAA1 .LT. 0.D0) TETAA1=2.D0*PI + TETAA1
      DLX = X(1,INTND(1,I+1)) - X(1,INTND(1,I-1))
      DYY = X(2,INTND(1,I+1)) - X(2,INTND(1,I-1))
      DL = DSQRT(DLX**2 + DYY**2)
      TETAA2 =(XLN(I-1)**2+XLN(I)**2-DL**2)/(2.D0*XLN(I-1)*XLN(I))
      TETAA2 =0.5*DACOS(TETAA2)
      ANGINT(I) =TETAA1 - TETAA2 -PI*0.50
      IF(ANGINT(I) .GE. 2.*PI) ANGINT(I)=ANGINT(I)-2.*PI
273   CONTINUE
      IF(NINTT .EQ. NINTND) THEN
      DXX = X(1,INTND(1,1)) - X(1,INTND(1,2))
      DYY = X(2,INTND(1,1)) - X(2,INTND(1,2))
      TETAA1 = DATAN2(DYY,DXX)
      IF (TETAA1 .LT. 0.D0) TETAA1=2.D0*PI + TETAA1
      DLX = X(1,INTND(1,NINTND)) - X(1,INTND(1,2))
      DYY = X(2,INTND(1,NINTND)) - X(2,INTND(1,2))
      DL = DSQRT(DLX**2 + DYY**2)
      TETAA2=(XLN(NINTND)**2+XLN(1)**2-DL**2)/(2.D0*XLN(NINTND)*XLN(1))
      TETAA2 =0.5*DACOS(TETAA2)
      ANGINT(1) =TETAA1 - TETAA2 -PI*0.50
      IF(ANGINT(1) .GE. 2.*PI) ANGINT(1)=ANGINT(1)-2.*PI
      ELSE
      DXX = X(1,INTND(1,1)) - X(1,INTND(1,2))
      DYY = X(2,INTND(1,1)) - X(2,INTND(1,2))
      TETAA1 = DATAN2(DYY,DXX)
      IF (TETAA1 .LT. 0.00001) TETAA1=2.D0*PI + TETAA1
      ANGINT(1) =TETAA1 + PI
      IF(ANGINT(1) .GE. 2.*PI) ANGINT(1)=ANGINT(1)-2.*PI
      DXX = X(1,INTND(1,NINTT)) - X(1,INTND(1,NINTND))
      DYY = X(2,INTND(1,NINTT)) - X(2,INTND(1,NINTND))
      TETAA1 = DATAN2(DYY,DXX)
      IF (TETAA1 .LT. 0.D0) TETAA1=2.D0*PI + TETAA1
      ANGINT(NINTND) =TETAA1 + PI
      IF(ANGINT(NINTND) .GE. 2.*PI) 
     1 ANGINT(NINTND)=ANGINT(NINTND)-2.*PI
      ENDIF
      WRITE(6,*)'INTERFACE NODES, ANGLES AND',
     1 'ASSOCIATED SEGMENT LENGTHS:'
C     DO 7777 I=1,NINTND
C     ANGINT(I) = 0.0
C7777  CONTINUE
      WRITE(6,2731) (INTND(1,I), (180.D0/PI)*ANGINT(I),
     1 XLNTND(I),I=1,NINTND)
2731  FORMAT(2X,I5,2X,E13.5,2X,E13.5)
      GO TO 10
C
C---- MACRO 'SURF'
C---- INPUT DATA FOR MODELLING WEDGE INDENTATION
28    CONTINUE
      READ(5,2872)NSURF,THTSRF,PENSRF,MSURF
2872  FORMAT(I5,2F10.0,I5)
      WRITE(6,2874)NSURF,THTSRF,PENSRF,MSURF
2874  FORMAT(1X,'NO. OF MASTER NODES = ',I5,/1X,
     1 'SEMI-ANGLE OF WEDGE = ',F10.2,/1X,'PENALTY NO. = ',
     2 E13.5/1X,'NO. OF MASTER NODES ON INCLINED SURFACE = ',I5)
      READ(5,2873)(ISURF(I),I=1,NSURF)
2873  FORMAT(15I5)
      WRITE(6,2875)(ISURF(I),I=1,NSURF)
2875  FORMAT(1X,'MASTER NODES ARE ',/1X,15I5)
      IF(MSURF .GT. 0) THEN
      DO 2876 I=2,MSURF
         NDSURF(I)=1
2876  CONTINUE
      ENDIF
      GO TO 10
C
C---- FORMATS
 1000 FORMAT(A4,75X,A1)
1002  FORMAT(16I5)
1003  FORMAT(15I5)
1004  FORMAT(I5,3F10.0)
1005  FORMAT(I5,F10.0,F10.0,F10.0,F10.0)
2003  FORMAT(/5X,12HMATERIAL SET,I3,17H FOR ELEMENT TYPE,I2,5X,//
     1   10X,49HDEGREE OF FREEDOM ASSIGNMENTS    LOCAL    GLOBAL   /
     2   42X, 6HNUMBER, 4X, 6HNUMBER/(36X,2I10))
2004  FORMAT(A1,20A4//5X,19HMATERIAL PROPERTIES)
2005  FORMAT(A1,20A4//5X,17HNODAL FORCE/DISPL//6X,4HNODE,9(I7,A4,A2))
2006  FORMAT(I10,9E13.3)
2010  FORMAT(5X,'** FATAL ERROR ** ATTEMPT TO INPUT',I6,'DISTRIBUTED',
     1 ' LOAD CARDS USING MACRO DIST')
2011  FORMAT(5X,'** FATAL ERROR ** ATTEMPT TO INPUT DISTRIBUTED LOADS'
     1 ,' FOR ELEMENT',I6,' TERMINATION OCCURED IN MACRO DIST')
2012  FORMAT(5X,'** FATAL ERROR ** ATTEMPT TO REFER TO PROP TABLE #',
     1 I6,' WHEN READING DISTRIBUTED LOADS IN MACRO DIST')
2014  FORMAT(//5X,'EQUIVALENT NODAL FORCES FOR DISTRIBUTED EDGE LOADS')
2015  FORMAT(//2X,'DIST LOAD#=',I5,3X,'ELEM#=',I5,3X,'PROP TABLE#=',I5
     1 ,3X,'MATL GROUP#=',I5,3X,'ELEM TYPE#=',I5/4X,'NODE #',6(8X,
     2 A3,I2))
2016  FORMAT(I10,6E13.5)
2017  FORMAT(/5X,'** FATAL ERROR ** NUMBER OF THE NODE OR ELEM SET',I4,
     1 ' EXCEEDS 9 OR'/5X,'NUMBER OF NODES OR ELEMS IN ABOVE SET',I4,
     2 ' EXCEEDS 200. TERMINATION OCCURED IN MACRO SET')
2018  FORMAT(/5X,' SET NUMBER',I3,' HAS TOTAL OF',I4,' ELEMENTS OR ',
     1 'NODES'/14(/15I6))
2019  FORMAT(//2X,'** NO. OF INTERFACE NODES =',I5,
     1        /2X,'** VALUE OF ALPAA         =',E13.5,
     2        /2X,'** VALUE OF DELTAA        =',E13.5,
     3        /2X,'** VALUE OF SIGMAMAX      =',E13.5)
2020  FORMAT(//2X,'INTERFACE NODE NUMBERS [MATRX MATERIAL]')
2129  FORMAT(//2X,'INTERFACE NODE NUMBERS [INCLUSION]')
3003  FORMAT(5X,'**WARNING 01** ELEMENT CONNECTIONS NECESSARY TO USE BLO
     1K IN MACRO PROGRAM')
      END
      SUBROUTINE GRID(IDL,IX,NEN1,PRT,ERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- READ AND/OR GENERATE ELEMENT CONNECTIVITY
C
      LOGICAL PRT,ERR
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
      DIMENSION IX(NEN1,*),IDL(*)
C
      IF (PRT) WRITE(6,2001) O,HEAD,(K,K=1,NEN)
      L = 0
      KOUNT = 1
      DO 206 N = 1,NUMEL
      IF (IX(NEN1,N).NE.0) GO TO 206
      IF (KOUNT.GE.50) THEN
      IF (PRT) WRITE(6,2001) O,HEAD,(K,K=1,NEN)
      KOUNT = 1
      ELSE
      KOUNT = KOUNT + 1
      END IF
      IF(L-N) 200,202,203
200   CONTINUE
      READ(5,1001) L,LK,(IDL(K),K=1,NEN),LX
      IF(L.EQ.0) L = NUMEL+1
      IF(LX.EQ.0) LX=1
      IF(L-N) 201,202,203
201   WRITE(6,3001) L,N
      ERR = .TRUE.
      GO TO 206
202   NX = LX
      DO 207 K = 1,NEN
      IF(IDL(K).GT.NUMNP.OR.IDL(K).LT.0) GO TO 208
207   IX(K,L) = IDL(K)
      IX(NEN1,L) = LK
      GO TO 205
203   IX(NEN1,N) = IX(NEN1,N-1)
      DO 204 K = 1,NEN
      IX(K,N) = IX(K,N-1) + NX
      IF(IX(K,N-1).EQ.0) IX(K,N) = 0
204   IF(IX(K,N).GT.NUMNP.OR.IX(K,N).LT.0) GO TO 208
205   IF(PRT) WRITE(6,2002) N,IX(NEN1,N),(IX(K,N),K=1,NEN)
      GO TO 206
208   WRITE(6,3002) N
      ERR = .TRUE.
206   CONTINUE
      RETURN
 1001 FORMAT(16I5)
2001  FORMAT(A1,20A4//5X,8HELEMENTS//3X,7HELEMENT,2X,8HMATERIAL,
     1   14(I3,5H NODE)/(20X,14(I3,5H NODE)))
2002  FORMAT(2I10,14I8/(20X,14I8))
3001  FORMAT(5X,20H**ERROR 03** ELEMENT,I5,22H APPEARS AFTER ELEMENT,I5)
3002  FORMAT(5X,'**ERROR 04** ELEMENT',I5,' HAS ILLEGAL NODES')
      END
      SUBROUTINE BOUND(IDL,ID,NDF,PRT,ERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- READ AND/OR GENERATE BOUNDARY CONDITIONS
C
      LOGICAL PRT,ERR
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
      DIMENSION ID(NDF,*),IDL(*),BC(2)
      DATA BC/4H B.C,2H. /
C
      IF(PRT) WRITE(6,2000) O,HEAD,(I,BC,I=1,NDF)
      N = 0
      NG = 0
402   L = N
      LG = NG
      READ(5,1001) N,NG,(IDL(I),I=1,NDF)
      IF(N.LE.0.OR.N.GT.NUMNP) GO TO 60
      DO 51 I = 1,NDF
      ID(I,N) = IDL(I)
51    IF(L.NE.0.AND.IDL(I).EQ.0.AND.ID(I,L).LT.0) ID(I,N) = -1
      IF((N-L).LT.0) LG = - LG
52    L = L + LG
      IF((N-L)*LG.LE.0) GO TO 402
      DO 53 I = 1,NDF
53    IF(ID(I,L-LG).LT.0) ID(I,L) = -1
      GO TO 52
60    DO 58 N = 1,NUMNP
      DO 56 I = 1,NDF
56    IF(ID(I,N).NE.0) GO TO 57
      GO TO 58
57    IF(PRT) WRITE(6,2007) N,(ID(I,N),I=1,NDF)
58    CONTINUE
      RETURN
 1001 FORMAT(16I5)
2000  FORMAT(A1,20A4//5X,17HNODAL B.C.       //6X,4HNODE,9(I7,A4,A2)/1X)
2007  FORMAT(I10,9I13)
      END
      SUBROUTINE GENEL(NDM,NQUAD,NL,TITLE,Y,PRT,ERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- GENERATE ELEMENT DATA ARRAYS BY LINEAR INTERPOLATION
C
      CHARACTER*18 TITLE
      CHARACTER*10 COMP
      LOGICAL PRT,ERR
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQD1
      DIMENSION Y(NL,*)
      DATA COMP/' COMPONENT'/
C---- INITIALIZE PRINTOUT COUNTER, ELMT NUMBER AND GENERATION
C---- INCREMENT
      MCT = 0
      N = 0
      NG = 0
C---- SET FIRST GENERATED ELEMENT TO ELMT NUMBER READ IN PREVIOUS
C---- CARD
  102 L = N
C---- LIKEWISE WITH GENERATION INCREMENT
      LG = NG
C---- READ NEW CARD
      READ(5,1000) N,NG
C---- CHECK FOR END OF INPUT PROCESS
      IF (N.LE.0) GO TO 109
      K0 = 0
      DO 103 K2 = 1,NQUAD
      READ(5,1005) (Y(K0 + K1,N),K1=1,NDM)
      K0 = K0 + NDM
  103 CONTINUE
      NG = IABS(NG)
      IF (LG) 104,102,104
C---- NONZERO GENERATION INCREMENT
  104 CONTINUE
      IF (N.LT.L) LG = - LG
C---- COMPUTE DATA INCREMENTS BY LINEAR INTERPOLATION
      LI = (IABS(N-L+LG) - 1)/IABS(LG)
      DO 105 I = 1,NL
  105 Y(I,N) = (Y(I,N) - Y(I,L))/LI
C---- MARCH THROUGH GENERATED ELEMENTS
  106 L = L + LG
      IF ((N-L)*LG.LT.0) GO TO 102
C---- CHECK IF ELEMENT NUMBER LIES WITHIN ADMISSIBLE RANGE
      IF ((L.LE.0).OR.(L.GT.NUMEL)) GO TO 900
C---- INCREMENT DATA
      DO 107 I = 1,NL
  107 Y(I,L) = Y(I,L-LG) + Y(I,N)
      GO TO 106
C---- PRINTOUT GENERATED DATA
  109 CONTINUE
      IF(.NOT.PRT) RETURN
      DO 113 J = 1,NUMEL
      MCT = MCT - 1
      IF (MCT.GT.0) GO TO 112
      MCT = 10
      WRITE(6,2000) O,HEAD,TITLE,(COMP,K,K=1,NDM)
  112 CONTINUE
      WRITE(6,2010) J,(Y(K1,J),K1=1,NDM)
      K0 = NDM
      DO 115 K2 = 2,NQUAD
      WRITE(6,2015) (Y(K0 + K1,J),K1=1,NDM)
      K0 = K0 + NDM
  115 CONTINUE
  113 CONTINUE
      RETURN
  900 CONTINUE
      WRITE(6,3000) L
 3000 FORMAT(///5X,'** FATAL ERROR 06 ** ',
     1 'ATTEMPT TO GENERATE ELEMENT ',I5)
      ERR = .TRUE.
      GO TO 102
 1000 FORMAT(16I5)
 1005 FORMAT(8F10.0)
 2000 FORMAT(A1,20A4//5X,A18//5X,'ELMT ',8(A10,I2),
     1 50(/10X,8(A10,I2)))
 2010 FORMAT(/I10,8E12.4,50(/10X,8E12.4))
 2015 FORMAT( 10X,8E12.4,50(/10X,8E12.4))
      END
      SUBROUTINE POLAR(X,NDM,PRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---- CONVERT POLAR TO CARTESIAN COORDINATES
      LOGICAL PRT
      DIMENSION X(NDM,*)
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      IF(NDM.EQ.1) RETURN
      MCT = 0
      TH = DATAN(1.0D0)/45.D0
100   CONTINUE
      READ(5,1000) NI,NE,INC,X0,Y0
      IF(NI.LE.0) RETURN
      IF(NI.GT.NUMNP.OR.NE.GT.NUMNP) GO TO 300
      INC = MAX0(IABS(INC),1)
      IF((NE-NI).LT.0) INC = - INC
      IF(NE.EQ.0) NE = NI
      N = NI
200   R = X(1,N)
      X(1,N) = X0 + R*DCOS(X(2,N)*TH)
      X(2,N) = Y0 + R*DSIN(X(2,N)*TH)
      IF(MCT.GT.0) GO TO 250
      IF(PRT) WRITE(6,2000) O,HEAD,X0,Y0,(I,I=1,NDM)
      MCT = 50
250   IF(PRT) WRITE(6,2001) N,(X(I,N),I=1,NDM)
      MCT = MCT - 1
      N = N + INC
      IF((NE-N)*INC.GE.0) GO TO 200
      IF(MOD(NE-NI,INC).EQ.0) GO TO 100
      NI = NE
      N = NE
      GO TO 200
C---- ERROR
300   WRITE(6,3000) NI,NE
      STOP
C---- FORMATS
1000  FORMAT(3I5,2F10.0)
2000  FORMAT(A1,20A4//'   CARTESIAN COORDINATES COMPUTED FROM POLAR INPU
     1T WITH X0 = ',E12.4,'    Y0 = ',E12.4/6X,'NODE',6(I7,'-COORD'))
2001  FORMAT(I10,6F13.4)
3000  FORMAT('  **FATAL ERROR 16** ATTEMPT TO CONVERT NDES NI = ',I6,
     1 ' TO NE = ',I6)
      END
C--------------------------------------------------------------------
C TEST SEGMENT----FOR PARABOLIC INTERPOLATION FOR NODE COORDINATE
      SUBROUTINE GENCOR(NDM,X,CD,PRT,ERR,PRTZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- GENERATE REAL DATA ARRAYS BY LINEAR INTERPOLATION OR
C---- PARABOLIC INTERPOLATION...TEST VERSION
C---- G. RAVICHANDRAN 10/8/86 AT CIT
C---- TO BE USED ONLY FOR GENERATING NODE COORDINATE
C---- DATA FOR QUARTER POINT ELEMENTS
C
      LOGICAL PRT,ERR,PCOMP,PRTZ
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      DIMENSION X(NDM,*),XL(3),CD(*),XG(3)
      DATA BL/4HBLAN/
      MCT = 0
      N = 0
      NG = 0
102   L = N
      LG = NG
      READ(5,1000) N,NG,LP,XL
C---- LP=1 (LINEAR)/LP=2 (PARABOLIC) INTERPOLATION FOR NODE COOR
      IF(N.GT.NUMNP) WRITE(6,3001) N,(CD(I),I=1,3)
      IF(N.LE.0.OR.N.GT.NUMNP) GO TO 109
      DO 103 I = 1,NDM
103   X(I,N) = XL(I)
C---- SET PATH FOR INTERPOLATION
      IF(LP.EQ.2) GO TO 201
C
      IF(LG) 104,102,104
104   CONTINUE
      IF((N-L).LT.0) LG = - LG
      LI =(IABS(N-L+LG)-1)/IABS(LG)
      DO 105 I = 1,NDM
105   XL(I) = (X(I,N)-X(I,L))/LI
106   L = L + LG
      IF((N-L)*LG.LE.0) GO TO 102
      IF(L.LE.0.OR.L.GT.NUMNP) GO TO 108
      DO 107 I = 1,NDM
107   X(I,L) = X(I,L-LG) + XL(I)
      GO TO 106
108   WRITE(6,3000) L,(CD(I),I=1,3)
      ERR = .TRUE.
      GO TO 102
109   IF(.NOT.PRT) RETURN
      DO 113 J = 1,NUMNP
      IF(PRTZ) GO TO 111
      DO 110 L = 1,NDM
      IF(X(L,J).NE.0.0D+00) GO TO 111
110   CONTINUE
      GO TO 113
111   MCT = MCT - 1
      IF(MCT.GT.0) GO TO 112
      MCT = 50
      WRITE(6,2000) O,HEAD,(CD(L),L=1,3),(L,CD(1),CD(2),L=1,NDM)
112   IF(PCOMP(X(1,J),BL)) WRITE(6,2008) J
      IF(.NOT.PCOMP(X(1,J),BL)) WRITE(6,2009) J,(X(L,J),L=1,NDM)
113   CONTINUE
      RETURN
C---- FOLLOWING SEGMENT DOES PARABOLIC INTERPOLATION
C
201   IF(LG) 204,102,204
204   CONTINUE
      IF((N-L).LT.0) LG = - LG
      LI =(IABS(N-L+LG)-1)/IABS(LG)
      DO 205 I = 1,NDM
205   XG(I) = (X(I,N)-X(I,L))
      ML=L
      L1=0
206   L1 = L1 + 1
      L=L+LG
      IF((N-L)*LG.LE.0) GO TO 102
      IF(L.LE.0.OR.L.GT.NUMNP) GO TO 208
      AL1=DFLOAT(L1)/DFLOAT(LI)
      AL=AL1*AL1
      DO 207 I = 1,NDM
207   X(I,L) = X(I,ML) + AL*XG(I)
      GO TO 206
208   WRITE(6,3000) L,(CD(I),I=1,3)
      ERR = .TRUE.
      GO TO 102
209   IF(.NOT.PRT) RETURN
      DO 213 J = 1,NUMNP
      IF(PRTZ) GO TO 211
      DO 210 L = 1,NDM
      IF(X(L,J).NE.0.0D+00) GO TO 211
210   CONTINUE
      GO TO 213
211   MCT = MCT - 1
      IF(MCT.GT.0) GO TO 212
      MCT = 50
      WRITE(6,2000) O,HEAD,(CD(L),L=1,3),(L,CD(1),CD(2),L=1,NDM)
212   IF(PCOMP(X(1,J),BL)) WRITE(6,2008) J
      IF(.NOT.PCOMP(X(1,J),BL)) WRITE(6,2009) J,(X(L,J),L=1,NDM)
213   CONTINUE
      RETURN
 1000 FORMAT(3I5,3E15.6)
2000  FORMAT(A1,20A4//5X, 5HNODAL,3A4//6X,4HNODE,9(I7,A4,A2))
2008  FORMAT(I10,32H HAS NOT BEEN INPUT OR GENERATED)
2009  FORMAT(I10,9F13.4)
3000  FORMAT(5X,43H**FATAL ERROR 02** ATTEMPT TO GENERATE NODE,I5,3H IN
     1  ,3A4)
3001  FORMAT(5X,'**FATAL ERROR 05** ATTEMPT TO INPUT NODE',I5,', TERMINA
     1TE INPUT OF NODES IN',3A4)
      END
C--------------------------------------------------------------------
      SUBROUTINE GENVEC(NDM,X,CD,PRT,ERR,PRTZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- GENERATE REAL DATA ARRAYS BY LINEAR INTERPOLATION
C
      LOGICAL PRT,ERR,PCOMP,PRTZ
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      DIMENSION X(NDM,*),XL(6),CD(*)
      DATA BL/4HBLAN/
      MCT = 0
      N = 0
      NG = 0
102   L = N
      LG = NG
      READ(5,1000) N,NG,XL
      IF(N.GT.NUMNP) WRITE(6,3001) N,(CD(I),I=1,3)
      IF(N.LE.0.OR.N.GT.NUMNP) GO TO 109
      DO 103 I = 1,NDM
103   X(I,N) = XL(I)
      IF(LG) 104,102,104
104   CONTINUE
      IF((N-L).LT.0) LG = - LG
      LI =(IABS(N-L+LG)-1)/IABS(LG)
      DO 105 I = 1,NDM
105   XL(I) = (X(I,N)-X(I,L))/LI
106   L = L + LG
      IF((N-L)*LG.LE.0) GO TO 102
      IF(L.LE.0.OR.L.GT.NUMNP) GO TO 108
      DO 107 I = 1,NDM
107   X(I,L) = X(I,L-LG) + XL(I)
      GO TO 106
108   WRITE(6,3000) L,(CD(I),I=1,3)
      ERR = .TRUE.
      GO TO 102
109   IF(.NOT.PRT) RETURN
      DO 113 J = 1,NUMNP
      IF(PRTZ) GO TO 111
      DO 110 L = 1,NDM
      IF(X(L,J).NE.0.0D+00) GO TO 111
110   CONTINUE
      GO TO 113
111   MCT = MCT - 1
      IF(MCT.GT.0) GO TO 112
      MCT = 50
      WRITE(6,2000) O,HEAD,(CD(L),L=1,3),(L,CD(1),CD(2),L=1,NDM)
112   IF(PCOMP(X(1,J),BL)) WRITE(6,2008) J
      IF(.NOT.PCOMP(X(1,J),BL)) WRITE(6,2009) J,(X(L,J),L=1,NDM)
113   CONTINUE
      RETURN
 1000 FORMAT(2I5,6F10.0)
2000  FORMAT(A1,20A4//5X, 5HNODAL,3A4//6X,4HNODE,9(I7,A4,A2))
2008  FORMAT(I10,32H HAS NOT BEEN INPUT OR GENERATED)
2009  FORMAT(I10,9F13.4)
3000  FORMAT(5X,43H**FATAL ERROR 02** ATTEMPT TO GENERATE NODE,I5,3H IN
     1  ,3A4)
3001  FORMAT(5X,'**FATAL ERROR 05** ATTEMPT TO INPUT NODE',I5,', TERMINA
     1TE INPUT OF NODES IN',3A4)
      END
C-------------------------------------------------------------------
C---- TEST SEGMENT FOR GENERATING FORCE AND IDPROP VALUES
      SUBROUTINE GENFRC(NDM,X,ID,CD,PD,PRT,ERR,PRTZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- GENERATE REAL DATA ARRAYS BY LINEAR INTERPOLATION
C
      LOGICAL PRT,ERR,PCOMP,PRTZ
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      DIMENSION X(NDM,*),XL(6),CD(*)
      DIMENSION ID(NDM,*),IDL(6),PD(*)
      DATA BL/4HBLAN/
      MCT = 0
      N = 0
      NG = 0
102   L = N
      LG = NG
      READ(5,1000) N,NG,IDL
      READ(5,1001) XL
      IF(N.GT.NUMNP) WRITE(6,3001) N,(CD(I),I=1,3)
      IF(N.LE.0.OR.N.GT.NUMNP) GO TO 109
      DO 103 I = 1,NDM
      ID(I,N)=IDL(I)
103   X(I,N) = XL(I)
      IF(LG) 104,102,104
104   CONTINUE
      IF((N-L).LT.0) LG = - LG
      LI =(IABS(N-L+LG)-1)/IABS(LG)
      DO 105 I = 1,NDM
      IDL(I)=ID(I,L)
105   XL(I) = (X(I,N)-X(I,L))/LI
106   L = L + LG
      IF((N-L)*LG.LE.0) GO TO 102
      IF(L.LE.0.OR.L.GT.NUMNP) GO TO 108
      DO 107 I = 1,NDM
C---- IDPROP VALUES FOR NODES BEING GENERATED ARE SET EQUAL TO THE
C---- CORRESPONDING VALUES ON THE MASTER CARD (CARD1).
      ID(I,L)=IDL(I)
107   X(I,L) = X(I,L-LG) + XL(I)
      GO TO 106
108   WRITE(6,3000) L,(CD(I),I=1,3)
      ERR = .TRUE.
      GO TO 102
109   IF(.NOT.PRT) RETURN
      DO 113 J = 1,NUMNP
      IF(PRTZ) GO TO 111
      DO 110 L = 1,NDM
      IF(X(L,J).NE.0.0D+00) GO TO 111
110   CONTINUE
      GO TO 113
111   MCT = MCT - 1
      IF(MCT.GT.0) GO TO 112
      MCT = 50
      WRITE(6,2000) O,HEAD,(CD(L),L=1,3),(L,CD(1),CD(2),L=1,NDM)
112   IF(PCOMP(X(1,J),BL)) WRITE(6,2008) J
      IF(.NOT.PCOMP(X(1,J),BL)) WRITE(6,2009) J,(X(L,J),L=1,NDM)
113   CONTINUE
      MCT=0
      DO 114 J=1,NUMNP
      IF(PRTZ) GO TO 115
      DO 121 L=1,NDM
      IF((ID(L,J).LT.1).OR.(ID(L,J).GT.7)) THEN
      WRITE(6,3002) J,L,(PD(I),I=1,3)
      ERR=.TRUE.
      END IF
121   CONTINUE
      DO 116 L=1,NDM
      IF(X(L,J).NE.0.0D0) GO TO 115
116   CONTINUE
      GO TO 114
115   MCT=MCT-1
      IF(MCT.GT.0) GO TO 117
      MCT=50
      WRITE(6,2000) O,HEAD,(PD(L),L=1,3),(L,PD(1),PD(2),L=1,NDM)
117   CONTINUE
      WRITE(6,2010) J,(ID(L,J),L=1,NDM)
114   CONTINUE
      RETURN
 1000 FORMAT(8I5)
 1001 FORMAT(6F10.0)
2000  FORMAT(A1,20A4//5X, 5HNODAL,3A4//6X,4HNODE,9(I7,A4,A2))
2008  FORMAT(I10,32H HAS NOT BEEN INPUT OR GENERATED)
2009  FORMAT(I10,9F13.4)
2010  FORMAT(I10,9I13)
3000  FORMAT(5X,43H**FATAL ERROR 02** ATTEMPT TO GENERATE NODE,I5,3H IN
     1  ,3A4)
3001  FORMAT(5X,'**FATAL ERROR 05** ATTEMPT TO INPUT NODE',I5,', TERMINA
     1TE INPUT OF NODES IN',3A4)
3002  FORMAT(5X,'**FATAL ERROR 5A** PROPID FOR NODE',I5,' AND DOF',I5,
     1 ' EXCEEDS 7 OR IS LESS THAN 1'/5X,'**TERMINATE INPUT IN',3A4)
      END
      SUBROUTINE BLKGEN(NDM,NDF,NEL,NEL1,X,IX,PRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- THIS SUBROUTINE GENERATES A WHOLE PATCH OF ELEMENTS AND NODES.
C---- THE INPUT IS AS FOLLOWS:
C
C---- NN = NUMBER OF NODES OF MACROELEMENT TO BE REFINED (4-9)
C---- NR = NUMBER OF SUBINTERVALS IN THE R-DIRECTION
C---- NS = NUMBER OF SUBINTERVALS IN THE S-DIRECTION
C---- NI = FIRST GLOBAL NODE TO BE GENERATED
C---- NE = FIRST GLOBAL ELMT TO BE GENERATED
C---- MA = MATERIAL GROUP OF THE ELEMENTS
C---- NGAP = GAP IN NODE NUMBERING  IN R-DIR
C
C---- THEN, FOR N = 1,NN, PROVIDE:
C
C---- L  = LOCAL NODE NUMBER (IN THE MACROELEMENT)
C---- R  = 1-COORDINATE
C---- S  = 2-COORDINATE
C
      LOGICAL PRT
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      DIMENSION X(NDM,*),IX(NEL1,*),XL(2,9),IXL(9),SHP(3,9)
      DATA XH/6H COORD/
C---- INPUT DATA
      READ(5,1000) NN,NR,NS,NI,NE,MA,NGAP
      NR = MAX0(NR,1)
      NS = MAX0(NS,1)
      NI = MAX0(NI,1)
      MA = MAX0(MA,1)
      IF (PRT) THEN
      WRITE(6,2000) O,HEAD,NR,NS,NI,NE,MA
      IF(NE.EQ.0) WRITE(6,2010)
      WRITE(6,2002) (I,XH,I=1,NDM)
      END IF
      DO 10 N = 1,9
      DO 10 J = 1,NDM
      XL(J,N) = 0.0
10    IXL(N) = 0
      NM = 0
C---- SET UP LOCAL CONNECTION ARRAY AND NODAL COORDINATES
      DO 20 N = 1,NN
      READ(5,1001) L,R,S
      IF(L.EQ.0) L = N
      NM = MAX0(NM,L)
      IXL(L) = L
      XL(1,L) = R
      XL(2,L) = S
      IF (PRT) WRITE(6,2001) L,(XL(I,L),I=1,NDM)
20    CONTINUE
C---- COMPUTE NATURAL COORDINATE INCREMENTS
      DR = 2./NR
      DS = 2./NS
C---- COMPUTE LAST ELEMENT NUMBER AND CHECK ADMISSIBILITY
      NF = NE + (NR+NGAP)*NS - (1+NGAP)
      IF(NF.GT.NUMEL.AND.NE.GT.0) GO TO 401
C---- COMPUTE LAST NODE NUMBER AND CHECK ADMISSIBILITY
      NR = NR + 1
      NS = NS + 1
      IF(NDM.EQ.1) NS = 1
      NG = NI + (NR+NGAP)*NS - (1+NGAP)
      IF(NG.GT.NUMNP) GO TO 400
      N = NI
      ME = NE - 1
      S = -1.D0
C---- LOOP OVER GENERATED NODES
      DO 200 J = 1,NS
      R = -1.D0
      DO 100 I = 1,NR
C---- COMPUTE INTERIOR NODE COORDINATES USING ISOPARAMETRIC
C---- INTERPOLATION FUNCTIONS
      CALL SHAP2D(R,S,XL,SHP,XSJ,NDM,NM,IXL,.FALSE.)
      DO 55 L = 1,NDM
      X(L,N) = 0.D0
      DO 50 K = 1,9
      M = IXL(K)
      IF(M.GT.0) X(L,N) = X(L,N) + SHP(3,M)*XL(L,M)
50    CONTINUE
55    CONTINUE
C---- LOOP OVER ELEMENTS AND SET UP CONNECTIVITY ARRAY
      N = N + 1
      IF(I.EQ.NR)THEN
      N = N+NGAP
      ME= ME+NGAP
      GO TO 100
      END IF
      IF(NE.LE.0.OR.(J.EQ.NS.AND.J.NE.1)) GO TO 100
      ME = ME + 1
      IX(NEL1,ME) = MA
      IX(1,ME)    = N - 1
      IX(2,ME)    = N
      IF(NDM.EQ.1) GO TO 100
      IX(3,ME)    = N + NR + NGAP
      IX(4,ME)    = N + NR + NGAP - 1
100   R = R + DR
200   S = S + DS
      IF(.NOT.PRT) RETURN
C---- PRINT NODAL LIST
      K=0
      N = NI
      NK= NI+NR-1
      DO 501 I2 = 1,NS
      DO 502 I1 = N,NK
      K=K+1
      IF(K.EQ.50)THEN
      WRITE(6,2003) 0,HEAD,(I,XH,I=1,NDM)
      K=0
      END IF
  502 WRITE(6,2004) I1,(X(K,I1),K=1,NDM)
      N = NK+NGAP+1
  501 NK= N+NR-1
C---- PRINT ELEMENT LIST
      IF(NE.LE.0) RETURN
      K=0
      NR = NR-1
      NS = NS-1
      N = NE
      NK= NE+NR-1
      DO 503 I2 = 1,NS
      DO 504 I1 = N,NK
      K=K+1
      IF(K.EQ.50)THEN
      WRITE(6,2005) O,HEAD,(I,I=1,NEL)
      K=0
      END IF
  504 WRITE(6,2006) I1,MA,(IX(K,I1),K=1,NEL)
      N = NK+NGAP+1
  503 NK= N+NR-1
      RETURN
400   WRITE(6,2030) NG,NUMNP
      RETURN
401   WRITE(6,2031) NF,NUMEL
      RETURN
 1000 FORMAT(16I5)
 1001 FORMAT(I5,3F10.0)
2000  FORMAT(A1,20A4//17H NODE GENERATIONS//
     1   10X,25HNUMBER OF R-INCREMENTS    ,I5/
     2   10X,25HNUMBER OF S-INCREMENTS    ,I5/
     3   10X,25HFIRST NODE NUMBER         ,I5/
     4   10X,25HFIRST ELEMENT NUMBER      ,I5/
     5   10X,25HELEMENT MATERIAL NUMBER   ,I5/)
2001  FORMAT(I9,3E12.3)
2002  FORMAT(5X,4HNODE,3(I6,A6))
2003  FORMAT(/A1,20A4//5X,17HNODAL COORDINATES//6X,4HNODE,3(I7,A6))
2004  FORMAT(I10,3F13.4)
2005  FORMAT(/A1,20A4//5X,19HELEMENT CONNECTIONS//10H   ELEMENT,7X,4HMATL
     1   ,9(I5,5H NODE))
2006  FORMAT(11I10)
2010  FORMAT(5X,'WARNING * * * NO ELEMENTS ARE GENERATED ')
2030  FORMAT(5X,'**FATAL ERROR 31** INSUFFICIENT STORAGE FOR NODES'/
     1   10X,12HFINAL NODE =,I5,5X,7HNUMNP =,I5)
2031  FORMAT(5X,'**FATAL ERROR 32** INSUFFICIENT STORAGE FOR ELEMENTS'/
     1   10X,'FINAL ELEMENT =',I5,5X,'NUMEL =',I5)
      END
C
      SUBROUTINE DOMCOR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL PRTN,PRT
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,
     1 IPR,NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
      COMMON /PMPRTN/ PRT
      COMMON /JINT/ NDOM(3),XC(8,4),DOM(4,10),DYNJ(8,7,10),
     1 DYNJM(8,7,10),TSWE,TKEE,TKEE1
C     DIMENSION IE(7,1),D(18,1),ID(NDF,1),X(NDM,1),IX(NEN1,1),
C    1 IDL(1),F(NDF,1),T(1),DR(NDF,1),STR(NSTR,1),EPS(NSTR,1),
C    2 Q(NQ,1),VEL(NDF,1),ACCEL(NDF,1)
C
C---- READ CRACK-TIP COORDINATES AND DOMAIN PARAMETERS
C---- READ NUMBER OF DOMAINS IN X-Y PLANE, SYMMETRIC, NUMBER OF CRACK TIPS
C
      READ(5,1005) NDOM(1),NDOM(2),NDOM(3)
      IF (NDOM(1).GT.10) THEN
       WRITE(6,2005)
       STOP
      END IF
      IF (NDOM(3).GT.8) THEN
       WRITE(6,2010)
       STOP
      END IF
      IF (NDOM(3).EQ.0) THEN
       NDOM(3) = 1
       WRITE(6,2015)
      ELSE
       WRITE(6,2020) NDOM(2),NDOM(3)
      END IF
C---- READ CRACK TIP COORDINATES AND VCE AREA
      DO 20 NDI = 1,NDOM(3)
      READ(5,1010) (XC(NDI,I),I=1,4)
      IF (XC(NDI,4).EQ.0.D0) XC(NDI,4) = 1.D0
  20  WRITE(6,2025) (XC(NDI,I),I=1,4)
      WRITE(6,2030) NDOM(1),NDOM(2)
      DO 30 I = 1,NDOM(1)
      READ(5,1010) (DOM(J,I),J=1,4)
 30   WRITE(6,2035) I,(DOM(J,I),J=1,4)
      RETURN
C
 1005 FORMAT (3I5)
 1010 FORMAT (4F10.0)
 2005 FORMAT(/5X,' ** FATAL ERROR MAXIMUM NUMBER (10) OF DOMAINS ',
     1      ' EXCEEDED IN DOMA ** ')
 2010 FORMAT(/5X,' ** FATAL ERROR MAXIMUM NUMBER (8) OF CRACK TIPS ',
     1      ' EXCEEDED IN DOMA ** ')
 2015 FORMAT(/5X,' 2 - D   J-INTEGRAL ANALYSIS'
     1   //5X,' CRACK TIP COORDINATES'
     2   /8X,'  XC',9X,' YC',9X,' ZC',7X,'VCE AREA'/)
 2020 FORMAT(/5X,' 3 - D   J-INTEGRAL ANALYSIS'
     1   /10X,'NUMBER OF LAYERS',I4,'    NUMBER OF CRACK TIPS',I4
     2   //5X,' CRACK TIP COORDINATES'
     3   /8X,'  XC',9X,' YC',9X,' ZC',7X,'VCE AREA'/)
 2025 FORMAT (5X,4G12.3)
 2030 FORMAT(/5X,' NUMBER OF DOMAINS : ',I5,'   SYMMETRIC  Y=1,N=0',
     1       I5//5X,
     2 ' DOMAIN       DX (DR)     DY          DW     RECT(1) CIRC(2)',/)
 2035 FORMAT (5X,I5,5X,4G12.3)
      END
C
      SUBROUTINE DUMREAD(Q,NQ)
C
C---- SUBROUTINE FOR READING INTERNAL VARIABLES FROM A SEPERATE DATA FILE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CDATA/ NUMNP, NUMEL, NUMMAT,NEN,NEQ,IPR
     1 ,NSDM,NQDM,NQUAD
      DIMENSION Q(NQ,*),A(2,5)
      DO 10 I=1,NUMEL
      DO 11 J=1,2
      READ(49,1000) (A(J,K),K=1,5)
11    CONTINUE
      DO 12 J=1,5
      DO 12 L=1,2
      Q((J-1)*2+L,I) = A(L,J)
12    CONTINUE
10    CONTINUE
      RETURN
1000  FORMAT(5E13.5)
      END
