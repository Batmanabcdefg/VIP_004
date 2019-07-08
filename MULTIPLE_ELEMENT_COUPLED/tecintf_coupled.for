      SUBROUTINE PLOT(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,STR,EPS,Q,
     1 B,DR,VELG,ACCELG,VEL,ACCEL,NDF,NDM,NEN1,NST,NSTR,NQ,NEND,KEY,
     2 IQ,PRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- CONVERT DATA INTO TECPLOT FORMAT
C
      LOGICAL PRT,DUMPCF1,DUMPCF2
      LOGICAL AFR,BFR,CFR,AFL,BFL,CFL,DFL,EFL,GFL,TFL,
     1 IFL,BLKFL,SRCLIN,HFL
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /MAINN/ M(20000000)
      COMMON /POINTR/ NE,NN,NA,NC,NV,NM,NA2,NA3 
      COMMON /FLAGS/ AFR,BFR,CFR,AFL,BFL,CFL,DFL,EFL,GFL,TFL,IFL,BLKFL
     1 ,SRCLIN,HFL
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /TDATA/ TIME,DT,BETA,GAMMA
      COMMON /TECDAT/ DD(3)
      common /hadapt/ enorm,unorm,ismth
      DIMENSION JDIAG(*),UL(*),XL(*),TL(*),LD(*),P(*),S(*),IE(*),D(*),
     2 IX(NEN1,*),F(*),T(*),B(NDF,*),DR(*),ID(*),X(NDM,*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(*),ACCELG(*),VEL(*),ACCEL(*)
      dimension dev(4),dev1(4)
      umag = 1.d0
      inquire(61,opened=dumpcf1)
      if(.not.dumpcf1) then
      open(unit=61,file='tecplt1.dat',form='formatted')                         !tecplt1.dat is unit 61
      end if
      if(iq.ne.0) then
      inquire(62,opened=dumpcf2)
      if(.not.dumpcf2) then
      open(unit=62,file='tecplt2.dat',form='formatted')                         !tecplt2.dat is unit 62
      end if
      end if
c---- form stiffness matrix for least-squares problem
      IF(GFL) CALL PSETM(NA,NE,JDIAG(NEQ)*IPR,GFL)
      CALL PZERO(M(NA),JDIAG(NEQ))
      CALL PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
     2 STR,EPS,Q,DR,M(NA),M(NC),
     3 VELG,ACCELG,VEL,ACCEL,NDF,NDM,
     4 NEN1,NST,NSTR,NQ,3,B,M(NV),.TRUE.,.FALSE.,.FALSE.,
     5 .FALSE.)
c---- smooth stresses  
      DD(1) = 1.d0
      DD(2) = 2.d0
      ismth = 1
c----------------------------------------------------------
c---- foll changes introduced to plot slip resistance on primary / conjugate slip systems 
c     if(iq.gt.200) then
c     iprslip = iq - 200
c     if(iprslip.le.2) then
c     isecr = iprslip + 2
c     else 
c     isecr = iprslip -2
c     end if
c     DD(1) =104.d0  + dfloat(iprslip)
c     DD(2) = 104.d0 + dfloat(isecr)
c     end if
c     if(iq.gt.100) then
c       dd(1) = dfloat(iq)
c       dd(2) = dfloat(iq+2)
c     end if
c-----------------------------------------------------------
      CALL PZERO(DR,NEQ)
      CALL PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
     1   STR,EPS,Q,DR,DR,DR,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,6,B,
     3   M(NV),.FALSE.,.TRUE.,.FALSE.,.FALSE.)
      CALL ACTCOL(M(NA),DR,JDIAG,NEQ,.TRUE.,.TRUE.)
      DD(1) = 3.d0
      DD(2) = 4.d0
      ismth = 2
c----------------------------------------------------------------
c---- foll change introduced to plot slip on primary / conjugate slip systems 
c     if(iq.gt.200) then
c     iprslip = iq - 200
c     if(iprslip.le.2) then
c     isecr = iprslip + 2
c     else 
c     isecr = iprslip -2
c     end if
c     DD(1) =108.d0  + dfloat(iprslip)
c     DD(2) = 108.d0 + dfloat(isecr)
c      if(iq.gt.100) then
c       dd(1) = dfloat(100+iq)
c       dd(2) = dfloat(100+iq+5)
c     end if     
c--------------------------------------------------------------
      CALL PZERO(VELG,NEQ)
      CALL PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
     1   STR,EPS,Q,VELG,VELG,VELG,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,6,B,
     3   M(NV),.FALSE.,.TRUE.,.FALSE.,.FALSE.)
      CALL ACTCOL(M(NA),VELG,JDIAG,NEQ,.FALSE.,.TRUE.)
c---- check smoothed stresses at center gauss points
c     CALL PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
c    1   STR,EPS,Q,VELG,VELG,VELG,VELG,ACCELG,VEL,ACCEL,
c    2   NDF,NDM,NEN1,NST,NSTR,NQ,4,B,
c    3   M(NV),.FALSE.,.FALSE.,.FALSE.,.FALSE.)
c---- write to tecplot file
      xkey1 = 4hstre
      write(61,'(''TITLE = " '',a4,'' "'')') xkey1
      do 15 i=1,numnp
c        dev(1) = dr(2*i)
c        dev(2) = velg(2*(i-1)+1)
c        dev(3) = velg(2*i)
c        dev(4) = dr(2*(i-1)+1)
         dev(1) = dr(4*i-3)
         dev(2) = dr(4*i-2)
         dev(3) = velg(4*i-3)
         dev(4) = velg(4*i-2)
c        dev(4) = dr(2*(i-1)+1)
c---- store hydrostatic stresses in dev(4)
         dev(3) = (dev(1)+dev(2)+dev(3))/3.d0
c     call pstress(dev,dev1)
      if (i.eq.1) then
      write(61,'('' VARIABLES = X Y S11 S22 SHYD TAU'')')
      write(61,'('' ZONE T = "ZONE ONE", I = '',i5,'', J = '',i5,
     1 '', F = FEPOINT'')') numnp,numel
      end if
      xdef = x(1,i) + umag*b(1,i)
      ydef = x(2,i) + umag*b(2,i)
      write(61,'(6e13.5)') xdef,ydef,(dev(k),k=1,4)
15    continue
c---- write i.vars
      if(iq.ne.0) then
c     DD(1) = 4.d0 + dfloat(iq)
c     DD(2) = 4.d0 + 2.d0
c     DD(2) = 4.d0 + dfloat(iq+1) 
c     dd(1) = dfloat(iq)
c     dd(2) = dfloat(iq+1)
      dd(1) = dfloat(iq+6)                                 !temperature t_n+1
      dd(2) = dfloat(iq+8)                                 !free volume xi_n+1
      CALL PZERO(DR,NEQ)
      CALL PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
     1   STR,EPS,Q,DR,DR,DR,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,6,B,
     3   M(NV),.FALSE.,.TRUE.,.FALSE.,.FALSE.)
      CALL ACTCOL(M(NA),DR,JDIAG,NEQ,.FALSE.,.TRUE.)
      dd(1) = dfloat(iq+104)                               !log(lambda1_p)
      dd(2) = dfloat(iq+122)                               !cohesion t_n+1
      CALL PZERO(VELG,NEQ)
      CALL PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
     1   STR,EPS,Q,VELG,VELG,VELG,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,6,B,
     3   M(NV),.FALSE.,.TRUE.,.FALSE.,.FALSE.)
      CALL ACTCOL(M(NA),VELG,JDIAG,NEQ,.FALSE.,.TRUE.)
c---- check smoothed internal variables at center gauss points
c     CALL PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,F,T,JDIAG,
c    1   STR,EPS,Q,VELG,VELG,VELG,VELG,ACCELG,VEL,ACCEL,
c    2   NDF,NDM,NEN1,NST,NSTR,NQ,4,B,
c    3   M(NV),.FALSE.,.FALSE.,.FALSE.,.FALSE.)
c
      xkey2 = 4hivar
      write(62,'(''TITLE = " '',a4,'' "'')') xkey2
      do 20 i=1,numnp
         dev(1) = dr(4*i-3)
         dev(2) = dr(4*i-2)
         dev(3) = velg(4*i-3)
         dev(4) = velg(4*i-2)
      if (i.eq.1) then
      write(62,'('' VARIABLES = X Y Q2 Q3 Q4 Q5'')')
      write(62,'('' ZONE T = "ZONE ONE", I = '',i5,'', J = '',i5,
     1 '', F = FEPOINT'')') numnp,numel
      end if
      xdef = x(1,i) + umag*b(1,i)
      ydef = x(2,i) + umag*b(2,i)
      write(62,'(6e13.5)') xdef,ydef,(dev(k),k=1,4)
20    continue
      end if
c---- write element connectivities
      do 21 n = 1,numel
      if(nen.eq.4) then
      ix4 = ix(4,n)
      if(ix(4,n).eq.0) ix4 = ix(3,n)
      if(ix4.ne.0) then
      write(61,'(4i5)') ix(1,n),ix(2,n),ix(3,n),ix4
      if(iq.ne.0)   
     1 write(62,'(4i5)') ix(1,n),ix(2,n),ix(3,n),ix4
      end if
      elseif(nen.eq.6) then
      write(61,'(4i5)') ix(1,n),ix(4,n),ix(6,n),ix(1,n)
      write(61,'(4i5)') ix(4,n),ix(2,n),ix(5,n),ix(4,n)
      write(61,'(4i5)') ix(5,n),ix(3,n),ix(6,n),ix(5,n)
      write(61,'(4i5)') ix(4,n),ix(5,n),ix(6,n),ix(4,n)
      if(iq.ne.0) then
      write(62,'(4i5)') ix(1,n),ix(4,n),ix(6,n),ix(1,n)
      write(62,'(4i5)') ix(4,n),ix(2,n),ix(5,n),ix(4,n)
      write(62,'(4i5)') ix(5,n),ix(3,n),ix(6,n),ix(5,n)
      write(62,'(4i5)') ix(4,n),ix(5,n),ix(6,n),ix(4,n)
      end if
      end if
   21 continue
      return
      end
      subroutine pstress(dev,dev1)
      implicit double precision (a-h,o-z)
c----  calculate principal stresses or strains
      dimension dev(4),dev1(4)
      devm = 0.5*(dev(1)+dev(2))
      devsq=dsqrt((dev(1)-dev(2))*(dev(1)-dev(2))*0.25 + dev(4)*dev(4))
      dev1(1) = devm + devsq
      dev1(2) = devm - devsq
      dev1(3) = dev(3)
      dev1(4) = 0.5*(dev1(1) - dev1(2))
      return
      end
c---------------------------------------------------------------------
      SUBROUTINE PFORMSM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,
     1   F,T,JDIAG,STR,EPS,Q,B,A,C,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,ISW,U,UD,AFL,BFL,CFL,DFL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C---- COMPUTE ELEMENT ARRAYS AND ASSEMBLE GLOBAL ARRAYS
C---- ASSEMBLING STIFFNESS MATRIX AND FORCE VECTORS FOR LEAST    
C---- SQUARES SMOOTHING OF STRESSES / STRAINS / I.VARS.
C---- NOTE: SMOOTH ELMT SUBROUTINE SHOULD BE ELMT03
C---- NO BCS ARE ENFORCED WHILE GLOBALLY SMOOTHING
C
      LOGICAL AFL,BFL,CFL,DFL,FLAG
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,N,MA,MCT,IEL,NEL
      COMMON /TECDAT/ DD(3)
      DIMENSION XL(NDM,*),LD(NDF,*),P(*),S(NST,*),IE(7,*),D(40,*),
     1 ID(NDF,*),X(NDM,*),IX(NEN1,*),F(NDF,*),JDIAG(*),B(*),A(*),C(*),
     2 UL(NDF,*),TL(*),T(*),U(NDF,*),UD(NDF,*)
      DIMENSION STR(NSTR,*),EPS(NSTR,*),Q(NQ,*)
      DIMENSION VELG(NDF,*),ACCELG(NDF,*),VEL(NDF,*),ACCEL(NDF,*)
C---- LOOP ON ELEMENTS
      IEL = 3
      DO 110 N = 1,NUMEL
C---- SET UP LOCAL ARRAYS
      MA = IX(NEN1,N)
C---- SKIP GAP ELEMENTS
      IF(IE(7,MA).EQ.2) GO TO 110
C
      DO 108 I = 1,NEN
      II = IX(I,N)
      IF(II.NE.0) GO TO 105
      TL(I) = 0.
      DO 103 J = 1,NDM
103   XL(J,I) = 0.
      DO 104 J = 1,NDF
      VEL(J,I) = 0.d0
      ACCEL(J,I) = 0.d0
      UL(J,I) = 0.d0
      UL(J,I+NEN)=0.d0
104   LD(J,I) = 0
      GO TO 108
105   CONTINUE
      NEL = I
      TL(I) = T(II)
      DO 106 J = 1,NDM
106   XL(J,I) = X(J,II)
      DO 107 J = 1,NDF
      VEL(J,I) = VELG(J,II)
      ACCEL(J,I) = ACCELG(J,II)
      UL(J,I) = U(J,II)
      UL(J,I+NEN)=UD(J,II)
      LD(J,I) = (II-1)*NDF + J
107   CONTINUE
108   CONTINUE
C---- FORM ELEMENT ARRAY
c---- put in mode switch
      dd(3) = d(3,ma)
      CALL ELMLIB(DD,UL,XL,STR(1,N),EPS(1,N),Q(1,N),IX(1,N),TL,
     1 S,P,VEL,ACCEL,NDF,NDM,NST,NSDM,NQDM,ISW)
C---- ASSEMBLE INTO GLOBAL ARRAY
      CALL ADDSTF(A,B,C,S,P,JDIAG,LD,NST,NEL*NDF,AFL,BFL,CFL)
110   CONTINUE
c---- put 1s on diagonal of both dofs of last 2 nodes
c     if(afl) then
c     do 200 k=1,2
c     kdf1 = (numnp-2)*ndf + k
c     kdf2 = (numnp-1)*ndf + k
c     a(jdiag(kdf1)) = 1.d0
c     a(jdiag(kdf2)) = 1.d0
      if(afl) then
      do 200 k=1,neq
      if(a(jdiag(k)).eq.0.d0) a(jdiag(k)) = 1.d0
200   continue
      end if
      RETURN
      END
   
