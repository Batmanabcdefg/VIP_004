      subroutine ma01(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f,t,jdiag,
     1 str,eps,q,b,dr,velg,accelg,vel,accel,ct,
     2 ndf,ndm,nen1,nst,nstr,nq,ne1,flag)
      implicit double precision (a-h,o-z)
c
c---- macro to solve coupled thermo-mechanical problems using
c---- explicit time integration method.
c---- this version assumes that ndf = ndm +1, where the first ndm
c---- dofs are the translational displacement dofs and the last one
c---- is temperature.
c
      logical flag
      LOGICAL AFR,BFR,CFR,AFL,BFL,CFL,DFL,EFL,GFL,TFL,
     1 IFL,BLKFL,SRCLIN,HFL
      COMMON /MAINN/ M(20000000)
      COMMON /POINTR/ NE,NN,NA,NC,NV,NM,NA2,NA3
      COMMON /FLAGS/ AFR,BFR,CFR,AFL,BFL,CFL,DFL,EFL,GFL,TFL,IFL,BLKFL
     1 ,SRCLIN,HFL
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /TDATA/ TIME,DT,BETA,GAMMA
      COMMON /MDATA/ NN1,N0,N1,N2,N3,N301,N302,N4,N5,N6,N7,N8,N9,N91,
     1 N10,N11,N12,N13,N14,N15,N16,N161,N162,N163,N164,N165,N166,
     1 N17,N18,N181,N182,N183,N184,N185,N186,N19,N20
      COMMON /PTABLE/ NPLD,NDIST,T1,P1,T2,P2
c
      dimension ul(*),xl(*),tl(*),ld(*),p(*),s(*),ie(*),d(*),
     1 id(ndf,*),x(ndm,*),ix(nen1,*),f(ndf,*),t(*),jdiag(*),
     2 str(*),eps(*),q(*),b(ndf,*),dr(*),velg(ndf,*),
     3 accelg(ndf,*),vel(*),accel(*),ct(*)
      if (efl) then
      write(6,'(///'' **error** detected by macro excdh'',
     1 '' mass matrix not computed'')')
      stop
      end if
      IF(.NOT.TFL) GAMMA = 0.5
      call excdh(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f,t,jdiag,
     1 str,eps,q,b,dr,velg,accelg,vel,accel,ct,m(n162),
     2 m(nn),m(n91),m(n165),m(n166),m(n301),m(n302),ndf,ndm,nen1,nst,
     3 nstr,nq,npld,ndist)
      return
      end
c
      subroutine excdh(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f,t,jdiag,
     1 str,eps,q,b,dr,velg,accelg,vel,accel,ct,unew,
     2 xmass,idprop,ttable,ptable,pdist,idist,ndf,ndm,nen1,nst,
     3 nstr,nq,npld,ndist)
      implicit double precision (a-h,o-z)
c
c---- assume ndf = ndm +1, where the last dof corresponds to temperature
c
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,
     1 NSDM,NQDM,NQUAD
      COMMON /TDATA/ TIME,DT,BETA,GAMMA
      COMMON /PRLOD/ PROP(7),PROPOL(7)
c
      dimension ul(*),xl(*),tl(*),ld(*),p(*),s(*),ie(*),d(*),
     1 id(ndf,*),x(ndm,*),ix(nen1,*),f(ndf,*),t(*),jdiag(*),
     2 str(nstr,*),eps(nstr,*),q(nq,*),b(ndf,*),dr(ndf,*),velg(ndf,*),
     3 accelg(ndf,*),vel(*),accel(*),ct(*),xmass(ndf,*),idprop(ndf,*),
     4 ptable(7,*),ttable(*),unew(ndf,*),pdist(nst,*),
     5 idist(2,*)
C---- LOAD/DISP TERMS
c     write(6,9001) ((xmass(j,i),j=1,ndf),i=1,numnp)
9001  format(6e13.5)
      IF(NPLD.GT.0) CALL PROPLD(TTABLE,PTABLE,PROP,TIME,NPLD)
C     CALL PLOAD(ID,IX,PDIST,IDIST,F,DR,NEQ,IDPROP,PROP,NEN,NEN1,NDF,
C    1 NST,NDIST)
C---- UPDATE DISPLACEMENTS
      DTH = DT*0.5
      DT2 = DT*DTH
      DO 10 I=1,NUMNP
      DO  5 J = 1,NDF
      B1 = B(J,I)
C      IF(J.EQ.NDF) THEN
      IF((NDF .GT.NDM) .AND. (J .GT. NDM)) THEN
C----- TEMPERATURE DOF -----------
      ACCELG(J,I) = 0.D0
c     T(I) = B1
      END IF
C
      IF(ID(J,I).EQ.0) THEN
      B(J,I) = B(J,I) + DT*VELG(J,I)+DT2*ACCELG(J,I)
      ELSE
      IND = IDPROP(J,I)
      B(J,I) = F(J,I)*PROP(IND)
      END IF
      UNEW(J,I) = B(J,I) - B1
5     CONTINUE
10    CONTINUE
C---- INTEGRATE CONSTITUTIVE EQNS
c---- First get material derivatives of free volume at element centroid
      CALL PFORM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,JDIAG,
     1 STR,EPS,Q,DR,DR,DR,VELG,ACCELG,VEL,ACCEL,
     2 NDF,NDM,NEN1,NST,NSTR,NQ,2,B,
     3 UNEW,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
c---- Get nodal values of free volume gradients by nodal averaging
       call nodal_fvol_grads(q,ix,nen1,nen,nq,nqdm,nquad,
     1 numnp,numel)
c
      CALL PFORM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,JDIAG,
     1 STR,EPS,Q,DR,DR,DR,VELG,ACCELG,VEL,ACCEL,
     2 NDF,NDM,NEN1,NST,NSTR,NQ,9,B,
     3 UNEW,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
C---- FORM RHS
      CALL PLOAD(ID,IX,PDIST,IDIST,F,DR,NEQ,IDPROP,PROP,NEN,NEN1,NDF,
     1 NST,NDIST)
C     CALL PFORM(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,JDIAG,
      CALL PFORM2(UL,XL,TL,LD,P,S,IE,D,ID,X,IX,IDPROP,PROP,F,T,
     1   JDIAG,STR,EPS,Q,DR,DR,DR,VELG,ACCELG,VEL,ACCEL,
     2   NDF,NDM,NEN1,NST,NSTR,NQ,6,B,
     3   UNEW,.FALSE.,.TRUE.,.FALSE.,.FALSE.)
c     write(6,9001) ((dr(j,i),j=1,ndf),i=1,numnp)
C---- SOLVE FOR ACCELERATIONS; UPDATE VELOCITIES (MECHANICAL DOFS)
      DO 20 I = 1,NUMNP
       IF (NDF .NE. 1) THEN 
      DO 25 J=1,NDM
      ACC = ACCELG(J,I)
      IF((ID(J,I).EQ.0).AND.(XMASS(J,I).GT.0.D0)) THEN
      ACCELG(J,I) = DR(J,I)/XMASS(J,I)
      VELG(J,I) = VELG(J,I) + DT*(GAMMA*ACCELG(J,I)+(1.D0-GAMMA)*ACC)
C     VELG(J,I)   = VELG(J,I)+DTH*(ACCELG(J,I)+ACC)
      ELSE
      VELPRV = VELG(J,I)
      VELG(J,I) = UNEW(J,I)/DT
      ACCELG(J,I) = (VELG(J,I)-VELPRV)/DT
      END IF
25    CONTINUE
       END IF
C---- SOLVE FOR TEMPERATURE INCREMENTS (THERMAL DOF)
       IF (NDF .GT.NDM) THEN
      do j=ndm+1,ndf
c     IF(ID(j,I).EQ.0) THEN
      IF((ID(J,I).EQ.0).AND.(XMASS(J,I).GT.0.D0)) THEN
      VELG(J,I) = DR(J,I)/XMASS(J,I)
      ELSE
      VELG(J,I) = UNEW(J,I)/DT
      END IF
      ACCELG(J,I) = 0.0
      end do
       END IF
20    CONTINUE
      RETURN
      END
c
      subroutine nodal_fvol_grads(q,ix,nen1,nen,nq,nqdm,nquad,
     1 numnp,numel)
      implicit double precision (a-h,o-z)
      common /freevolgrad/ xigrad(2,5000)
      dimension q(nq,*),ix(nen1,*),near_elems(5000)
c
      if((numel.gt.5000).or.(numnp.gt.5000)) then
       write(6,2001) numnp,numel
       stop
      end if 
c---- Initialize nodal values of material gradients of free volume
      do i = 1,numnp
       xigrad(1,i) = 0.d0
       xigrad(2,i) = 0.d0
       near_elems(i) = 0
      end do
c---- Determine nodal values of free volume gradients by nodal
c---- averaging based on nearby element centroidal values
      do  n=1,numel
      do  j=1,nen
        ii = ix(j,n)
        xigrad(1,ii) = xigrad(1,ii) + q((nquad-1)*nqdm + 23,n)
        xigrad(2,ii) = xigrad(2,ii) + q((nquad-1)*nqdm + 24,n)
        near_elems(ii) = near_elems(ii) + 1
      end do
      end do
      do i = 1,numnp
        if(near_elems(i).eq.0) then
          write(6,2002) i
          stop
        end if
       xigrad(1,i) = xigrad(1,i)/dfloat(near_elems(i))
       xigrad(2,i) = xigrad(2,i)/dfloat(near_elems(i))
      end do
      return
2001  format(5x,'** Fatal Error in EXCDH. #nodes or elems exceeds 5000',
     1 2i6)
2002  format(5x,'** Fatal Error in EXCDH. #nearby elems for node #',
     1 i6,' is zero')
      end
