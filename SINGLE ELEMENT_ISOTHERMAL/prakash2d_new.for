      subroutine elmt01(d,ul,xl,strl,epsl,ql,ix,tl,s,p,
     1 vel,accel,ndf,ndm,nst,ns,nq,isw)
      implicit double precision (a-h,o-z)
c
c---- 8-node, finite plasticity element
c---- Fe-Fp decomposition. Explicit update
c---- Assumed strain formulation for near-incompressibility.
c---- Coupled mechanics-free volume evolution analysis: NDF = 3
c
c---- strl :  Cauchy stress tensor 
c---- epsl :  2nd P-K stress based on intermediate configuration
c---- ql   :  internal variables on intermediate configuration
c---- ul   :  total displacements
c
c---- Both initial and updated state variables are stored
c---- Foll model specific subs are to be coded separately :
c     ivarout,dtang, update, elworkdens, totworkdens, pkstrs1
c---- Also needed from library file are foll subs:
c     shap2d & djfn
c
      logical prt,option
      CHARACTER*8 HEAD
      CHARACTER*1 O 
      common /tdata/ time,dtime,bet1,gam1
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,
     1 NEN,NEQ,IPR,NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL
      COMMON /JINT/ NDOM(3),XC(8,4),DOM(4,10),DYNJ(8,7,10)
     1 ,DYNJM(8,7,10),TSWE,TKEE,tkee1
      COMMON /HYDSMAC/ s11_mac,s22_mac,shyd_mac,tot_vol,ncount
      COMMON /CNVG1/ ICONV,IDIVG,IPATH
      common /print/ prt
      COMMON /PRLOD/ PROP(7),PROPOL(7)
      dimension d(*),ul(ndf,*),xl(ndm,*),ix(*),tl(*),s(nst,*),
     1          p(*),shp(3,9),sg(16),tg(16),wg(16),xlc(2,9)
      dimension strl(ns,*),epsl(ns,*),ql(nq,*),vel(ndf,*),accel(ndf,*)
      dimension b(4,18),x(2),f(3,3),finv(3,3),bvol(18),shpb(3,9)
      dimension sig(6),fp(9),fpnew(9),snew(9),qnew(100)
      dimension squad1(8,8,4),pquad1(8,4),pk1(3,3),pkst(5),
     1 squad2(4,4,4),pquad2(4,4)
      DIMENSION QW(9),DJ(7,10),DJM(7,10),p1(8),p2(4),s1(8,8),s2(4,4)
      DATA PI/3.141592654/
c
      ndfm = 2
      nstm = nen*ndfm
c---- go to correct array processor
      go to (1,2,3,4,5,3,2,2,9,10,11), isw
c
c---- input material properties
  1   continue
c
c---- d(1)      x-gravity acceleration
c---- d(2)      y-gravity acceleration
c---- d(3)      mode : 1 (plane strain); 2(axisymm)
c---- d(4)      stabilization parameter
c---- d(5)      mass density    
c
      read(5,1001)  (d(i),i=1,2), mode, (d(i),i=4,5)
      d(3) = mode
      write(6,2000) d(1),d(2),mode,d(4),d(5)
c---- input model dependent parameters
      call input(d(6))
      return
c
c----  calculate gradients of free volume at element centroid
  2   continue
      ql(23,5) = 0.d0
      ql(24,5) = 0.d0
      CALL SHAP2D(0.D0,0.D0,XL,SHPB,XSJ,NDM,NEL,IX,.TRUE.)
      do i = 1,nel
       ql(23,5) = ql(23,5) + shpb(1,i)*ul(3,i)
       ql(24,5) = ql(24,5) + shpb(2,i)*ul(3,i)
      end do
       return
c
c
  3   continue
c---- compute tangent stiffness matrix and out of balance forces
c---- stabilization factors 
      aa = d(4)
      bb = 1.0d0 - d(4)
      mode = d(3)
c---- compute current nodal coordinates
      do  i = 1,ndm
       do j= 1, nen
         xlc(i,j) = xl(i,j) + ul(i,j)
       end do
      end do
c---- retrieve shape functions at centroid
      CALL SHAP2D(0.D0,0.D0,XL,SHPB,XSJ,NDM,NEL,IX,.TRUE.)
c---- compute jacobian of deformation at centroid
      call defg(xl,xlc,ul,shpb,f,finv,theta,ndm,ndf,nel,mode,.true.)
c---- compute spatial shape functions
      do 353 i = 1,nel
         x(1) = shpb(1,i)*finv(1,1) + shpb(2,i)*finv(2,1) 
         x(2) = shpb(1,i)*finv(1,2) + shpb(2,i)*finv(2,2) 
         shpb(1,i) = x(1)
         shpb(2,i) = x(2)
353   continue
c---- compute volumetric compatibility matrix
      CALL BMAT(B,XJC,SHPB,XLC,NDM,NDF,NEL,NSTM,MODE)
      do 355 i = 1,nstm
         bvol(i) = b(1,i) + b(2,i) + b(3,i)
355   continue
c
c---- Main integration loop over regular quadrature points
      if ((isw.eq.3).or.(isw.eq.7)) then
         do  l = 1,4
            do  j = 1,8
               do  i = 1,8
                  squad1(i,j,l) = 0.d0
               end do
            end do
            do j=1,4
              do i=1,4
               squad2(i,j,l) = 0.d0
              end do
            end do
         end do
       do i = 1,8
        do j = 1,8
         s1(i,j) = 0.d0
        end do
       end do
       do i=1,4
        do j=1,4
         s2(i,j) = 0.d0
        end do
       end do
      end if
      if ((isw.eq.6).or.(isw.eq.7)) then
         do  l = 1,4
            do  i = 1,8
               pquad1(i,l) = 0.d0
            end do
            do i=1,4
              pquad2(i,l) = 0.d0
            end do
         end do
      do i=1,8
        p1(i) = 0.d0
       end do
      do i = 1,4
       p2(i) = 0.d0
      end do
      end if
      CALL PGAUS2(2,LINT,SG,TG,WG)
      do 300 l = 1,lint
c----  retrieve shape functions
      CALL SHAP2D(SG(L),TG(L),XL,SHP,XSJ,NDM,NEL,IX,.TRUE.)
         call sloop(d,ul,xl,xlc,strl(1,l),epsl(1,l),ql(1,l),ix,tl,
     1              squad1(1,1,l),pquad1(1,l),pquad2(1,l),
     2              vel,accel,theta,shpb,bvol,aa,bb,shp,
     3              xsj,ndf,ndfm,ndm,nst,nstm,nel,ns,nq,mode,isw)
300   continue
      if ((isw.eq.3).or.(isw.eq.7)) then
         do 395 l = 1,4
            do 395 j = 1,8
               do 395 i = 1,8
                  s1(i,j) = s1(i,j) + squad1(i,j,l)*wg(l)
395      continue
c------------------------------------------
c     write(6,9001) ((s(i,j),i=1,8),j=1,8)
9001  format(8e13.5)
c----------------------------------------------
      CALL STFSTORE(S,S1,S2,P,P1,P2,NDFM,NDF,NEL,NST,1)
      end if
      if ((isw.eq.6).or.(isw.eq.7)) then
         do  l = 1,4
            do  i = 1,8
               p1(i) = p1(i) + pquad1(i,l)*wg(l)
            end do
            do i=1,4
              p2(i) = p2(i) + pquad2(i,l)*wg(l)
            end do
        end do
        CALL STFSTORE(S,S1,S2,P,P1,P2,NDFM,NDF,NEL,NST,2)
c     write(6,9007) (p(i),i=1,12)
9007  format(2x,'elem force',6e13.5)
      end if
      return
c
C---- OUTPUT ELEMENT VARIABLES
    4 CONTINUE
C---- PRINT HEADER EVERY 25 ELEMENTS
      mode = d(3)
c---- compute current nodal coordinates
      do  i = 1,ndm
       do j= 1, nen
         xlc(i,j) = xl(i,j) + ul(i,j)
       end do
      end do
      MCT = MCT - 1
      IF (MCT.LE.0) THEN
      WRITE(6,2001) O,HEAD
      WRITE(6,2010)
      MCT = 25
      END IF
C---- SET UP OUTPUT POINTS AND SHAPE FUNCTIONS
      NINT = 1
C     NINT = 2
      IF (NEL.GT.4) NINT = 2
      CALL PGAUS2(NINT,LINT,SG,TG,WG)
      LL = 4
C     LL = 0
      IF(NEL.GT.4) LL = 9
      DO 400 L = 1,LINT
      LL = LL + 1
      CALL SHAP2D(SG(L),TG(L),XL,SHP,XSJ,NDM,NEL,IX,.false.)
c---- compute jacobian of deformation 
c     call defg(xl,xlc,ul,shp,f,finv,theta,ndm,ndf,nel,mode,.false.)
C---- COMPUTE OUTPUT POINT COORDINATES
      DO 410 I = 1,NDM
      X(I) = 0.D0
      DO 410 J = 1,NEL
  410 X(I) = X(I) + SHP(3,J)*XLC(I,J)
C---- PRINT VARIABLES
      WRITE(6,2002) NELMT,MA,(X(I),I=1,NDM),(STRL(I,LL),I=1,4),
     1 QL(3,LL),QL(4,LL),QL(5,LL),QL(7,LL)
c---- model dependent output of ivars (model specific subroutine) 
c      call ivarout(d(6),f,ql(1,ll),X,NDM)
  400 CONTINUE
      RETURN
c
  5   continue
      mode = d(3)
        do i=1,8
        p1(i) = 0.d0
       end do
      do i = 1,4
       p2(i) = 0.d0
      end do
C---- SET UP INTEGRATION POINTS
      CALL PGAUS2(2,LINT,SG,TG,WG)
C---- LOOP OVER INTEGRATION POINTS
      DO 500 L = 1,LINT
C---- COMPUTE SHAPE FUNCTIONS
      CALL SHAP2D(SG(L),TG(L),XL,SHP,XSJ,NDM,NEL,IX,.TRUE.)
      IF (MODE.EQ.2) THEN
      RADIUS = 0.D0
      DO 507 I = 1,NEL
  507 RADIUS = RADIUS + XL(1,I)*SHP(3,I)
      XSJ = 2.d0*PI*RADIUS*XSJ
      END IF
      TERM1 = WG(L)*XSJ*D(5)
      TERM2 = WG(l)*XSJ
      II = 1
      DO 505 I = 1,NEL
      W1 = SHP(3,I)*TERM1
      W2 = SHP(3,I)*TERM2
C---- LUMPED MASS APPROXIMATION
      P1(II) = P1(II) + W1
      P2(I) = P2(I) + W2
C---- CONSISTENT MASS (UPPER TRIANGULAR PART)
c     JJ = II
c     DO 510 J = I,NEL
c     S(II,JJ) = S(II,JJ) + SHP(4,J)*W
c 510 JJ = JJ + NDF
  505 II = II + NDFM
  500 CONTINUE
C---- COMPUTE MISSING PARTS BY SYMMETRY
      NSL = NEL*NDFM
      DO 515 I = 1,NSL,NDFM
      P1(I+1) = P1(I)
c     DO 515 J = I,NSL,NDF
c     S(I+1,J+1) = S(I,J)
c     S(J,I) = S(I,J)
c 515 S(J+1,I+1) = S(I,J)
515   continue
C---- STORE MECH & FREE VOLUME MASS MXS IN A SINGLE ARRAY P(NST)
      CALL STFSTORE(S,S1,S2,P,P1,P2,NDFM,NDF,NEL,NST,2)
c     write(6,9002) (p(i),i=1,12)
9002  format(6e13.5)
      return
c
  8   continue
      return
c
c---- update element variables at integration and output points
  9   continue
      aa = d(4)
      bb = 1.0d0 - d(4)
      mode = d(3)
c---- Initialize macroscopic stresses & body volume
      if(nelmt.eq.1) then
c       if(iconv.eq.1) then
             ncount = ncount + 1
             s11_mac = 0.d0
             s22_mac = 0.d0
             shyd_mac = 0.d0
             tot_vol = 0.d0
c       end if
      end if
c---- compute current nodal coordinates
      do  i = 1,ndm
       do j= 1, nen
         xlc(i,j) = xl(i,j) + ul(i,j)
       end do
      end do
c---- retrieve shape functions at centroid
      CALL SHAP2D(0.D0,0.D0,XL,SHPB,XSJ,NDM,NEL,IX,.TRUE.)
c---- compute deformation jacobian at centroid
      call defg(xl,xlc,ul,shpb,f,finv,theta,ndm,ndf,nel,mode,.false.)
c---- loop over integration points
      ll = 0
      do 900 kkk=1,2
       if(kkk.eq.1) nqpts = 2
       if(kkk.eq.2) nqpts = 1
       CALL PGAUS2(nqpts,LINT,SG,TG,WG)
      do 910 l = 1,lint
       ll = ll+1
c----  compute shape functions & derivatives
      CALL SHAP2D(SG(L),TG(L),XL,SHP,XSJ,NDM,NEL,IX,.TRUE.)
         call uloop(d(6),ul,xl,xlc,strl(1,ll),epsl(1,ll),ql(1,ll),ix,tl,
     & vel,accel,theta,aa,bb,shp,xsj,xsj1,ndf,ndm,nst,nel,ns,nq,mode,ll)
      if(kkk.eq.1) then
        xjc = 1.d0
        if(mode.eq.2) then
          xjc = 0.d0
           do i=1,nel
            xjc = xjc + shp(3,i)*xlc(1,i)
           end do
           xjc = xjc*2.d0*pi
         end if
          xsj1 = xsj1*xjc*wg(l)
      s11_mac = s11_mac + strl(1,ll)*xsj1
      s22_mac = s22_mac + strl(2,ll)*xsj1
      shyd_mac = shyd_mac +(strl(1,ll)+strl(2,ll)+strl(3,ll))*xsj1/3.d0
      tot_vol = tot_vol + xsj1
      end if
910   continue
900   continue
c---- Output, macroscopic stresses and total body volume to file fpseff.dat
c     if((nelmt.eq.numel).and.(iconv.eq.1)) then
      if((nelmt.eq.numel).and.(ncount.eq.10000)) then
         ncount = 0
         WRITE(22,'(2X,6E13.5)') TIME,PROP(1),s11_mac,s22_mac,shyd_mac,
     1   tot_vol
      end if
      return
c
 10   continue
      return
c
11    CONTINUE
C
      RETURN
c
c---- formats
 1000 format(8f10.0)
 1001 format(2f10.0,i10,2f10.0)
 2000 format(/5x,' FD MULTIPLICATIVE VISCOPLASTICITY ELEMENT'
     1 /5x,'MIXED FORMULATION FOR NEAR-INCOMPRESSIBILITY'
     2 /10x,'X-GRAVITY ACCELERATION        'e15.5
     3 /10x,'Y-GRAVITY ACCELERATION        'e15.5
     3 /10x,'MODE SWITCH                   'i10
     5 /10x,'STABILIZATION PARAMETER       'e15.5
     4 /10x,'MASS DENSITY                  'e15.5)
 2001 FORMAT(A1,20A4//5X,'ELEMENT VARIABLES'//
     1 '  ELEM MATL',2X,'1-COORD',2X,'2-COORD',
     2 5X,'11-STRESS',5X,'22-STRESS',5X,'33-STRESS',
     3 5X,'12-STRESS')
2010  FORMAT(23X,'EQV STRES',5X,'HYD STRES',5X,'TEMPERATURE',5X,
     1 'FREE VOL')
 2002 FORMAT(I5,1x,i3,2e13.5,4e13.5,/, 4e13.5)
      end
c
c***********************************************************************
c
      subroutine defg(xl,xlc,ul,shp,f,finv,xjac,ndm,ndf,nel,mode,flag)
      implicit double precision (a-h,o-z)
c
c---- compute deformation gradients and jacobian
c
      logical flag
      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL1
      dimension xl(ndm,*),xlc(ndm,*),shp(3,*),f(3,3),finv(3,3),
     1          ul(ndf,*)
c
c---- fij = Na,j uia + deltaij
      do i=1,3
       do j=1,3
        f(i,j) = 0.d0
        finv(i,j) = 0.d0
       end do
        f(i,i) = 1.d0
        finv(i,i) = 1.d0
      end do
      do i=1,2
       do j=1,2
        do l=1,nel
         f(i,j) = f(i,j) + ul(i,l)*shp(j,l)
        end do
       end do
      end do
      IF(MODE.EQ.2) THEN
      RAD = 0.d0
      RADC = 0.d0
      DO  I=1,nel
      RAD = RAD + XL(1,I)*SHP(3,I)
      RADC = RADC + XLC(1,I)*SHP(3,I)
      end do
      F(3,3) = RADC/RAD
      END IF
c---------------------------------------------
c     write(6,9001) ((f(i,j),j=1,3),i=1,3)
9001  format(5x,'f',/,3(3e20.12))
c-----------------------------------
c---- invert deformation gradients
      xjac = f(1,1)*f(2,2) - f(1,2)*f(2,1)
      if (flag) then
         finv(1,1) =   f(2,2)/xjac
         finv(2,2) =   f(1,1)/xjac
         finv(1,2) = - f(1,2)/xjac
         finv(2,1) = - f(2,1)/xjac
      finv(3,3) = 1.d0/f(3,3)
      end if
      xjac = xjac*f(3,3)
C
      if(xjac.le.0.d0) then
       write(6,2001) xjac,nelmt
       stop
      end if
c
      return
2001  format(5x,'xjac ',e13.5,' le 0 in element',i6)
      end
c**********************************************************
      SUBROUTINE BMAT(B,XJC,SHP,XLC,NDM,NDF,NEL,NST,MODE)
      implicit double precision (a-h,o-z)
c
c---- set up B-matrix
c
      dimension xlc(ndm,*),shp(3,*),b(4,*)
      DATA PI/3.141592654/
c
      do 12 i = 1,nst
         do 11 j = 1,4
            b(j,i) = 0.d0
11       continue
12    continue
      i1 = -1
      do 20 i = 1,nel
         i1 = i1 + 2
         i2 = i1 + 1
         b(1,i1) = shp(1,i)
         b(2,i2) = shp(2,i)
         b(4,i1) = shp(2,i)
         b(4,i2) = shp(1,i)
20    continue
      xjc = 1.d0
      IF (MODE.EQ.2) THEN
      RADIUS = 0.d0
      DO 30 I = 1,nel
  30  RADIUS = RADIUS + XLC(1,I)*SHP(3,I)
      XJC = 2.d0*PI*RADIUS
      II = 1
      DO 35 I = 1,nel
      B(3,II) = SHP(3,I)/RADIUS
  35  II = II + 2
      END IF
      return
      end
c*******************************************************************
      subroutine sloop(d,ul,xl,xlc,strl,epsl,ql,ix,tl,s1,p1,p2,vel,
     1                 accel,theta,shpb,bvol,aa,bb,shpu,xsj,
     2                 ndf,ndfm,ndm,nst,nstm,nel,ns,nq,mode,isw)
c
      implicit double precision (a-h,o-z)
c
      dimension d(*),ul(ndf,*),xl(ndm,*),ix(*),tl(*),shpu(3,*),
     1     xlc(ndm,*),strl(*),epsl(*),ql(*),vel(ndf,*),accel(ndf,*)
      dimension b(4,18),dt(6,6),q(4,18),x(2),f(3,3),finv(3,3),
     1          bvol(*),dbvol(18),shpb(3,*),sig(6),pl(18),fp(9),
     2          dtnum(6,6),bbar(4,18),shp(3,9)
      dimension fpnew(9),snew(9),qnew(100),s1(8,*),p1(*),p2(*)
      DATA PI/3.141592654/
c
c---- compute deformation gradients and inverse
      call defg(xl,xlc,ul,shpu,f,finv,xjac,ndm,ndf,nel,mode,.true.)
c---- compute spatial shape functions
      do  i = 1,nel
         shp(1,i) = shpu(1,i)*finv(1,1) + shpu(2,i)*finv(2,1) 
         shp(2,i) = shpu(1,i)*finv(1,2) + shpu(2,i)*finv(2,2) 
         shp(3,i) = shpu(3,i)
      end do
      do i=1,6
       do j=1,6
        dt(i,j) = 0.d0
       end do
      end do
c---- compute shape functions for volumetric deformation
      ratio = (theta/xjac)**(1.d0/3.d0)
      factor = aa + bb*ratio
c---- reduce deformation gradients
      f33 = f(3,3)
      do 307 i = 1,3
         do 307 j = 1,3
            f(i,j) = factor*f(i,j)
307   continue
c---- compute compability matrix (B-bar)
      CALL BMAT(B,XJC,SHP,XLC,NDM,NDF,NEL,NSTM,MODE)
      xsj1 = xsj*xjac
      if(mode.eq.2) xsj1 = xsj*(xjc/f33)*xjac
      xsj2 = xsj
      if(mode.eq.2) xsj2 = xsj*(xjc/f33)
       do i=1,nstm
        do j=1,4
         bbar(j,i) = b(j,i)
        end do
       end do
      factor = bb*ratio/factor
      do 310 i = 1,nstm
         trace = b(1,i) + b(2,i) + b(3,i)
c---- dbvol : (Bbar^{vol} - B^{vol})
         dbvol(i) = (bvol(i) - trace)/3.d0
         prod = factor*dbvol(i)
         do 314 j = 1,3
            bbar(j,i) = b(j,i) + prod
314      continue
310   continue
c
c---- compute spatial tangent stiffness  (dt)
      if ((isw.eq.3).or.(isw.eq.7)) then
c---- dtang is model specific subroutine 
c      call update(d(6),f,strl,epsl,ql,dt,nq,ns,lquad,.true.)
c       call dtang(d(6),f,dt,ql,strl,epsl)
      do 330 i = 1,6
         sig(i) = xsj1*strl(i)
330   continue
c---------------------------------------------------------------
c        write(6,9003) ((dt(i,j),i=1,6),j=1,6)
9003     format(5x,'dt',6e13.5)
c        write(6,9004) xsj1
9004     format(5x,'xsj1=',e13.5)
c-----------------------------------------------------------
         do  i = 1,4
          do j=1,4
           dt(i,j) = dt(i,j)*xsj1
          end do
         end do
c
c---- compute Bbar^T*DT*Bbar
         do  i = 1,nstm
          do j=1,4
           q(j,i) = 0.d0
           do l=1,4
            q(j,i) = q(j,i) + dt(j,l)*bbar(l,i)
           end do
          end do
         end do
c
         do i=1,nstm
          do j=1,nstm
           do l=1,4
            s1(i,j) = s1(i,j) + bbar(l,i)*q(l,j)
           end do
          end do
         end do
c
c---- add regular geometric stiffness term (Moran et.al., 1990) :
c     tau: Kirchhoff stress; del: Kronecker delta; (a,b): Typical nodes
c     (i,j,k,l): 1,2,3 ;  Na,j : Spatial gradients of shape function
c     K_{iakb} = int_Vo^e {del_{ik}*tau_{jl}*Na,j*Nb,l} dVo^e
c---- Note order of sig : s11,s22,s33,s12,s13,s23
c
c       if(factor.gt.0.d0) then
c-----------------------------------------------
c       write(6,9005)
9005    format(5x,'***Entering geometric nonlinear calculn in isw=3')
c----------------------------------------------------------------------
      IF(MODE.EQ.2)THEN
      RAD = 0.d0
      rad0 = 0.d0
      DO 326 I=1,nel
      RAD = RAD + XLC(1,I)*SHP(3,I)
        rad0 = rad0 + xlc(1,i)*shpb(3,i)
 326    CONTINUE
      END IF
         ii = 1
         sb3 = 0.d0
         do 335 i = 1,4
            i1 = ii + 1
            sb1 = sig(1)*shp(1,i) + sig(4)*shp(2,i)
            sb2 = sig(4)*shp(1,i) + sig(2)*shp(2,i)
      IF (MODE.EQ.2)SB3 = SIG(3)*SHP(3,I)/(RAD*RAD)
            bsb = shp(1,1)*sb1 + shp(2,1)*sb2
            s1(ii,1)   = s1(ii,1)   + bsb  + sb3*shp(3,1)
            s1(i1,2)   = s1(i1,2)   + bsb
            bsb = shp(1,2)*sb1 + shp(2,2)*sb2
            s1(ii,3)   = s1(ii,3)   + bsb  + sb3*shp(3,2)
            s1(i1,4)   = s1(i1,4)   + bsb
            bsb = shp(1,3)*sb1 + shp(2,3)*sb2
            s1(ii,5)   = s1(ii,5)   + bsb  + sb3*shp(3,3) 
            s1(i1,6)   = s1(i1,6)   + bsb
            bsb = shp(1,4)*sb1 + shp(2,4)*sb2
            s1(ii,7)   = s1(ii,7)   + bsb  + sb3*shp(3,4)
            s1(i1,8)   = s1(i1,8)   + bsb
            ii = ii + 2
335      continue
c        end if
c        ii = 1
c        sb3 = 0.d0
c        do  i = 1,nel
c           sb1 = sig(1)*shp(1,i) + sig(4)*shp(2,i) 
c           sb2 = sig(4)*shp(1,i) + sig(2)*shp(2,i)
c     IF(MODE.EQ.2) SB3 = SIG(3)*SHP(3,I)/(RAD*RAD)
c           JJ = 1
c         DO  J=1,NEL
c            BSJ = SHP(1,J)*SB1 + SHP(2,J)*SB2 
c          S(II,JJ)  = S(II,JJ)  + BSJ + SB3*SHP(3,J)
c            S(II+1,JJ+1)  = S(II+1,JJ+1)  + BSJ 
c            JJ = JJ+NDF
c         end do
c           II = II + NDF
c        end do
c
c---- add geometric terms due to reduction (assumed strain)
c---- Note factor =1 for stabilization parameter epsilon = 0. In this
c     case (Moran et.al.(1990)) :
c---- K = 2/3{int_Vo^e {(B^T*tau)x(Bbar^{vol}-B^{vol}) +(Bbar^{vol}-B^{vol})x(B^T*tau)
c                        - tau_{kk}/3*(Bbar^{vol}-B^{vol})x(Bbar^{vol}-B^{vol})}dVo^e}
c
        if(factor.gt.0.d0) then
c-----------------------------------------------
c       write(6,9001)
9001    format(5x,'***Entering assumed strain calculation in isw=3')
c----------------------------------------------------------------------
         xm1 = 2.d0*factor
         xm2 = factor*(3.d0*factor - 1.d0)
         trace = sig(1) + sig(2) + sig(3)
c---- Form pl = Bbar^T*tau
         do  i = 1,nstm
          pl(i) = 0.d0
          do j=1,4
           pl(i) = pl(i) + bbar(j,i)*sig(j)
c          pl(i) = pl(i) + b(j,i)*sig(j)
          end do
         end do
         do i=1,nstm
          do j=1,nstm
           s1(i,j) = s1(i,j) + xm1*(pl(i)*dbvol(j) + dbvol(i)*pl(j))
     1                     - xm2*trace*dbvol(i)*dbvol(j) 
          end do
         end do
c
c---- Add foll term : K_{iakb} = int_{Vo^e} {tau_{kk}/3.0*(Na,k*Nb,i 
c                                            - Na,k*Nb,i|_{reduced GP})} dVo^e
         pp = factor*trace/3.d0
         ii = 0
         do i=1,nel
          do k=1,ndfm
          jj = 0
           do j=1,nel
c----------- foll are trial statements ----------------
          if(mode.eq.2) then
           dum1 = shp(3,i)*shp(3,j)/(rad*rad)
           dum2 = shpb(3,i)*shpb(3,j)/(rad0*rad0)
          end if
c----------------------------------------------------------
            do l=1,ndfm
             ssb = shpb(k,j)*shpb(l,i)
             s1(ii+k,jj+l) = s1(ii+k,jj+l) + pp*(shp(k,j)*shp(l,i) - 
     1       ssb)
c----------- foll are trial statements ----------------
             if(mode.eq.2) then
             if((k.eq.1).and.(l.eq.1)) then
              s1(ii+k,jj+l) = s1(ii+k,jj+l) + pp*(dum1-dum2)
             end if
             end if
c--------------------------------------------------------------
            end do
            jj = jj+ndfm
           end do
          end do
          ii = ii+ndfm
         end do
c
      end if
      end if
c---- compute internal forces : strl : Cauchy stress
      if ((isw.eq.6).or.(isw.eq.7)) then
      do  i = 1,6
         sig(i) = xsj1*strl(i)
         end do
c----   compute gravity loads : D(1), D(2) are body forces / ref volume
         if ((d(1).ne.0.d0).or.(d(2).ne.0.d0)) then
            ii = 1
            do 605 i = 1,nel
               i1 = ii + 1
               w = shp(3,i)*xsj2
               p1(ii) = p1(ii) + d(1)*w
               p1(i1) = p1(i1) + d(2)*w
               ii = ii + ndfm
605         continue
         end if
c---- subtract internal forces
         do  i = 1,nstm
          do j=1,4
            p1(i) = p1(i) - bbar(j,i)*sig(j)
          end do
         end do
c
c----- obtain force vector associated with free volume evolution
        call pxeta(d(6),ql,p2,shpu,ul,xsj2,ndf,nel)
      end if
      return
      end
c
c***********************************************************************
c
      subroutine uloop(d,ul,xl,xlc,strl,epsl,ql,ix,tl,vel,accel,
     & theta,aa,bb,shp,xsj,xsj1,ndf,ndm,nst,nel,ns,nq,mode,lquad)
      implicit double precision (a-h,o-z)
c
      common /freevolgrad/ xigrad(2,5000)
      dimension d(*),ul(ndf,*),xl(ndm,*),ix(*),tl(*),shp(3,*),
     1          xlc(ndm,*),strl(*),epsl(*),ql(*),vel(ndf,*),accel(ndf,*)
      dimension f(3,3),finv(2,2),sig(6),fp(9),dt(6,6)
c
c---- compute deformation gradients and jacobian
      call defg(xl,xlc,ul,shp,f,finv,xjac,ndm,ndf,nel,mode,.false.)
      xsj1 = xsj*xjac
      if(mode.eq.2) xsj1 = xsj*(xjac/f(3,3))
c---- compute shape functions for volumetric deformation
      ratio = (theta/xjac)**(1.d0/3.d0)
      factor = aa + bb*ratio
c---- reduce deformation gradients
      do 962 i = 1,3
         do 961 j = 1,3
            f(i,j) = factor*f(i,j)
961      continue
962   continue
        do i=1,4
         sig(i) = strl(i)
        end do
c---- store xi_n1 in ql(8) & transfer xi_n to ql(7)
      ql(7) = ql(8)
      ql(8) = 0.d0
c---- Compute Del^2 (xi) - Laplacian of xi wrt material coords
      del2xi = 0.d0
      do i=1,nel
       ii = ix(i)
       ql(8) = ql(8) + ul(3,i)*shp(3,i)
       del2xi = del2xi + shp(1,i)*xigrad(1,ii) + shp(2,i)*xigrad(2,ii)
      end do
c---- update is model specific subroutine 
c     call update(d,f,strl,epsl,ql,dt,nq,ns,lquad,.false.)
c     call update1(d,f,strl,epsl,ql,nq,ns)
 	call update(d,f,strl,epsl,ql,del2xi,nq,ns)
      return
      end
c---------------------------------------------------------------
      SUBROUTINE STFSTORE(S,S1,S2,P,P1,P2,NDFM,NDF,NEL,NST,IOP)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C---- STORES MECHANICAL AND THERMAL STIFFNESS MATRIX AND FORCE VECTORS
C---- IN GLOBAL MATRIX S AND GLOBAL ARRAY P
      DIMENSION S(NST,*), S1(8,8),S2(4,4),
     1 P1(8),P2(4),P(*)
      IF((IOP.EQ.1).OR.(IOP.EQ.3)) THEN
      IF (NDF .EQ. 3) THEN
      DO 15 I=1,NEL
      DO 15 K=1,NEL
      DO 10 J=1,2
      DO 10 L=1,2
      S(3*(I-1)+J,3*(K-1)+L) = S1(2*(I-1)+J,2*(K-1)+L)
10    CONTINUE
      S(3*I,3*K) = S2(I,K)
15    CONTINUE
      ELSEIF (NDF .EQ. 2) THEN
      DO 20 I=1,NST
      DO 20 J=1,NST
      S(I,J) = S1(I,J)
20    CONTINUE
      ELSE
      DO 40 I=1,NST
      DO 40 J=1,NST
      S(I,J) = S2(I,J)
40    CONTINUE
      END IF
      END IF
      IF((IOP.EQ.2).OR.(IOP.EQ.3)) THEN
      IF (NDF .EQ. 3) THEN
      DO 30 I=1,NEL
      DO 25 J=1,2
      P(3*(I-1)+J) = P1(2*(I-1)+J)
25    CONTINUE
      P(3*I) = P2(I)
30    CONTINUE
      ELSEIF (NDF .EQ. 2) THEN
      DO 35 I=1,NST
      P(I) = P1(I)
35    CONTINUE
      ELSE
      DO 45 I=1,NST
      P(I) = P2(I)
45    CONTINUE
      END IF
      END IF
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine pxeta(d,ql,p2,shpu,ul,xsj2,ndf,nel)
      implicit double precision (a-h,o-z)
      dimension d(*),ql(*),p2(*),shpu(3,*),ul(ndf,*),sxi(4,4)
c----- following are material parameters d(1) to d(20)  in E-T model
      ekappa = d(1)
      emu = d(2)
      alfa_th = d(3)
      theta_o = d(4)
      f_o = d(5)
      ek_c = d(6)
      gama_dot_o = d(7)
      a_rexp = d(8)
      phi = d(9)
      xi_o = d(10)
      xeta = d(11)      
      theta_r = d(12)
      xi_r = d(13)
      ek_theta = d(14)
      s1 = d(15)
      s2 = d(16)
      s3 = d(17)
c
c---- following are the internal variables
      xi_dot_o = ql(1)
      gama_dot = ql(2)
      tau = ql(3)
      p = ql(4)
      theta_n = ql(5)
      del_theta_n = ql(6) - ql(5)
      theta_n1 = ql(6)
      xi_n = ql(7)
      del_xi_n = ql(8) - ql(7)
      xi_n1 = ql(8)
      tree = ql(9)
      del_trE_e = ql(10)
      xi_T = ql(20)
      coh = ql(22)
c     
c---- form [K_xi] matrix
      do i =1,4
        do j=1,4
         sxi(i,j) = 0.d0
        end do
      end do
      const1 = s1*xi_dot_o/s3
      const2 = s2*xi_dot_o/s3
      do i=1,nel
       do j=i,nel
        sxi(i,j) = sxi(i,j) + const1*(shpu(1,i)*shpu(1,j) +
     1   shpu(2,i)*shpu(2,j)) + shpu(3,i)*shpu(3,j)*const2 
       end do
      end do
      do i=1,nel
       do j=i,nel
        sxi(j,i) = sxi(i,j)
       end do
      end do
c
c---- form [K_xi]*{XI}
      do i=1,nel
       do j=1,nel
        p2(i) = p2(i) - sxi(i,j)*ul(3,j)
       end do
      end do
c
c---- form {F_xi^int}
      const3 = gama_dot*xeta - xi_dot_o*(p-s2*xi_T)/s3
      do i=1,nel
        p2(i) = p2(i) + shpu(3,i)*const3
      end do
      do i=1,nel
       p2(i) = p2(i)*xsj2
      end do
c
      return
      end

