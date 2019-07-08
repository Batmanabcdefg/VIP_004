      subroutine elmt03(d,ul,xl,strl,epsl,ql,ix,tl,s,p,
     1 vel,accel,ndf,ndm,nst,ns,nq,isw)
      implicit double precision (a-h,o-z)
c
c---- this element routine is written to do leat squares smoothing
c---- for the element variables that are obtained as an output from
c---- a finite element program.
c
c---- **LEAST SQUARES POST-PROCESSING ROUTINE FOR FEAP**
c
c---- The program smooths the variable ql(1,i) which are obtained
c---- at quadrature points from some finite element calculations.
c---- eg. ql(1,i)= E33(i)--- out of plane displacement
c---- see nodal force section isw=6
c---- best way to input data for this program is to use the macro ivar
c---- in the input data file--- write the variable to be smoothed, say
c---- the E33 strain at all quadrature points in the format for ivar
c---- and input this in the new data file for this element routine.
c---- mesh data in the input file same as the original data file
c---- used to create the ivar data
c
c---- Reference: E. Hinton and J. S. Cambell., `` Local and Global
c----            Smoothing of Discontinous Finite Element Functions Using a
c----            Least Squares Method,'' Intenational Journal of Numerical
c----            Methods in Engineering, Vol. 8, 1974, pp. 461-480.
c
      character*7 label
      logical flag
      common /chead/ o,head(20)
      common /cdata/ numnp,numel,nummat,
     1 nen,neq,ipr,nsdm,nqdm,nquad
      common /eldata/ dm,nelmt,ma,mct,iel,nel
c     common /slipvecs/ xm(3,12),xs(3,12),schmid(3,3,12),tauo(12)
      dimension d(*),ul(ndf,*),xl(ndm,*),ix(*),tl(*),s(nst,*),
     1 p(*),vel(*),accel(*),shp(3,9),sg(16),tg(16),wg(16),xlc(2,4)
      dimension strl(ns,*),epsl(ns,*),ql(nq,*),x(3),b(2,18),qvar(2)
     1 ,svar(4),shpb(3,9),xsc(3),fe(3,3),s1(8,8),s2(4,4),p1(8),p2(4)
      dimension voigt1(6),fp(3,3),fpi(3,3),cp(3,3),f(3,3),finv(3,3)
      data pi/3.141592653/
c
      ndfm = 2
      nstm = ndfm*nen
c---- go to correct array processor
      go to (1,2,3,4,5,6,2,8,9,10,11), isw
c
C---- input material properties
c---- it is essential that you have  mesh macro 'MATE'
c----MATE
c----    1    1
c---- id = any number, or you can read some material data you need to use
    1 continue
      return
c
    2 continue
      return
c
c---- compute stiffness matrix
    3 continue
      mode = d(3)
c     mode = 1
      do i = 1,8
       do j=1,8
        s1(i,j) = 0.d0
       end do
      end do
      do i=1,4
       do j=1,4
        s2(i,j) = 0.d0
       end do
       s2(i,i) = 1.d0
      end do
c---- find coords in current config
      do 301 i=1,2
      do 301 j=1,nel
      xlc(i,j) = xl(i,j) + ul(i,j) 
301   continue
c---- set up integration points and weights
      nint = 2
      if (nel.gt.4) nint = 3
      call pgaus2(nint,lint,sg,tg,wg)
c---- loop over integration points
      do 300 l = 1,lint
c---- compute shape functions and derivatives
      call shap2d(sg(l),tg(l),xlc,shp,xsj,ndm,nel,ix,.false.)
      xsj = xsj*wg(l)
      if(mode.eq.2) then
      rad = 0.d0
      do 302 i=1,4
      rad = rad + shp(3,i)*xlc(1,i)
302   continue
      xsj = xsj*rad
      end if
      do 310 i=1,nst
      b(1,i) = 0.d0
      b(2,i) = 0.d0
310   continue
      do 315 i=1,nel
      b(1,2*i-1) = shp(3,i)
      b(2,2*i)   = shp(3,i)
315   continue
c
c---- compute NT*N
      do 326 i = 1,nstm
      do 326 j = 1,nstm
      do 325 k = 1,2
      s1(i,j) = s1(i,j) + b(k,i)*b(k,j)*xsj
  325 continue
  326 continue
  300 continue
      CALL STFSTORE(S,S1,S2,P,P1,P2,NDFM,NDF,NEL,NST,1)
      return
c
c---- output element variables
    4 continue
c---- set up output points and shape functions
      id1 = d(1)
      id2 = d(2)
c     nint = 2
c     if (nel.gt.4) nint = 3
c     call pgaus2(nint,lint,sg,tg,wg)
c     ll = 0
      ll=4
c     do 400 l = 1,lint
      ll = ll + 1
c     call shap2d(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.true.)
      call shap2d(0.d0,0.d0,xl,shp,xsj,ndm,nel,ix,.true.)
c---- compute output point coordinates
      do 410 i = 1,ndm
      x(i) = 0.d0
      qvar(i) = 0.d0
      svar(i) = 0.d0
      svar(2+i) = 0.d0
      do 410 j = 1,nel
      x(i) = x(i) + shp(3,j)*xl(i,j)
      if(id1.le.4) then
      svar(i) = svar(i) + shp(3,j)*ul(i,j+nen)
      svar(2+i) = svar(2+i) + shp(3,j)*vel(2*(j-1)+i)
      else
      qvar(i) = qvar(i) + shp(3,j)*vel(2*(j-1)+i)
      end if
410   continue
c---- sv is the value of the variable that has been smoothed.
c---- svx and svy are the x,y derivatives of the value sv.
c---- compute output variable
c     sv =0.d0
c     svx = 0.d0
c     svy = 0.d0
c     do 420 j = 1,nel
c     sv = sv + shp(3,j)*ul(j)
c     svx = svx + shp(1,j)*ul(j)
c20   svy = svy + shp(2,j)*ul(j)
c---- print variables
c     write(7,2020) (x(i),i=1,ndm),sv,svx,svy
      if(id1.le.4) then
      snrm1 = 0.d0
      snrm2 = 0.d0
      do 411 k=1,3
      snrm1 = snrm1 + svar(k)*svar(k)
      snrm2 = snrm2 + strl(k,ll)*strl(k,ll)
411   continue
      snrm1 = snrm1 + 2.d0*svar(4)*svar(4)
      snrm2 = snrm2 + 2.d0*strl(4,ll)*strl(4,ll)
      snrm1 = dsqrt(snrm1)
      snrm2 = dsqrt(snrm2)
      write(6,2020) (x(i),i=1,ndm),snrm1,snrm2
      end if
      if((id1.gt.4).and.(id2.gt.4)) then
      write(6,2020) (x(i),i=1,ndm),qvar(1),qvar(2),ql(id1-4,ll),
     1 ql(id2-4,ll)
      end if
  400 continue
      return
c
c---- compute consistent and lumped mass matrices
    5 continue
      return
c
c---- compute internal forces
    6 continue
      id1 = d(1)
      id2 = d(2)
      mode = d(3)
c     mode = 1
c---- find coord in current config.
      do 601 i=1,2
      do 601 j=1,nel
      xlc(i,j) = xl(i,j) + ul(i,j)
601   continue 
c
c---- set up integration points
      nint = 2
c     nint = 1
      if (nel.gt.4) nint = 3
      call pgaus2(nint,lint,sg,tg,wg)
c---- loop over integration points
      do 600 l = 1,lint
c---- compute shape functions and derivatives
      call shap2d(sg(l),tg(l),xlc,shp,xsj,ndm,nel,ix,.false.)
      xsj = xsj*wg(l)
      if(mode.eq.2) then
      rad = 0.d0
      do 602 i=1,4
      rad = rad + shp(3,i)*xlc(1,i)
602   continue
      xsj = xsj*rad
      end if
c
c---- determine fp and maximum plastic log strain
c---- retrieve Fn^p{inv} from eps (single crystal plasticity routines)
      if(id1.gt.200) then
c
c     call invoigtv(fpi,voigt1,epsl(1,l),2)
c     call invm3x3(fpi,fp,flag)
c      if(.not.flag) then
c       write(6,2030) nelmt
c       stop
c      end if
c
      call invoigtv(fp,voigt1,ql(11,l),2)
c---- Form Cp = Fp^T*Fp 
      do i=1,3
       do j=i,3
        cp(i,j) = 0.d0
         do k=1,3
          cp(i,j) = cp(i,j) + fp(k,i)*fp(k,j)
         end do
        end do
       end do
      do i=1,3
       do j=i,3
        cp(j,i) = cp(i,j)
       end do
      end do
c---- find eigenvalues of cp (2D case)
      detcp = cp(1,1)*cp(2,2) - cp(1,2)*cp(2,1)
      trcp = cp(1,1) + cp(2,2)
      discr = dsqrt(trcp*trcp - 4.d0*detcp)
      eig1 = 0.5*(trcp + discr)
      eig2 = 0.5*(trcp - discr)
      if((eig1.lt.0.d0).or.(eig2.lt.0.d0)) then
       write(6,2040) nelmt
      stop
      end if
      pstr1 = dsqrt(eig1)
      pstr2 = dsqrt(eig2)
      peplog1 = dlog(pstr1)
c
c---- compute deformation jacobian
c     call shap2d(sg(l),tg(l),xl,shpb,xsjo,ndm,nel,ix,.true.)
c     call defg(xl,xlc,ul,shpb,f,finv,xjac,ndm,ndf,nel,mode,.false.)
c----  obtain fe = f*fpinv
c     do i = 1,3
c      do j= 1,3
c       fe(i,j) = 0.d0
c        do k=1,3
c         fe(i,j) = fe(i,j) + f(i,k)*fpi(k,j)
c        end do
c       end do
c      end do
c---- Find current orientation of slip vector for slip system #id3 (2Dcase) sc = Fe*s
c      id3 = id1 - 200
c     do k=1,2
c      xsc(k) = 0.d0
c       do j=1,2
c        xsc(k) = xsc(k) + fe(k,j)*xs(j,id3)
c       end do
c      end do
c      xscmag = dsqrt(xsc(1)*xsc(1) + xsc(2)*xsc(2))
c      if(xscmag.eq.0.d0) then
c       write(6,2050) nelmt
c       stop
c      end if
c      xsc(1) = xsc(1)/xscmag
c      xsc(2) = xsc(2)/xscmag
c----- Find angle (in degrees) beween current slip vector & the undeformed confign.
c      ang = dacos(xsc(1)*xs(1,id3) + xsc(2)*xs(2,id3))
c      ang = ang*180.d0/pi
c
       end if
c
c---- compute applied force vector
      do 605 i = 1,nel
      w = shp(3,i)*xsj
      if(id1.ne.0) then
      if(id1.le.100) then
      p(3*i-2) = p(3*i-2) + strl(id1,l)*w
      end if
      if(id1.gt.100) then
      if(id1.le.200) then
      id3 = id1 - 100
c---- internal variables to be plotted
      p(3*i-2) = p(3*i-2) + ql(id3,l)*w
      else
c---- max log plastic strain to be plotted (finite strain case)
      p(3*i-2) = p(3*i-2) + peplog1*w
c     id3 = id1 - 13
c     p(2*i-1) = p(2*i-1) + epsl(id3,l)*w
      end if
      end if
      end if
      if(id2.ne.0) then
      if(id2.le.100) then
      p(3*i-1) = p(3*i-1) + strl(id2,l)*w
      end if
      if(id2.gt.100) then
      if(id2.le.200) then
      id3 = id2 - 100 
      p(3*i-1) = p(3*i-1) + ql(id3,l)*w
      else
       id3 = id2 - 200
      p(3*i-1) = p(3*i-1) + ql(id3,l)*w
c     id3 = id2 - 13
c
c     p(2*i) = p(2*i) + ang*w
      end if
      end if
      end if
  605 continue
  600 continue
c
c---- unused isw
    8 continue
c---- compute coords of 1x1 g.p. and store in eps aray
c---- find coords in current config
c     do 801 i=1,2
c     do 801 j=1,nel
c     xl(i,j) = xl(i,j) + ul(i,j) 
801   continue
      call shap2d(0.d0,0.d0,xl,shp,xsj,ndm,nel,ix,.false.)
c---- compute output point coordinates
      do 810 i = 1,ndm
      x(i) = 0.d0
      do 810 j = 1,nel
      x(i) = x(i) + shp(3,j)*xl(i,j)
810   continue
      epsl(1,5) = x(1)
      epsl(2,5) = x(2)
c     do 811 i=1,2
c     do 811 j=1,nel
c     xl(i,j) = xl(i,j) - ul(i,j) 
811   continue
      return
c
c---- unused isw
    9 continue
      return
c
c---- unused isw
   10 continue
      return
c
c---- unused isw
   11 continue
      return
c
c---- formats
2010  format(i5)
2020  format(6(1x,1pe13.5))
2030  format(5x,'Cannot invert Fpinv in smth for elem',i6)
2040  format(5x,'Negative e-value for CP in smth for elem',i6)
2050  format(5x,'Magnitude of slip vector zero in smth for elem',i6)
      end