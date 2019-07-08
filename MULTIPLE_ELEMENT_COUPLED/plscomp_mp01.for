c---- MP01 macro file for plane strain compression specimen
      subroutine mp01(idl,ie,d,id,x,ix,f,t,dr,str,eps,ql,vel,accel,
     1 ndf,ndm,nen1,nstr,nq,iii,prtn)
      implicit double precision (a-h,o-z)
c
      logical prtn,prt
      character*8 head
      common /cdata/ numnp,numel,nummat,nen,neq,ipr,
     1 nsdm,nqdm,nquad                                                          !control card data
      common /bndrcntr/ nbncpt, ibnc(20)
      common /pmprtn/ prt
      common /prtset/ nset(9,201)                                               !set card data
c
      dimension ie(7,*),d(40,*),id(ndf,*),x(ndm,*),ix(nen1,*),f(ndf,*),
     1 t(*),dr(ndf,*),str(nstr,*),eps(nstr,*),ql(nq,*),
     2 vel(ndf,*),idl(*),accel(ndf,*)
      dimension coh(3000)
      data small /0.4d-11/,temp0/295.d0/,xi0/0.0006/
      data xmax /10.d-09/, ymax /20.d-09/, scale /1.d0/
c
c---- read cohesion data
      open(unit=38,file='cohes.dat',form='formatted')
      nel =4
      read(38,1001) (coh(i),i=1,numel)
	write(6,1001) (coh(i),i=1,numel)                                          !read the values of cohesion for each of the 800 elements
1001  format(8f10.4)
c
c---- scale up x, y coordinates
      do i = 1,numnp                                                            !numnp is total number of nodes
       x(1,i) = x(1,i)*scale                                                    
       x(2,i) = x(2,i)*scale                                                    !x= coordinates of all nodes(dimension ndm*numnp)
      end do
      xmax1 = xmax*scale
      ymax1 = ymax*scale
c---- set internal variables for each element
      do 20 i=1,numel                                                           !numel is the number of elements
      do 21 j = 1,nquad                                                         !nquad is the number of Gauss points
c---- first initialize all internal vars at the output Gauss point to be zero
      do 22 l = 1,nqdm                                                          !nqdm is the number of internal variables per Gauss point
      ql((j-1)*nqdm+l,i) = 0.d0
22    continue 
      ql((j-1)*nqdm+5,i) = temp0 
      ql((j-1)*nqdm+6,i) = temp0 
      ql((j-1)*nqdm+7,i) = xi0                                         !???? isingh : Why two Xi0 ? 
      ql((j-1)*nqdm+8,i) = xi0 
      ql((j-1)*nqdm+11,i) = 1.d0 
      ql((j-1)*nqdm+12,i) = 1.d0 
      ql((j-1)*nqdm+13,i) = 1.d0 
      ql((j-1)*nqdm+20,i) = xi0 
      ql((j-1)*nqdm+21,i) = coh(i)*1.d+09 
      ql((j-1)*nqdm+22,i) = coh(i)*1.d+09                                        !set the initial conditions for the internal variables at the output GP 
21    continue
20    continue
c---- set up boundary conditions and forces
      call bounplscomp(id,ix,x,f,xmax1,ymax1,numnp,numel,ndf,ndm,nen1)
c
      return
1000  format(i9,3f15.6,2i9)
1100  format(5i9)
1110  format(5x,'No of boundary points exceeds 20', i5)
1112  format(15i5)
      end
c
      subroutine norder(x,ix,mset,xstart,idof,kset,
     1 ndm,nen1,numel)
      implicit double precision (a-h,o-z)
      common /prtset/ nset(9,201)                                               !total of 9 element/node sets with 200 entities each
      dimension ix(nen1,*),x(ndm,*)
      dimension idum(200),theta(200)
      data tol /1.0d-06/
      npset = nset(mset,1)                                                      !npset is the number of nodes in the current node set, mset is the
      xmin = 1.0d+12                                                            !node set number
      do 100 i=1,npset
      nn = nset(mset,i+1)                                                       !nn is the entity number of the current entity
c     if(abs(x(idof,nn) - xstart).lt.tol) then
      if(x(idof,nn).lt.xmin) then
      idum(1) = nn                                                              !at the end of do loop, idum(1) contains the entity number of the entity
      xmin = x(idof,nn)                                                         !with the least x coordinate
      end if                                                                    
c     ndum = 1                                                                  
c     go to 110
c     end if
100    continue
c     write(6,2001) mset
c     stop
      ndum = 1
110    continue
      xmax = 1.0d+12
      nprev = idum(ndum)
      do 150 i=1,npset
      nn = nset(mset,i+1)
      if((x(idof,nn).lt.xmax).and.(x(idof,nn).gt.x(idof,nprev))) then
      idum(ndum+1) = nn
      xmax = x(idof,nn)
      end if
150    continue                                                                 !idum(2) contains the entity number of the entity with the second least
      ndum = ndum+1                                                             !x coordinate and so on...
      if(ndum.lt.npset) go to 110                                               
      do 200 i=1,npset
      nset(mset,i+1) = idum(i)                                                  !ordering the entity numbers in ascending order of x coordinate values
200    continue
      if(kset.eq.0) return
c---- locate elements that lie along the node line given by node set #mset
c---- and store these in kset.
      kpset = 0
      do 30 i = 1,npset
      nn = nset(mset,i+1)                                                       !nn is entity number of the current entity of the current node set mset
      do 35 n = 1,numel
      do 40 k=1,4
      if(ix(k,n).eq.nn) then                                                    !if the global node number of kth node of nth element is equal to nn
      if(i.ne.npset) then
      do 41 l = 1,4
      if(ix(l,n).eq.nset(mset,i+2)) then                                        !if the global node number of lth node of nth element is equal to the 
      kpset = kpset + 1                                                         !next entity number
      nset(kset,kpset+1) = n                                                    !node set number takes the value of the element number
      go to 35
      end if
41    continue
      else
      do 42 l = 1,4
      if(ix(l,n).eq.nset(mset,i)) then
      kpset = kpset + 1
      nset(kset,kpset+1) = n                                                    !same statements as above repeated
      go to 35
      end if
42    continue
      end if
      end if
40    continue
35    continue
30    continue
      nset(kset,1) = kpset
      return
2001  format(5x,'cannot find starting node of set',i5)
      end
c
      subroutine bounplscomp(id,ix,x,f,xmax,ymax,numnp,numel,ndf,
     1 ndm,nen1)
c-----Set up boundary conditions for plane strain compression specimen
      implicit double precision (a-h,o-z)
c
      logical prt
      common /pmprtn/ prt
      common /prtset/ nset(9,201)
       dimension id(ndf,*),ix(nen1,*),x(ndm,*),f(ndf,*)
      data small /0.4d-11/,temp0/295.d0/
c     data xmax /10.d0/, ymax /20.d0/
       IF (PRT)
     1 write(6,2000)
        iset1 = 1
	do 10 i=1,numnp
        temp1 = dabs(x(2,i) - 0.d0)
        temp2 = dabs(x(2,i) - ymax)
        temp3 = dabs(x(1,i) - 0.d0)
        temp4 = dabs(x(1,i) - xmax)
c----- symmetry bcs on nodes along x-axis
      if(temp1.lt.small) then
       id(2,i) = 1
c	 id(4,i) = 1
c	 f(4,i) = temp0                                                              !nodes very close to y=0.d0 have y-displacement specified
c     if(prt)
c    1 write(6,2002) i,x(1,i),x(2,i),id(1,i),id(2,i),f(1,i),f(2,i)
c      end if
c----- fix node at x=0.5*xmax, y =0
      if(temp3.lt.small) then
       id(1,i) = 1                                                              !nodes very close to x=0.d0 have x-displacement specified
       end if
      if(prt)
     1 write(6,2002) i,x(1,i),x(2,i),id(1,i),id(2,i),id(3,i),id(4,i),
     2 f(1,i),f(2,i),f(3,i),f(4,i)
       end if
c----- specified u_2 bcs on nodes along top boundary
      if(temp2.lt.small) then
       id(2,i) = 1                                                              !nodes very close to y=ymax have y-displacement specified and
       f(2,i)= -1.d0
c	 id(4,i) =1
c	 f(4,i) = temp0                                                            !y displacements are specified to be -1.d0
       iset1 = iset1 + 1
c----- fix node at x=0.5*xmax, y =0
c      if(temp3.lt.small) then
c      id(1,i) = 1                                                              !nodes very close to x=xmax/2.d0 have x-displacement specified
c      end if
c---- store nodes on top boundary in set 1
      nset(1,iset1) = i                                                         !the first index refers to the node set number
      if(prt)
     1 write(6,2002) i,x(1,i),x(2,i),id(1,i),id(2,i),id(3,i),id(4,i),
     2 f(1,i),f(2,i),f(3,i),f(4,i)
       end if
c	if(temp3.lt.small) then
c      id(4,i) = 1
c	f(4,i) = temp0
c	if(prt)
c     1 write(6,2002) i,x(1,i),x(2,i),id(1,i),id(2,i),id(3,i),id(4,i),
c     2 f(1,i),f(2,i),f(3,i),f(4,i)
c       end if
c	if(temp4.lt.small) then
c      id(4,i) = 1
c	f(4,i) = temp0
c	if(prt)
c     1 write(6,2002) i,x(1,i),x(2,i),id(1,i),id(2,i),id(3,i),id(4,i),
c     2 f(1,i),f(2,i),f(3,i),f(4,i)
c       end if
10    continue
      nset(1,1) = iset1 - 1                                                     !(iset1-1) is the number of nodes along the top boundary
c---- order nodes  along x-axis
c---- in increasing order with respect to x-coordinate
      call norder(x,ix,1,0.d0,1,0,ndm,nen1,numel)
c---- print set information
c     do 25 mset=1,6
      mset = 1
      npset = nset(mset,1)
      IF(PRT) THEN
      WRITE(6,2003) mset,(NSET(mset,NN),NN=1,npset+1)
      end if
25    continue
      return
2000  format(5x,'BOUNDARY CONDITIONS AND FORCES',/,
     1 2x,'Node',7x,'X',12x,'Y',3x,'ID(1)',2x,'ID(2)',2x,'ID(3)',2x,
     2 'ID(4)',10x,'FX',10x,'FY',10x,'FXI',10x,'FTHETA')
2002  format(i5,2e13.5,4i3,4e13.5)
2003  FORMAT(/5X,' SET NUMBER',I3,' HAS TOTAL OF',I4,' ELEMENTS OR ',
     1 'NODES'/14(/15I6))
2004  FORMAT(i5,2e13.5)
      end
c
      SUBROUTINE SHAP2D1(SS,TT,X,SHP,XSJ,NDM,NEL,IX,FLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FLAG
C
C---- SHAPE FUNCTION ROUTINE FOR TWO DIMENSIONAL ELEMENTS
C
      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL1
      DIMENSION SHP(3,*),X(NDM,*),S(4),T(4),XS(2,2),SX(2,2),IX(*)
      DATA S/-0.5,0.5,0.5,-0.5/,T/-0.5,-0.5,0.5,0.5/
C---- FORM 4-NODE QUADRILATERAL SHAPE FUNCTIONS
      DO 100 I = 1,4
      SHP(3,I) = (0.5+S(I)*SS)*(0.5+T(I)*TT)
      SHP(1,I) = S(I)*(0.5+T(I)*TT)
100   SHP(2,I) = T(I)*(0.5+S(I)*SS)
      IF(NEL.GE.4) GO TO 120
C---- FORM TRIANGLE BY ADDING THIRD AND FOURTH TOGETHER
      DO 110 I = 1,3
110   SHP(I,3) = SHP(I,3)+SHP(I,4)
C---- ADD QUADRATIC TERMS IF NECESSARY
120   IF(NEL.GT.4) CALL QUAD(SS,TT,SHP,IX,NEL)
C---- CONSTRUCT JACOBIAN AND ITS INVERSE
      DO 130 I = 1,NDM
      DO 130 J = 1,2
      XS(I,J) = 0.D0
      DO 130 K = 1,NEL
130   XS(I,J) = XS(I,J) + X(I,K)*SHP(J,K)
      XSJ = XS(1,1)*XS(2,2)-XS(1,2)*XS(2,1)
c     IF(XSJ.LE.0.D0) THEN
c     WRITE(6,2001) NELMT,XSJ
c     STOP
c     END IF
      IF (.NOT.FLAG) RETURN
      SX(1,1) = XS(2,2)/XSJ
      SX(2,2) = XS(1,1)/XSJ
      SX(1,2) =-XS(1,2)/XSJ
      SX(2,1) =-XS(2,1)/XSJ
C---- FORM GLOBAL DERIVATIVES
      DO 140 I = 1,NEL
      TP        = SHP(1,I)*SX(1,1)+SHP(2,I)*SX(2,1)
      SHP(2,I)  = SHP(1,I)*SX(1,2)+SHP(2,I)*SX(2,2)
140   SHP(1,I) = TP
      RETURN
2001  FORMAT(//5X,'** FATAL ERROR ** JACOBIAN OF TRANSFORMATION',
     1 ' FOR ELEM',I5,' IS ',E15.7)
      END
