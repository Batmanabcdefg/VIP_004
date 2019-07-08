	Subroutine input(d)
	implicit real*8(a-h,o-z)
	dimension d(*)
c-----  d(1) : Young's modulus E   -->  Bulk modulus kappa
c-----  d(2) : Poission's ratio nu ---> shear modulus mu
c-----  d(3) : Coeff thermal expansion : alfa_th
c-----  d(4) : Initial temp theta_o
c-----  d(5) : Frequency of atomic vibration f_o
c-----  d(6) : Constant k in cohesion evolution law
c-----  d(7) : Reference shear strain rate gamma^dot_o
c-----  d(8) : Rate exponent a
c-----  d(9) : Shape factor phi
c-----  d(10) : Initial value xi_o of free volume
c-----  d(11) : Free volume creation coeff xeta
c-----  d(12) : Reference temperatute theta^*
c-----  d(13) : Free volume concn at temp theta^* : xi^*
c-----  d(14) : Coeff for evolution of xi_T :  k_theta
c-----  d(15) : constant s_1 in evolution of free volume xi
c-----  d(16) : constant s_2 in evolution of free volume xi
c-----  d(17) : constant s_3 in evolution of free volume xi
c
        read(5,1000) (d(i),i=1,17)
        write(6,2000) (d(i),i=1,17)
        e = d(1)
        pois = d(2)
c---- Store shear modulus & bulk modulus in D(1) & D(2)
c     D(1) = E*D(2)/((1.d0-2.d0*D(2))*(1.d0+D(2)))
c     d(1) = e/(3.d0*(1.d0 - 2.d0*pois))
c     D(2) = 0.5*E/(1.d0 + pois)
	return
1000  format(8f10.0)
 2000 FORMAT(/5X,'EKAMBARAM-THAMBURAJA  BMG MODEL',/
     1 /10X,'YOUNG*S MODULUS -> Bulk Modulus  ',E10.2
     2 /10X,'POISSON*S RATIO -> Shear Modulus ',E10.2
     3 /10X,'COEFF OF THERMAL EXP             ',E10.2
     4 /10X,'INITIAL TEMP THETA_O             ',E10.2
     5 /10X,'FREQUENCY OF ATOMIC VIBRN f_O    ',E10.2
     6 /10X,'CONSTANT k in FREE VOL EVOL LAW  ',E10.2
     7 /10X,'REF SHEAR STRAIN RATE            ',E10.2
     8 /10X,'RATE EXP a                       ',E10.2
     9 /10X,'SHAPE FACTOR PHI                 ',E10.2
     9 /10X,'INITIAL VALUE xi_o OF FREE VOL   ',E10.2
     1 /10X,'FREE VOL CREATION COEFF xeta     ',E10.2
     2 /10X,'REFERENCE TEMP THETA^*           ',E10.2
     3 /10X,'FREE VOL CONCN AT THETA^*:  xi^* ',E10.2
     4 /10X,'COEFF FOR EVOLN OF xi_t: k_theta ',E10.2
     5 /10X,'CONST IN EVOLN OF FREE VOL: S1   ',E10.2
     6 /10X,'CONST IN EVOLN OF FREE VOL: S2   ',E10.2
     7 /10X,'CONST IN EVOLN OF FREE VOL: S3   ',E10.2)
	end

C	================================================================	
 	subroutine update(d,fn1,sig,eps,ql,del2xi,nq,ns)
C	SUBROUTINE TO UPDATE Ekambaram-Thamburaja constitutive model for BMGs USING EXPLICIT METHOD
c
c------ Fn1 : Deformation gradient tensor at time step (n+1)
c------ eps <--> se : 2nd P-K stress corresp to intermediate confign at t_n (updated to t_n+1) 
c------ sig : Cauchy stress at t_n (updated by this sub to t_n+1)  
c------ Fn^p : Plastic part of deformation gradient 
c------ ql : Internal variables at t_n (updated by this sub to t_n+1)
c------
c------ ql(1)  :  xi^dot_o
c------ ql(2)  :  gamma^dot
c------ ql(3)  :  tau - equivalent shear stress
c------ ql(4)  :  p^bar : Hydrostatic pressure 
c------ ql(5) :   temperature at t_n : theta_n
c------ ql(6) :  temperature at t_n+1 : theta_n+1
c------ ql(7) :  free volume xi at t_n : xi_n
c------ ql(8) :  free volume xi at t_n+1 : xi_n+1
c------ ql(9) : Trace(E^e)
c------ ql(10) : Trace(E^dot^e) or Delta Trace(E^e)
c------ ql(11) - ql(19): Plastic part of deformation gradient Fn^p
c------ ql(20) : xi_T at t_n+1
c------ ql(21) : Initial value of cohesion c_o
c------ ql(22) : Current value of cohesion c
c------ ql(23) : Maximum principal log plastic strain
c------ ql(24) : Minimum principal log plastic strain
c
c------ d(1) - d(17) : material parameters as above
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	common /tdata/ time,dtime,bet1,gam1
        dimension d(*),fn1(3,*),sig(*),eps(*),ql(*),voigt1(6)
	DIMENSION Fn_p(3,3),se(3,3),temp(3,3),temp1(3,3),voigt2(9)
        dimension fn1_e(3,3),fn1_p(3,3),fn1_p_inv(3,3),sigc(3,3)
C
C	INITIALIZE PARAMETERS
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
      xi_dot_on = ql(1)
      gama_dot_n = ql(2)
      tau_n = ql(3)
      p_n = ql(4)
      theta_n = ql(5)
      del_theta_n = ql(6) - ql(5)
      theta_n1 = ql(6)
      xi_n = ql(7)
      del_xi_n = ql(8) - ql(7)
      xi_n1 = ql(8)
      tree = ql(9)
      del_trE_e = ql(10)
      xi_T = ql(20)
      coh_o = ql(21)
      coh = ql(22)
c     tau_y = ek_b*theta_n1/omega
c     tau_y = 2.d+08
c
c      if(time.le.dtime) then 
       if(tau_n.le.(0.5*coh_o)) then 
c---- Alternately elastic update can be performed if tau_n < Initial value of cohesion
c---- perform purely elastic update
c
c---- retrieve Fn^p
      call invoigtv(fn1_p,voigt1,ql(11),2)
c
c---- Invert plastic deformation gradient
	call inverse(Fn1_p,Fn1_p_inv,det_fn1_p)
c----  form fn1_e
	do i=1,3
	  do j=1,3
	    fn1_e(i,j)=0.d0
		do k=1,3
		  fn1_e(i,j)=fn1_e(i,j)+fn1(i,k)*fn1_p_inv(k,j)
		enddo
	  enddo
	enddo
        call elas_update(d,fn1_e,se,p_n1,tau_n1,tree,theta_n1)
c
      xi_dot_on1 = 0.d0 
      gama_dot_n1 = 0.d0
      xi_T = xi_r + ek_theta*(theta_n1 - theta_r)
c 
       else
c
c---- perform elastic-viscoplastic update
c---- retrieve Fn^p
      call invoigtv(fn_p,voigt1,ql(11),2)
c
c---- Retrieve 2nd P-K stress based on intermediate config
      call invoigtv(se,eps,voigt2,1)
c
c----	Step 1: form flow tensor N
	do i=1,3
	  do j=i,3
		if(j.eq.i)then
                temp(i,j)=(se(i,j)+p_n)/(dsqrt(2.d0)*tau_n)
		else
                temp(i,j)=se(i,j)/(dsqrt(2.d0)*tau_n)
		endif
	  enddo
	enddo
c
c---- Step 2: form Delta D^p
	fnpcoeff1=dtime*gama_dot_n/dsqrt(2.d0)
        fnpcoeff2=del_xi_n/3.d0
c       fnpcoeff2=0.d0
        do i=1,3
         do j=i,3
          if(j.eq.i) then
               temp(i,j) = temp(i,j)*fnpcoeff1 + fnpcoeff2
           else
               temp(i,j) = temp(i,j)*fnpcoeff1
           end if
          end do
        end do
        temp(2,1) = temp(1,2)
        temp(3,2) = temp(2,3)
        temp(3,1) = temp(1,3)
c
c---- Step 3: update F^p
	call mat_mult(temp,Fn_p,temp1,3)
	do i=1,3
	  do j=1,3
		Fn1_p(i,j)=Fn_p(i,j) + temp1(i,j)
	  enddo
	enddo

c---- Step 4(a): invert plastic deformation gradient
	call inverse(Fn1_p,Fn1_p_inv,det_fn1_p)
c---- Step 4(b): form fn1_e
	do i=1,3
	  do j=1,3
	    fn1_e(i,j)=0.d0
		do k=1,3
		  fn1_e(i,j)=fn1_e(i,j)+fn1(i,k)*fn1_p_inv(k,j)
		enddo
	  enddo
	enddo
c
c---- Step 5 to Step 7 : Update  2nd PK stress tensor se, tau, pe, tree
        call elas_update(d,fn1_e,se,p_n1,tau_n1,tree,theta_n1)
c
c---- Store F(n+1)^p in  ql
	  call voigtv(fn1_p,voigt1,ql(11),2)
c
c     end if
c
c---- Step 8: Compute xi_dot_(n+1), gama_dot_o(n+1)
      call gdot(d,tau_n1,p_n1,theta_n1,xi_n1,xi_T,xi_dot_on1,
     1 gama_dot_n1,coh_o,coh,del2xi,tau_v1,tau_v2)
c
      end if
c
C----- update internal variables
      ql(1) = xi_dot_on1 
      ql(2) = gama_dot_n1 
      ql(3) = tau_n1
      ql(4) = p_n1
      ql(10) = tree - ql(9)
      ql(9) = tree
      ql(20) = xi_T
      ql(22) = coh
c
      ql(25) = tau_v1
      ql(26) = tau_v2
c---- update 2nd P-K stress (eps) & Cauchy stress (sig)
      call cauchy (se,sigc,fn1_e)
      call voigtv(se,eps,voigt2,1)
      call voigtv(sigc,sig,voigt2,1)
c---- find maximum and minimum principal log plastic strains
c     call plogpstrain(fn1_p,peplog1,peplog2)
c     ql(23) = peplog1
c     ql(24) = peplog2
c
	return
	end

c	****************************************************************
        subroutine elas_update(d,fe,se,pe,tau,tree,theta)
        implicit double precision (a-h,o-z)
        dimension d(*),se(3,3),ce(3,3),ee(3,3),
     1 fe(3,3),soe(3,3),eoe(3,3)
      ekappa = d(1)
      emu = d(2)
      alfa_th = d(3)
      theta_o = d(4)      
c
c---- Form ce and ee
      do i=1,3
       do j=i,3
        ce(i,j) = 0.d0
         do l=1,3
          ce(i,j) = ce(i,j) + fe(l,i)*fe(l,j)
         end do
        end do
       end do
c
      do i = 1,3
         do j=i,3
           if(i.eq.j) then
               ee(i,j) = (ce(i,j) - 1.d0)/2.d0
                else
                  ee(i,j) = ce(i,j)/2.d0
            end if
         end do
      end do
c
c---- get trace(Ee), hydrostatic & deviatoric parts of Se
        tree = (ee(1,1) + ee(2,2) + ee(3,3))/3.d0
         do i=1,3
          do j=i,3
          if(i.eq.j) then
           eoe(i,j) = ee(i,j) - tree
            else
             eoe(i,j) = ee(i,j)
            end if
          end do
         end do
          tree = 3.d0*tree
         pe = -ekappa*(tree - 3.d0*alfa_th*(theta - theta_o))
         do i=1,3
          do j=i,3
            soe(i,j) = 2.d0*emu*eoe(i,j)
          end do
         end do
       do i=1,3
       do j=i,3
        ce(j,i) = ce(i,j)
        ee(j,i) = ee(i,j)
        soe(j,i) = soe(i,j)
       end do
      end do
      do i=1,3
       do j=1,3
        if(i.eq.j) then
           se(i,j) = soe(i,j) - pe
        else
           se(i,j) = soe(i,j)
        end if
       end do
      end do
c---- calculate equivalent stress tau
      tau = soe(1,1)**2 + soe(2,2)**2 + soe(3,3)**2 + 2.d0*(
     1 soe(1,2)**2 + soe(2,3)**2 + soe(1,3)**2)
      tau = dsqrt(tau/2.d0)
c
      return
      end
c-------------------------------------------------------------------------
      subroutine gdot(d,tau,p,theta,xi,xi_T,xi_dot_on1, gama_dot_n1,
     1 coh_o,coh,del2xi,tau_v1,tau_v2)
      implicit double precision (a-h,o-z)
      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL
      common /tdata/ time,dtime,bet1,gam1
      dimension d(*)
      data tol/1.d-06/,uplimit/12.d0/
c
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
      if(dabs(xi).gt.tol) then
         xi_dot_on1 = f_o*dexp(-phi/xi)
      else
         xi_dot_on1 = 0.d0
      end if
         xi_T = xi_r + ek_theta*(theta - theta_r)
         tau_v1 = s2*(xi - xi_T)
         tau_v2 =  - s1*del2xi
c        tau_v = tau_v1 + tau_v2
         tau_v = tau_v1 
c---- Note we have ignored -s1*del^2(xi) contribution to viscous stress tau_v 
         f1 = tau + xeta*(-tau_v - p)
c
c        arg = tau/2.d+08
         if(f1.gt.0.d0) then
         coh = coh_o*dexp(ek_c*(xi - xi_o))
         arg = f1/coh
c     if(dabs(arg).gt.uplimit) then
c        write(6,2001) nelmt, time, arg 
c        arg = uplimit
c        end if
c
         arg1 = dlog(arg)/a_rexp
         gama_dot_n1 = gama_dot_o*dexp(arg1)
c        arg1 = 10.d0*dlog(dabs(arg))
c        gama_dot_n1 = gama_dot_on1*dexp(arg1)
       else
c        xi_dot_on1 = 0.d0
         gama_dot_n1 = 0.d0
      end if
      return
2001  format(5x,'Arg of sinh exceeds 12 in elem=',i5,' at time=',
     1 e13.5, ' Arg=',e13.5)
      end
c----------------------------------------------------------------------
      subroutine cauchy(se,sigc,fe)
      implicit double precision (a-h,o-z)
      dimension se(3,3),sigc(3,3),fe(3,3),fet(3,3),temp1(3,3)
c
c---- get transpose of fe
      do i=1,3
       do j=1,3
        fet(i,j) = fe(j,i)
       end do
      end do
c
c---- get det(Fe) = J^e
	xje=fe(1,1)*(fe(2,2)*fe(3,3)-fe(3,2)*fe(2,3))+fe(1,2)*(fe(3,1)
     1 *fe(2,3)-fe(2,1)*fe(3,3))+fe(1,3)*(fe(2,1)*fe(3,2)
     2 -fe(3,1)*fe(2,2))
c
c---- form se*fe^T
      do i=1,3
       do j=1,3
        temp1(i,j) = 0.d0
         do k=1,3
          temp1(i,j) = temp1(i,j) + se(i,k)*fet(k,j)
         end do
        end do
       end do
c
c---- form fe*se*fe^T
      do i=1,3
       do j=i,3
        sigc(i,j) = 0.d0
         do k=1,3
          sigc(i,j) = sigc(i,j) + fe(i,k)*temp1(k,j)
         end do
         sigc(i,j) = sigc(i,j)/xje
        end do
       end do
       sigc(2,1) = sigc(1,2)
       sigc(3,1) = sigc(1,3)
       sigc(3,2) = sigc(2,3)
      return
      end
c-------------------------------------------------------------------
c
      subroutine elas(d,dt)
      implicit double precision (a-h,o-z)
c
c---- isotropic elasticity matrix
c     d(1) : 1st Lame const lamda
c     d(2) : 2nd Lame const mu
c
      dimension d(*),dt(6,6)
      do i = 1,6
         do j = 1,6
            dt(i,j) = 0.d0
         end do
      end do
      do i = 1,3
         do j = 1,3
            dt(i,j) = d(1)
         end do
         dt(i,i) = dt(i,i) + 2.d0*d(2)
      end do
      do i = 4,6
         dt(i,i) = d(2)
      end do
      return
      end
c------------------------------------------------------------------------
      subroutine voigtv(a,voigt1,voigt2,isym)
      implicit double precision (a-h,o-z)
c
c---- converts a 2nd order tensor to column vector using chi2 mapping
c---- Note order : {11,22,33,12,13,23,21,31,32}
c
      dimension a(3,*), voigt1(*),voigt2(*)
      if(isym.eq.1) then
       do i=1,3
        voigt1(i) = a(i,i)
       end do
        voigt1(4) = a(1,2)
        voigt1(5) = a(1,3)
        voigt1(6) = a(2,3)
      else
       do i=1,3
        voigt2(i) = a(i,i)
       end do
        voigt2(4) = a(1,2)
        voigt2(5) = a(1,3)
        voigt2(6) = a(2,3)
        voigt2(7) = a(2,1)
        voigt2(8) = a(3,1)
        voigt2(9) = a(3,2)
      end if
      return
      end       
c
c--------------------------------------------------------------------------------------
      subroutine invoigtv(a,voigt1,voigt2,isym)
      implicit double precision (a-h,o-z)
c
c---- obtains a 2nd order tensor from column vector using inverse chi2 mapping
c---- Note order : {11,22,33,12,13,23,21,31,32}
c
      dimension a(3,*), voigt1(*),voigt2(*)
      if(isym.eq.1) then
       do i=1,3
         a(i,i) = voigt1(i) 
       end do
         a(1,2) = voigt1(4) 
         a(1,3) = voigt1(5) 
         a(2,3) = voigt1(6) 
         a(2,1) = a(1,2)
         a(3,1) = a(1,3)
         a(3,2) = a(2,3)
      else
       do i=1,3
         a(i,i) = voigt2(i) 
       end do
         a(1,2) = voigt2(4) 
         a(1,3) = voigt2(5) 
         a(2,3) = voigt2(6)
         a(2,1) = voigt2(7) 
         a(3,1) = voigt2(8) 
         a(3,2) = voigt2(9)
      end if
      return
      end
c
C	================================================================
	Subroutine inverse(A,A_inv,det_A)
	implicit double precision (a-h,o-z)
	Dimension A(3,*),A_inv(3,*)
c
	det_A=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))+A(1,2)*(A(3,1)*A(2,3)
     1-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
c
	A_inv(1,1)=(A(2,2)*A(3,3)-A(3,2)*A(2,3))/det_A
     	A_inv(2,1)=(A(3,1)*A(2,3)-A(2,1)*A(3,3))/det_A
     	A_inv(3,1)=(A(2,1)*A(3,2)-A(3,1)*A(2,2))/det_A
     	A_inv(1,2)=(A(1,3)*A(3,2)-A(1,2)*A(3,3))/det_A
     	A_inv(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))/det_A
     	A_inv(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/det_A
     	A_inv(1,3)=(A(1,2)*A(2,3)-A(2,2)*A(1,3))/det_A
     	A_inv(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))/det_A
     	A_inv(3,3)=(A(1,1)*A(2,2)-A(2,1)*A(1,2))/det_A
c
	Return
	End
C	*************************************************************
	SUBROUTINE MAT_MULT (A,B,C,SIZE)
	
	implicit double precision (a-h,o-z)
	DIMENSION A(3,*),B(3,*),C(3,*)
	INTEGER SIZE
c
	TEMP=0.0d0
	DO 11 I=1,SIZE
	  DO 11 J=1,SIZE
	    TEMP=0.0d0
		DO K=1,SIZE
		  TEMP=TEMP+A(I,K)*B(K,J)
		ENDDO
		C(I,J)=TEMP
11	CONTINUE
c
	RETURN
	END
c
C	*************************************************************
	SUBROUTINE MAT_TRANSPOSE (A,B,SIZE)
c
	implicit double precision (a-h,o-z)
	DIMENSION A(3,*),B(3,*)
	INTEGER SIZE
c
	DO I=1,SIZE
	  DO J=1,SIZE
		B(I,J)=A(J,I)
	  ENDDO
	ENDDO
c
	RETURN
	END
c------------------------------------------------------------------
      subroutine plogpstrain(fp,peplog1,peplog2)
      implicit double precision (a-h,o-z)
c
      logical flag
      common /eldata/DM,NELMT,MA,MCT,IEL,NEL 
      dimension  fp(3,3),cp(3,3)
      data pi/3.141592653/
c---- Compute maximum and minimum principal log plastic strains
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
      peplog2 = dlog(pstr2)
       return
2040  format(5x,'Negative e-value for CP in ivarout for elem',i6)
       end
c------------------------------------------------------------