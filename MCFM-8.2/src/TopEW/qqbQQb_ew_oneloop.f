      subroutine qqbQQb_ew_oneloop(corr,s,beta,z)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'masses.f'
      include 'epinv.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'anomcoup.f'
      real(dp):: corr(nf),s,t,u,beta,z,alpha,sigma0,gvt,gat,
     .     xI1,xI2,xI3,xI4,D6,sw2,cw2,mw,mz,mh,born,f1,f2,
     .     rz,rw,rb,rh,yphi,ys,ini_corr(nf),dFD(nf),T3(nf),gvq(nf),
     .     gaq(nf),qa(5),bxew(nf,-1:0),bxqcd(nf,-1:0),BB(nf,-1:0),fac
      real(dp) :: db0,gvt_sq,gat_sq,gw_sq,g_rest,MARKUS_qa(5),myDB0
      real(dp) :: C33phiq3,C33phiu,vol2,vol4,kap,kapT
      integer ep
      complex(dp) :: test
      logical,save :: FirstTime=.true.
      

!          s=   173308.20441330347d0
!          beta = 0.55464632643236811d0
!          z = -0.61536288270299611d0
         mb = 1d-16



      ep = 0
c      musq = (2._dp*mt)**2
      alpha = esq/fourpi
c      alpha = 1._dp/126.3_dp
!      as = gsq/fourpi

      mw = wmass
      mz = zmass
      mh = hmass



      sw2 = xw
! on-shell condition cross-check w/ Doreen
!      sw2 = 1._dp - mw**2/mz**2
      cw2 = 1._dp - sw2

      T3 = (/-0.5_dp,0.5_dp,-0.5_dp,0.5_dp,-0.5_dp/)

      gvq(1:5) = (T3(1:5)-2._dp*sw2*Q(1:5))/2._dp/sqrt(sw2*cw2)
      gaq(1:5) = T3(1:5)/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)

      gvt = gvq(2)
      gat = gaq(2)


c**********************************************************************************
c     MARKUS: add dim-6 operator contributions ( variables are in common block of anomcoup.f and set in reader_input.f )       


        C33phiq3 = C_phiq_333 
        C33phiu = C_phiu_33
        kap = kappa_htt
        kapT= kappa_tilde_htt
        
        call ResetEWCouplings(sw2,gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4)!   defined in qqb_QQb_mix.f
        

      
      
        if( FirstTime ) then
          print *, "SM  top couplings:"
          print *, "gvt   =",(T3(2)-2._dp*sw2*Q(2))/2._dp/sqrt(sw2*cw2)
          print *, "gat   =",T3(2)/2._dp/sqrt(sw2*cw2)
          print *, "gw    =",0.5_dp/sqrt(2._dp*sw2)
          print *, ""
          print *, "New top couplings:"
          print *, "gvt   =",gvt
          print *, "gat   =",gat
          print *, "gw    =",gw
          print *, "gvt_sq=",gvt_sq
          print *, "gat_sq=",gat_sq
          print *, "gw_sq =",gw_sq
          print *, "g_rest=",g_rest
          print *, ""
          print *, "mtop, alphas, alpha, mu",mt, as, alpha,dsqrt(musq)
          FirstTime=.false.
        endif
        
                
c     END MARKUS      
c**********************************************************************************

    
!      as = 1._dp
      sigma0 = pi*as**2/8._dp*(Nc**2-1._dp)/Nc**2*beta/s  !  MARKUS: Note that the tree is proportional to it and all loop corrections are proportional to it. Therefore the return value vew/born is independent of it. 
                                        !          Hence, the Jacobian from dz is also dropping out. And the Jacobian is whatever is in the born of m2 outside of this routine. 

      t = -s*(1._dp-beta*z)/2._dp+mt**2
      u = -s*(1._dp+beta*z)/2._dp+mt**2

      
      
      rz = mz**2/s
      rb = mb**2/s
      rw = mw**2/s
      rh = mh**2/s

      yphi = mb**2/mw**2
      ys = 1._dp - beta**2 + 4._dp*rb

C --- set BB  = 0 for now
C      BB = 0._dp
      fac = pi*alpha*as*(Nc**2-1._dp)*8._dp*s/(s - mz**2)
      BB(:,0) = - fac*(
     .     + (2._dp - beta**2*(1._dp - z**2))*gvq(:)*gvt
     .     + 2._dp*beta*z*gaq(:)*gat )

!      write(6,*) 'BB: ', BB(1:2,0)
!      write(6,*) 'born: ', 4._dp*s/(s - mz**2)*(
!     .     + (2._dp - beta**2*(1._dp - z**2))*gvq(1:2)*gvt
!     .     - 2._dp*beta*z*gaq(1:2)*gat )
!      stop
!      BB(:,-1) = - 2._dp*fac*(-gvq(:)*gvt + 3._dp*beta*z*gaq(:)*gat)


!      MARKUS: this is the weak vertex correction at the initial state
!       dFD = - alpha/8._dp/pi*((gvq**2+gaq**2)*f1(rz)*g_rest+2._dp*gw_sq*f1(rw))
      dFD = - alpha/8._dp/pi*( (gvq**2+gaq**2)*f1(rz) * g_rest
     .      +2._dp*(0.5_dp/sqrt(2._dp*sw2))**2*f1(rw) * g_rest ) 

      born = sigma0*(2._dp - beta**2 + beta**2*z**2)

      ini_corr = 2._dp*born*dFD

      qa(1) = 
     .     (0.25_dp*alpha*sigma0*(-2._dp*(gat_sq + gvt_sq)*(1._dp + z**2) 
     .     + 2._dp*((1._dp - beta**2)*(-3._dp*gat_sq + gvt_sq) 
     .     + 2._dp*(gat_sq + gvt_sq)*rz)*s*(2._dp - beta**2*(1._dp 
     .     - z**2))*db0(mt**2,mz**2,mt**2) + (8._dp*(gat_sq 
     .     + gvt_sq)*(1._dp + z**2)*(-xI1(mt**2,musq,ep) 
     .     + xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s) 
     .     - 4._dp*(gat_sq*(1._dp - 5._dp*z**2 - 3._dp*beta**2*(1._dp 
     .     - z**2)) - gvt_sq*(3._dp + z**2 - beta**2*(1._dp - z**2)) 
     .     + ((gat_sq + gvt_sq)*rz*(1._dp - beta**2 - 3._dp*z**2 
     .     + 7._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)))/(beta**2*
     .     (1._dp - beta**2)))*xI2(mt**2,mz**2,mt**2,musq,ep) 
     .     + 2._dp*(-4._dp*gvt_sq*(1._dp + z**2) - (-3._dp*gat_sq 
     .     + gvt_sq)*f2(z,beta) + (2._dp*(gat_sq + gvt_sq)*rz*f2(z,beta))
     .     /beta**2)*xI2(s,mt**2,mt**2,musq,ep) + 2._dp*s*(-(1._dp 
     .     + beta**2)*gvt_sq*(2._dp - beta**2*(1._dp - z**2)) 
     .     - gat_sq*(2._dp - 3._dp*beta**2 + 5._dp*beta**2*z**2 
     .     + 3._dp*beta**4*(1._dp - z**2)) + (2._dp*(gat_sq 
     .     + gvt_sq)*rz**2*f2(z,beta))/beta**2 + 4._dp*rz*(-gvt_sq*
     .     (1._dp + z**2) + gat_sq*f2(z,beta)))*xI3(mt**2,mt**2,s,mt**2,
     .     mz**2,mt**2,musq,ep)))/pi

      qa(2) = 
     .     (0.5_dp*alpha*gw_sq*sigma0*(-2._dp*(1._dp + z**2) 
     .     + (-1._dp + beta**2 + 4._dp*(-rb + rw))*s*(2._dp - beta**2*(1._dp 
     .     - z**2))*db0(mt**2,mb**2,mw**2) + (8._dp*(1._dp 
     .     + z**2)*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))/((1._dp 
     .     - beta**2)*s) - (1._dp*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp 
     .     - z**2) - 5._dp*beta**2*(1._dp + z**2) + (4._dp*(-rb + rw)*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2)))/(1._dp - beta**2))*xI2(mt**2,mb**2,mw**2,musq,
     .     ep))/beta**2 + ((1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp 
     .     - z**2) - 5._dp*beta**2*(1._dp + z**2) + 4._dp*(-rb + rw)*
     .     f2(z,beta))*xI2(s,mb**2,mb**2,musq,ep))/beta**2 + (0.25_dp*s*
     .     (1._dp - 16._dp*beta**2 + beta**4 - 3._dp*z**2 - 4._dp*beta**2*
     .     z**2 - 11._dp*beta**4*z**2 - 2._dp*beta**6*(1._dp - z**2) 
     .     + 8._dp*rw*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp - z**2) 
     .     - 3._dp*beta**2*(1._dp + z**2)) + 8._dp*(-(1._dp + beta**2)*rb 
     .     + 2._dp*(-rb + rw)**2)*f2(z,beta))*xI3(mt**2,mt**2,s,mb**2,
     .     mw**2,mb**2,musq,ep))/beta**2))/pi

     
     
!      MARKUS: chi_0 boson exchange
!      MARKUS: replaced gat_sq coupl. by explicit expression  1d0/16._dp/sw2/cw2
      qa(3) = (0.5_dp*alpha/16._dp/sw2/cw2*mt**2*sigma0*(
     .     (-2._dp*(1._dp+z**2))
     .     +4._dp*rz*s*(2._dp - beta**2*(1._dp - z**2))*db0(mt**2,mz**2,
     .     mt**2)
     .     - (8._dp*(1._dp + z**2)*(xI1(mt**2,musq,ep)
     .     - xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s) 
     .     - (4._dp*rz*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2))*xI2(mt**2,mz**2,mt**2,musq,ep))/(beta**2*(1._dp 
     .     - beta**2)) 
     .     + (2._dp*(beta**2*(1._dp + z**2) + 2._dp*rz*
     .     f2(z,beta))*xI2(s,mt**2,mt**2,musq,ep))/beta**2 
     .     + (4._dp*s*(beta**2*(1._dp - beta**2)*rz*(1._dp - z**2) 
     .     + rz**2*f2(z,beta))*xI3(mt**2,mt**2,s,mt**2,mz**2,mt**2,musq,
     .     ep))/beta**2)
     .     )/(mz**2*pi)

                
   
!      MARKUS: phi^+- boson exchange
!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
      qa(4)=(0.5_dp*alpha*0.125_dp/sw2*sigma0*( (-0.25_dp*ys*(1._dp+z**2))/rw
     .     + 0.125_dp*s*((-(1._dp - beta**2)**2)/rw + 8._dp*(1._dp 
     .     - beta**2)*yphi - 16._dp*rb*yphi + 4._dp*ys)*(2._dp - beta**2*
     .     (1._dp - z**2))*db0(mt**2,mb**2,mw**2) + (ys*(1._dp 
     .     + z**2)*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))/((1._dp 
     .     - beta**2)*mw**2) - (0.125_dp*(-16._dp*beta**2*
     .     (1._dp - beta**2)**2*yphi*(1._dp - z**2) - 16._dp*rb*yphi*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2)) + 4._dp*ys*(1._dp - beta**2 - 3._dp*z**2 
     .     + 7._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)) 
     .     + ((1._dp - beta**2)**2*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp 
     .     - z**2) + 3._dp*beta**2*(1._dp + z**2)))/rw)*xI2(mt**2,mb**2,
     .     mw**2,musq,ep))/(beta**2*(1._dp - beta**2)) + (0.125_dp*(((1._dp 
     .     - beta**2)*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp - z**2) 
     .     + 3._dp*beta**2*(1._dp + z**2)))/rw + 4._dp*(-2._dp*beta**2*yphi 
     .     - 4._dp*rb*yphi + ys)*f2(z,beta))*xI2(s,mb**2,mb**2,musq,ep))/
     .     beta**2 + (0.03125_dp*s*(8._dp*(1._dp - beta**2)**2*(1._dp 
     .     - 3._dp*z**2 + 2._dp*beta**2*(1._dp - z**2)) - 4._dp*(1._dp 
     .     - beta**2)*yphi*(1._dp - 11._dp*beta**2 - 3._dp*z**2 
     .     - 3._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)) 
     .     + ((1._dp - beta**2)**2*(1._dp - 7._dp*beta**2 - 3._dp*z**2 
     .     + beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)))/rw 
     .     - 16._dp*rb*yphi*(1._dp + beta**2 - 3._dp*z**2 
     .     + 9._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)) 
     .     + 16._dp*(-8._dp*rb**2 + 4._dp*rb**2*yphi + rw*ys)*f2(z,beta))*
     .     xI3(mt**2,mt**2,s,mb**2,mw**2,mb**2,musq,ep))/beta**2))/pi

     
           
!      MARKUS: Higgs boson exchange
!      MARKUS: replaced gw_sq coupl. by vev
      qa(5) = 
     .     (0.5_dp/vev**2*mt**2*sigma0*(-1._dp - z**2 + 2._dp*(-1._dp 
     .     + beta**2 + rh)*s*(2._dp - beta**2*(1._dp - z**2))*db0(mt**2,
     .     mt**2,mh**2) + (4._dp*(1._dp + z**2)*
     .     (xI1(mh**2,musq,ep) - xI1(mt**2,musq,ep)))/((1._dp - beta**2)*
     .     s) - 2._dp*(2._dp*(1._dp - beta**2)*(1._dp - z**2) + (rh*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2)))/(beta**2*(1._dp - beta**2)))*xI2(mt**2,mt**2,
     .     mh**2,musq,ep) + ((beta**2*(5._dp - 3._dp*z**2 - 4._dp*beta**2*
     .     (1._dp - z**2)) + 2._dp*rh*f2(z,beta))*xI2(s,mt**2,mt**2,musq,
     .     ep))/beta**2 + (2._dp*s*(3._dp*beta**2*(1._dp - beta**2)*rh*
     .     (1._dp - z**2) - beta**2*(1._dp - beta**2)*(2._dp - beta**2*
     .     (1._dp - z**2)) + rh**2*f2(z,beta))*xI3(mt**2,mt**2,s,mt**2,
     .     mh**2,mt**2,musq,ep))/beta**2))/(8*pi**2)
    
         

                       
     
     
     
      do ep = -1,0
         bxew(:,ep) = 
     .     (-0.015625_dp*as*BB(:,0)*beta*(-(1._dp + beta*z)*
     .     xI3(0._dp,u,mt**2,0._dp,0._dp,mt**2,musq,ep) 
     .     + (1._dp - beta*z)*
     .     xI3(0._dp,t,mt**2,0._dp,0._dp,mt**2,musq,ep)))
     .     /(Nc**2*pi) + (alpha*sigma0*((-2._dp*beta**2*s*(1._dp - z**2)*
     .     (2._dp*gaq*gat*(-rz**2 + beta*z + rz*(1._dp + beta*z)) 
     .     + gvq*gvt*(3._dp - beta**2 + 2._dp*rz**2 + 2._dp*beta*z*(1._dp 
     .     + beta*z) - rz*(1._dp + beta**2 + 2._dp*beta*z)))*D6(0._dp,0._dp
     .     ,mt**2,mt**2,s,u,0._dp,0._dp,mz**2,mt**2,musq,ep))/
     .     ((1._dp - rz)**2*(1._dp + beta*z)) + (2._dp*beta**2*s*(1._dp 
     .     - z**2)*(2._dp*gaq*gat*(rz**2 + beta*z - rz*(1._dp - beta*z)) 
     .     + gvq*gvt*(3._dp - beta**2 + 2._dp*rz**2 - rz*(1._dp + beta**2 
     .     - 2._dp*beta*z) - 2._dp*beta*z*(1._dp - beta*z)))*D6(0._dp,0._dp,
     .     mt**2,mt**2,s,t,0._dp,0._dp,mz**2,mt**2,musq,ep))/
     .     ((1._dp - rz)**2*(1._dp - beta*z)) - (2._dp*(1._dp - beta**2)*
     .     (2._dp*beta*gaq*gat + gvq*gvt*z*(1._dp + 2._dp*beta**2 
     .     - beta**2*z**2))*xI2(mt**2,0._dp,mt**2,musq,ep))/(beta*(1._dp 
     .     - beta**2*z**2)) - (2._dp*(1._dp - beta**2)*(2._dp*beta*gaq*gat 
     .     + gvq*gvt*z*(1._dp + 2._dp*beta**2 - beta**2*z**2))*xI2(mt**2,
     .     mz**2,mt**2,musq,ep))/(beta*(1._dp - beta**2*z**2)) 
     .     + (4._dp*(beta*gaq*gat + gvq*gvt*z)*xI2(s,0._dp,mz**2,musq,ep))
     .     /beta - (2._dp*(-gaq*gat + gvq*gvt)*(1._dp - 2._dp*beta**2 
     .     + beta**2*z**2)*xI2(u,0._dp,mt**2,musq,ep))/(1._dp + beta*z) 
     .     + (2._dp*(gaq*gat + gvq*gvt)*(1._dp - 2._dp*beta**2 
     .     + beta**2*z**2)*xI2(t,0._dp,mt**2,musq,ep))/(1._dp - beta*z) 
     .     - (1._dp*s*(2._dp*gaq*gat*(1._dp - beta**2 + 2._dp*beta*z 
     .     - beta**3*z + 2._dp*beta**2*z**2 + beta**3*z**3 + rz*(-3._dp 
     .     + 3._dp*beta**2 - 3._dp*beta*z + beta**3*z - 2._dp*beta**2*z**2) 
     .     - rz**3*(1._dp - 2._dp*beta**2 + beta**2*z**2) + rz**2*(3._dp 
     .     - 4._dp*beta**2 + beta*z - 2._dp*beta**3*z + beta**2*z**2 
     .     + beta**3*z**3)) + gvq*gvt*(2._dp - beta**2 + 4._dp*beta*z 
     .     - 2._dp*beta**3*z + 3._dp*beta**2*z**2 - beta**4*z**2 
     .     + 2._dp*beta**3*z**3 + beta**4*z**4 + 2._dp*rz**3*(1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2) - rz**2*(4._dp - 5._dp*beta**2 
     .     - beta**4 - 2._dp*beta**3*z + beta**2*z**2 + beta**4*z**2 
     .     + 2._dp*beta**3*z**3) + beta*rz*(-4._dp*beta + beta**3 - 4._dp*z 
     .     - 2._dp*beta**3*z**2 + beta**3*z**4)))*xI3(0._dp,mt**2,u,0._dp,
     .     mz**2,mt**2,musq,ep))/((1._dp - rz)**2*(1._dp + beta*z)) 
     .     + (s*(2._dp*gaq*gat*(-1._dp + beta**2 + 2._dp*beta*z - beta**3*z 
     .     - 2._dp*beta**2*z**2 + beta**3*z**3 + rz**3*(1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2) + rz*(3._dp - 3._dp*beta**2 
     .     - 3._dp*beta*z + beta**3*z + 2._dp*beta**2*z**2) + rz**2*(-3._dp 
     .     + 4._dp*beta**2 + beta*z - 2._dp*beta**3*z - beta**2*z**2 
     .     + beta**3*z**3)) + gvq*gvt*(2._dp - beta**2 - 4._dp*beta*z 
     .     + 2._dp*beta**3*z + 3._dp*beta**2*z**2 - beta**4*z**2 
     .     - 2._dp*beta**3*z**3 + beta**4*z**4 + 2._dp*rz**3*(1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2) - rz**2*(4._dp - 5._dp*beta**2 
     .     - beta**4 + 2._dp*beta**3*z + beta**2*z**2 + beta**4*z**2 
     .     - 2._dp*beta**3*z**3) + beta*rz*(-4._dp*beta + beta**3 + 4._dp*z 
     .     - 2._dp*beta**3*z**2 + beta**3*z**4)))*xI3(0._dp,mt**2,t,0._dp,
     .     mz**2,mt**2,musq,ep))/((1._dp - rz)**2*(1._dp - beta*z)) 
     .     + (2._dp*(1._dp - beta**2)*s*(2._dp*beta*gaq*gat*(rz**2 
     .     + beta**2*z**2 - rz*(1._dp - beta**2*z**2)) + gvq*gvt*z*(1._dp 
     .     + 2._dp*beta**2 - beta**4 - beta**2*z**2 + beta**4*z**2 
     .     - rz**2*(-1._dp - 2._dp*beta**2 + beta**2*z**2) + rz*(-2._dp 
     .     - beta**4 + 2._dp*beta**2*z**2 + beta**4*z**2)))*xI3(s,mt**2,
     .     mt**2,0._dp,mz**2,mt**2,musq,ep))/(beta*(1._dp - rz)*(1._dp 
     .     - beta**2*z**2))))/pi

         bxqcd(:,ep) = 
     .     (-0.015625_dp*as*BB(:,0)*beta*(-(1._dp + beta*z)*
     .     xI3(0._dp,u,mt**2,0._dp,0._dp,mt**2,musq,ep)
     .     + (1._dp - beta*z)*
     .     xI3(0._dp,t,mt**2,0._dp,0._dp,mt**2,musq,ep)))
     .     /(Nc**2*pi) + (alpha*sigma0*((-beta**2*s*(1._dp - z**2)*
     .     (2._dp*beta*gaq*gat*z + gvq*gvt*(3._dp - beta**2 
     .     + 2._dp*beta*z*(1._dp + beta*z)))*D6(0._dp,0._dp,mt**2,mt**2,s,
     .     u,0._dp,0._dp,0._dp,mt**2,musq,ep))/(1._dp + beta*z) 
     .     + (beta**2*s*(1._dp - z**2)*(2._dp*beta*gaq*gat*z 
     .     + gvq*gvt*(3._dp - beta**2 - 2._dp*beta*z*(1._dp 
     .     - beta*z)))*D6(0._dp,0._dp,mt**2,mt**2,s,t,0._dp,0._dp,0._dp,
     .     mt**2,musq,ep))/(1._dp - beta*z) + (2._dp*(1._dp - beta**2)*
     .     (-beta*gaq*gat*(1._dp + beta**2*z**2) + gvq*gvt*z*(-1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2))*
     .     xI2(mt**2,0._dp,mt**2,musq,ep))/(beta*(1._dp - beta**2*z**2)) 
     .     + (2._dp*(beta*gaq*gat + gvq*gvt*z)*xI2(s,0._dp,0._dp,musq,ep))
     .     /beta - (1._dp*(beta*gaq*gat*(beta + z)*(1._dp - beta*z) 
     .     + gvq*gvt*(1._dp - 2._dp*beta**2 + beta**2*z**2))*xI2(u,0._dp,
     .     mt**2,musq,ep))/(1._dp + beta*z) + ((-beta*gaq*gat*(beta - z)*
     .     (1._dp + beta*z) + gvq*gvt*(1._dp - 2._dp*beta**2 
     .     + beta**2*z**2))*xI2(t,0._dp,mt**2,musq,ep))/(1._dp - beta*z) 
     .     + ((1._dp - beta**2)*s*(beta*gaq*gat*(1._dp + beta**2*z**2) 
     .     + gvq*gvt*z*(1._dp + beta**2 + beta**2*(1._dp - beta**2)*(1._dp 
     .     - z**2)))*xI3(s,mt**2,mt**2,0._dp,0._dp,mt**2,musq,ep))/
     .     (beta*(1._dp - beta**2*z**2))))/(pi*(1._dp - rz))


      end do 
      
      ep = 0 
     
!----------------------------------- anomalous Goldstone boson contributions replacing the ones above -------------------------------------------------------------






!      MARKUS: chi0 boson exchange       
         MARKUS_qa(3)= 
     -  (4*pi*alpha)*(4*pi*as)**2*2d0/9d0/(64d0*Pi**2*SW2*MW**2*beta**2)* (     
     -     -2*beta**2*(-1 + beta**2)*MZ**2*s*
     - (1+4*C33phiq3*vol2 + 2*C33phiu*vol2 +(2*C33phiq3 + C33phiu)**2*vol4)*
     -   (2 + beta**2*(-1 + z**2))*DB0(MT**2,MT**2,MZ**2) - 
     - (1+4*C33phiq3*vol2 + 2*C33phiu*vol2 +(2*C33phiq3 + C33phiu)**2*vol4)*
     -   (-(beta**2*(-1 + beta**2)*s*(1 + z**2)) + 
     -     4*beta**2*(1 + z**2)*xI1(MT**2,musq,ep) - 
     -     4*beta**2*(1 + z**2)*xI1(MZ**2,musq,ep) + 
     -     2*MZ**2*xI2(MT**2,MT**2,MZ**2,musq,ep) - 
     -     2*beta**2*MZ**2*xI2(MT**2,MT**2,MZ**2,musq,ep) + 
     -     4*beta**4*MZ**2*xI2(MT**2,MT**2,MZ**2,musq,ep) - 
     -     6*MZ**2*z**2*xI2(MT**2,MT**2,MZ**2,musq,ep) + 
     -     14*beta**2*MZ**2*z**2*xI2(MT**2,MT**2,MZ**2,musq,ep) - 
     -     4*beta**4*MZ**2*z**2*xI2(MT**2,MT**2,MZ**2,musq,ep) - 
     -     2*MZ**2*xI2(s,MT**2,MT**2,musq,ep) + 
     -     6*beta**2*MZ**2*xI2(s,MT**2,MT**2,musq,ep) - 
     -     4*beta**4*MZ**2*xI2(s,MT**2,MT**2,musq,ep) - 
     -     beta**2*s*xI2(s,MT**2,MT**2,musq,ep) + 
     -     beta**4*s*xI2(s,MT**2,MT**2,musq,ep) + 
     -     6*MZ**2*z**2*xI2(s,MT**2,MT**2,musq,ep) - 
     -     10*beta**2*MZ**2*z**2*xI2(s,MT**2,MT**2,musq,ep) + 
     -     4*beta**4*MZ**2*z**2*xI2(s,MT**2,MT**2,musq,ep) - 
     -     beta**2*s*z**2*xI2(s,MT**2,MT**2,musq,ep) + 
     -     beta**4*s*z**2*xI2(s,MT**2,MT**2,musq,ep) + 
     -     2*(-1 + beta**2)*MZ**2*(beta**2*(-1 + beta**2)*s*(-1 + z**2) + 
     -        MZ**2*(1 - 3*z**2 + 2*beta**2*(-1 + z**2)))*
     -      xI3(MT**2,MT**2,s,MT**2,MZ**2,MT**2,musq,ep)) )
     -   *beta/128d0/Pi/s       



       myDB0 = -1d0/MT**2-MW**2/MT**4*dlog(MT**2/MW**2-1d0)! this is   d/dpsq  B0(psq,0,MW**2)  at psq=MT**2, equivalent to db0(mt**2,1e-4,mw**2)

!      MARKUS: phi^+- boson exchange           
         MARKUS_qa(4)=    
     -  (4*pi*alpha)*(4*pi*as)**2*2d0/9d0/(64d0*Pi**2*SW2*MW**2*beta**2)* (          
     -   -(beta**2*(-1 + beta**2)*s*(4*MW**2 + (-1 + beta**2)*s)*
     -      (1 + 2*C33phiq3*vol2 + C33phiq3**2*vol4)*(2 + beta**2*(-1 + z**2))*
     -      myDB0)/2d0 - 
     -  ((1 + 2*C33phiq3*vol2 + C33phiq3**2*vol4)*
     -     (-32*beta**2*(1 + z**2)*xI1(MW**2,musq,ep) - 
     -       4*(4*MW**2*(-1 + beta**2 - 2*beta**4 + 
     -             (3 - 7*beta**2 + 2*beta**4)*z**2) + 
     -          (-1 + beta**2)*s*(1 - 3*z**2 + 2*beta**4*(-1 + z**2) + 
     -             3*beta**2*(1 + z**2)))*xI2(MT**2,0d0,MW**2,musq,ep) + 
     -       (-1 + beta**2)*(-8*beta**2*s*(1 + z**2) + 
     -          4*(4*MW**2*(1 - 3*z**2 + 2*beta**2*(-1 + z**2)) + 
     -         s*(1 - 3*z**2 + 2*beta**4*(-1 + z**2) + 3*beta**2*(1 + z**2)))*
     -           xI2(s,0d0,0d0,musq,ep) + 
     -          ((-1 + beta**2)*s**2*
     -         (-1 + 7*beta**2 - 2*beta**4 + (3 - beta**2 + 2*beta**4)*z**2) + 
     -          16*MW**4*(1 - 3*z**2 + 2*beta**2*(-1 + z**2)) + 
     -        8*(-1 + beta**2)*MW**2*s*(-1 + 3*z**2 + 2*beta**2*(-1 + z**2)))*
     -         xI3(MT**2,MT**2,s,0d0,MW**2,0d0,musq,ep))))/8d0 ) 
     -   *beta/128d0/Pi/s       


!       MARKUS: Higgs boson exchange, so far SM only
       MARKUS_qa(5) = 
     -  (4*pi*alpha)*(4*pi*as)**2*2d0/9d0/(64d0*Pi**2*SW2*MW**2*beta**2)*(
     -   beta**2*(-1 + beta**2)*(kap**2 + kapT**2)*s*(1 + z**2) - 
     -  2*beta**2*(-1 + beta**2)*s*
     -   (kapT**2*MH**2 + kap**2*(MH**2 + (-1 + beta**2)*s))*
     -   (2 + beta**2*(-1 + z**2))*DB0(MT**2,MT**2,MH**2) + 
     -  4*beta**2*(kap**2 + kapT**2)*(1 + z**2)*xI1(MH**2,musq,ep) - 
     -  4*beta**2*(kap**2 + kapT**2)*(1 + z**2)*xI1(MT**2,musq,ep) + 
     -  2*(-1 + beta**2)*(kapT**2*MH**2*(1 - 3*z**2 + 4*beta**2*(-1 + z**2)) + 
     -     kap**2*(2*beta**2*(-1 + beta**2)*s*(-1 + z**2) + 
     -        MH**2*(1 - 3*z**2 + 4*beta**2*(-1 + z**2))))*
     -   xI2(MT**2,MH**2,MT**2,musq,ep) - 
     -  4*beta**2*(kap**2 + kapT**2)*MH**2*(2 + beta**2*(-1 + z**2))*
     -   xI2(MT**2,MT**2,MH**2,musq,ep) - 
     -  (-1 + beta**2)*(kapT**2*(beta**2*s*(1 + z**2) + 
     -        MH**2*(2 - 6*z**2 + 4*beta**2*(-1 + z**2))) + 
     -     kap**2*(MH**2*(2 - 6*z**2 + 4*beta**2*(-1 + z**2)) + 
     -        beta**2*s*(5 - 3*z**2 + 4*beta**2*(-1 + z**2))))*
     -   xI2(s,MT**2,MT**2,musq,ep) - 
     -  2*(-1 + beta**2)*(kapT**2*MH**2*
     -      (beta**2*(-1 + beta**2)*s*(-1 + z**2) + 
     -        MH**2*(1 - 3*z**2 + 2*beta**2*(-1 + z**2))) + 
     -     kap**2*(3*beta**2*(-1 + beta**2)*MH**2*s*(-1 + z**2) + 
     -        beta**2*(-1 + beta**2)*s**2*(2 + beta**2*(-1 + z**2)) + 
     -        MH**4*(1 - 3*z**2 + 2*beta**2*(-1 + z**2))))*
     -   xI3(MT**2,MT**2,s,MT**2,MH**2,MT**2,musq,ep) )
     -   *beta/128d0/Pi/s       


        
!          print *, "old qa3",qa(3)
!          print *, "old qa4",qa(4)
!          print *, "old qa5",qa(5)
!          print *, "new qa5",MARKUS_qa(5)
        
        
        
        
        
        
!       over-writing the exising ones
        qa(3) = MARKUS_qa(3)
        qa(4) = MARKUS_qa(4)
        qa(5) = MARKUS_qa(5)


!         print *, "new qa1",qa(1)
!         print *, "new qa2",qa(2)
!         print *, "new qa3",qa(3)
!         print *, "new qa4",qa(4)
!         print *, "new qa5",qa(5)
!         print *, "new bxew",bxew(:,0)
!         print *, "new bxqcd",bxqcd(:,0)
!         print *, ""
!          print *, "rat qa1",qa(1)/(6.8983663068961290d-011 )    ! comparison with Till, see email from 7/5/19, 4:11 PM
!          print *, "rat qa2",qa(2)/(-6.2956972667896152d-010 )
!          print *, "rat qa3",qa(3)/(-1.6030271106979848d-010 )
!          print *, "rat qa4",qa(4)/(-1.0033401929399480d-009 )
!          print *, "rat bxew",bxew(1,0)/(-9.7712326901855281d-010 )
!          print *, "rat bxqcd",bxqcd(1,0)/(-3.6416923669956070d-011 )
!          print *, "rat bxew",bxew(2,0)/(1.0346766433422984d-009)
!          print *, "rat bxqcd",bxqcd(2,0)/(6.0573496881493194d-011)
!         pause

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------


! check the coefficient prop to epinv
!      write(6,*) bxqcd(1:2,-1)*BB(1:2,0)/BB(1:2,-1)
!      write(6,*) as/32._dp/pi/Nc**2*beta/s*BB(1:2,0)*log((t-mt**2)
!     .     /(u-mt**2))
! check the finite piece from epinv/epinv
!      write(6,*) bxqcd(1:2,-1)
!      write(6,*) as/32._dp/pi/Nc**2*beta/s*BB(1:2,-1)*log((t-mt**2)
!     .     /(u-mt**2))
!      bxew(:,-1) = bxew(:,-1)*(BB(:,0)/BB(:,-1)*epinv + 1._dp)
!      bxqcd(:,-1) = bxqcd(:,-1)*(BB(:,0)/BB(:,-1)*epinv + 1._dp)

      bxew(:,-1)  = bxew(:,-1)*epinv
      bxqcd(:,-1) = bxqcd(:,-1)*epinv

c--- use the value of "tevscale" (passed via anomcoup.f)
c--- as an anomalous top Yukawa coupling:
c--- g(top Y) = (tevscale) x g(SM, top Y)
c--- it therefore affects Higgs diagrams as the square
      qa(5)=qa(5)*tevscale**2


!     MARKUS: marking contributions that do not scale like gvt,gat,gw
      qa(3:5)     = qa(3:5)     * g_rest

!       print*,"QUESTION: why arent't qqb s-channels multiplied with 1/2?!"

      corr = qa(1) + qa(2) + qa(3) + qa(4) + qa(5) 
     .     + bxew(:,0) + bxew(:,-1) + bxqcd(:,0) + bxqcd(:,-1)
      corr = corr + ini_corr



      corr = corr/born




      end subroutine qqbQQb_ew_oneloop

      
      
      
      
      
      
      
      
      
      
      
      
      
      

      function f1(x)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'cplx.h'
      real(dp):: x,cli2,f1

      f1 = real(1._dp + 2._dp*(
     .     + (1._dp+log(cplx1(x)))*(2._dp*x+3._dp)
     .     - 2._dp*(1._dp+x)**2*(cli2(cplx1(1._dp+1._dp/x))
     .     -pi**2/6._dp)),dp)

      end function f1


      function f2(z,beta)
      implicit none
      include 'types.f'
      real(dp):: z,beta,f2

      f2 = 1._dp - 3._dp*z**2 - 2._dp*beta**2*(1._dp - z**2)

      end function f2


      function D6(p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,
     .     m2sq,m3sq,m4sq,musq,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,
     .     m4sq,musq,xI3,xI4,p1dp2,p1dp3,p2dp3,m1,m2,m3,m4,D6
      integer ep
      external xI3,xI4

      p1dp2 = (s12 - p1sq - p2sq)/2._dp
      p2dp3 = (s23 - p2sq - p3sq)/2._dp
      p1dp3 = (p2sq + p4sq - s12 -s23)/2._dp
      m1 = sqrt(m1sq)
      m2 = sqrt(m2sq)
      m3 = sqrt(m3sq)
      m4 = sqrt(m4sq)

      D6 = 
     .     (0.25_dp*(-((-m1**2 + m2**2)*p1dp3*p2sq + p1dp2*((m1**2 
     .     - m2**2 - p1sq)*p2dp3 + p1dp3*(m2**2 - m3**2 + p2sq)) 
     .     - p1dp2**2*(m3**2 - m4**2 + 2._dp*p2dp3 + p3sq) 
     .     + p1sq*((-m2**2 + m3**2)*p2dp3 + p2sq*(m3**2 - m4**2 + p1dp3 
     .     + p2dp3 + p3sq)))*xI3(p1sq,p2sq,s12,m1**2,m2**2,m3**2,musq,
     .     ep) + (p1dp3**2*(m2**2 - m3**2 + p2sq) - p1dp3*(m1**2 - m2**2 
     .     - p1sq)*(p2dp3 + p2sq) - p1dp2**2*(m3**2 - m4**2 + 2._dp*p2dp3 
     .     + p3sq) + p1sq*(2._dp*p2dp3**2 + (m3**2 - m4**2 + p2dp3)*p2sq 
     .     - p2dp3*(m2**2 - 2._dp*m3**2 + m4**2 - p3sq) 
     .     + (-m2**2 + m3**2)*p3sq) + p1dp2*(p1dp3*(m2**2 - 2._dp*m3**2 
     .     + m4**2 - 2._dp*p2dp3 + p2sq - p3sq) + (m1**2 - m2**2 - p1sq)*
     .     (p2dp3 + p3sq)))*xI3(p1sq,s23,p4sq,m1**2,m2**2,m4**2,musq,ep) 
     .     - ((m1**2 - m2**2 - p1sq)*(p2dp3**2 - p2sq*p3sq) 
     .     + p1dp2*(-2._dp*p2dp3**2 + (m2**2 - m3**2 + p2sq)*p3sq 
     .     - p2dp3*(m3**2 - m4**2 + p3sq)) + p1dp3*((-m2**2 + m3**2)*
     .     p2dp3 + p2sq*(m3**2 - m4**2 + p2dp3 + p3sq)))*xI3(p2sq,p3sq,
     .     s23,m2**2,m3**2,m4**2,musq,ep) + (-m3**2*p1sq*p2dp3 
     .     + m4**2*p1sq*p2dp3 + m1**2*p2dp3**2 - m2**2*p2dp3**2 
     .     - p1sq*p2dp3**2 + p1dp3**2*(-m2**2 + m3**2 + p2sq) 
     .     + 2._dp*p1dp2**2*p3sq + m2**2*p1sq*p3sq - m3**2*p1sq*p3sq 
     .     - p1sq*p2dp3*p3sq - m1**2*p2sq*p3sq + m2**2*p2sq*p3sq 
     .     + p1dp2*(-2._dp*p2dp3**2 + (-m1**2 + 2._dp*m2**2 - m3**2 + p1sq 
     .     + p2sq)*p3sq - p2dp3*(m3**2 - m4**2 + p3sq) + p1dp3*(m3**2 
     .     - m4**2 - 2._dp*p2dp3 + p3sq)) + p1dp3*((m1**2 - 2._dp*m2**2 
     .     + m3**2 - p1sq)*p2dp3 + p2sq*(m3**2 - m4**2 + p2dp3 + p3sq))
     .     )*xI3(s12,p3sq,p4sq,m1**2,m3**2,m4**2,musq,ep) 
     .     + (-2._dp*m2**2*m3**2*p1sq*p2dp3 + 2._dp*m3**4*p1sq*p2dp3 
     .     + 2._dp*m2**2*m4**2*p1sq*p2dp3 - 2._dp*m3**2*m4**2*p1sq*p2dp3 
     .     - m1**4*p2dp3**2 + 2._dp*m1**2*m2**2*p2dp3**2 - m2**4*p2dp3**2 
     .     + 2._dp*m1**2*p1sq*p2dp3**2 - 2._dp*m2**2*p1sq*p2dp3**2 
     .     + 4._dp*m3**2*p1sq*p2dp3**2 - p1sq**2*p2dp3**2 
     .     + m3**4*p1sq*p2sq - 2._dp*m3**2*m4**2*p1sq*p2sq 
     .     + m4**4*p1sq*p2sq + 2._dp*m3**2*p1sq*p2dp3*p2sq 
     .     - 2._dp*m4**2*p1sq*p2dp3*p2sq - p1dp3**2*((m2**2 - m3**2)**2 
     .     - 2._dp*(m2**2 + m3**2)*p2sq + p2sq**2) + m2**4*p1sq*p3sq 
     .     - 2._dp*m2**2*m3**2*p1sq*p3sq + m3**4*p1sq*p3sq 
     .     - 2._dp*m2**2*p1sq*p2dp3*p3sq + 2._dp*m3**2*p1sq*p2dp3*p3sq 
     .     + m1**4*p2sq*p3sq - 2._dp*m1**2*m2**2*p2sq*p3sq 
     .     + m2**4*p2sq*p3sq - 2._dp*m1**2*p1sq*p2sq*p3sq 
     .     - 2._dp*m4**2*p1sq*p2sq*p3sq + p1sq**2*p2sq*p3sq 
     .     + 2._dp*p1sq*p2dp3*p2sq*p3sq + p1sq*p2sq**2*p3sq 
     .     + p1sq*p2sq*p3sq**2 - p1dp2**2*((m3**2 - m4**2)**2 
     .     + 4._dp*p2dp3**2 - 2._dp*(2._dp*m2**2 - m3**2 + m4**2)*p3sq 
     .     + p3sq**2 + 4._dp*p2dp3*(m3**2 - m4**2 + p3sq)) 
     .     - 2._dp*p1dp3*(m1**2 - m2**2 - p1sq)*((-m2**2 + m3**2)*p2dp3 
     .     + p2sq*(m3**2 - m4**2 + p2dp3 + p3sq)) + 2._dp*p1dp2*((m1**2 
     .     - m2**2 - p1sq)*(2._dp*p2dp3**2 - (m2**2 - m3**2 + p2sq)*p3sq 
     .     + p2dp3*(m3**2 - m4**2 + p3sq)) + p1dp3*(-2._dp*(m2**2 
     .     + m3**2)*p2dp3 + (m2**2 - m3**2)*(m3**2 - m4**2 + p3sq) 
     .     + p2sq*(m3**2 - m4**2 + 2._dp*p2dp3 + p3sq))))*xI4(p1sq,p2sq,
     .     p3sq,p4sq,s12,s23,m1**2,m2**2,m3**2,m4**2,musq,ep)))/
     .     (-2._dp*p1dp2*p1dp3*p2dp3 + p1dp3**2*p2sq + p1dp2**2*p3sq 
     .     + p1sq*(p2dp3**2 - p2sq*p3sq))

      D6 = -2._dp*pi*D6

      end function D6


