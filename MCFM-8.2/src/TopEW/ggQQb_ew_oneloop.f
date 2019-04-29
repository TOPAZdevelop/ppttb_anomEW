      subroutine ggQQb_ew_oneloop(vew,s,beta,z)
C --- the analytical result is taken from [arXiv:hep-ph/0610335] 
C --- by J. H. Kuehn, A. Scharf and P. Uwer
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'scale.f'
      include 'anomcoup.f'
      real(dp):: vew,s,t,beta,z,vrts(5),slf(5),vrt(5),bx(5),alpha,
     .     genfacs,genfacself,genfacvert,genfacbox,genfacschi,sigma0,
     .     chifacs,chifac,gvt,gat,gvt_sq,gat_sq,xI1,xI2,xI3,xI4,sw2,cw2,
     .     mw,mz,mh,born,aemmz,trih,trizx,gw_sq,gab,g_rest
      real(dp) :: MARKUS_vrts3,MARKUS_slf3,MARKUS_vrt3,MARKUS_bx3
      real(dp) :: MARKUS_vrts4,MARKUS_slf4,MARKUS_vrt4,MARKUS_bx4
      real(dp) :: vol2,vol4,C33phiq3,C33phiu
      real(dp) :: db0
      integer ep,ierr
      complex(dp) :: cdb0,cdb0p
      common/em/aemmz

      ep = 0
c      musq = (2._dp*mt)**2
c      alpha = 1._dp/126.3_dp
c--- alpha -> standard MCFM value
      alpha=aemmz
      sw2 = xw

         
!          s=   173308.20441330347d0
!          beta=  0.55464632643236811d0
!          z=  0.61536288270299611d0

        mb=0.0001d0


      mw = wmass
      mz = zmass
      mh = hmass

c      sw2 = 1._dp - mw**2/mz**2
      cw2 = 1._dp-sw2

      gvt = (0.5_dp-4._dp*sw2/3._dp)/2._dp/sqrt(sw2*cw2)
      gat = 0.5_dp/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)

      gab = -gat 

c**********************************************************************************
c     MARKUS: add dim-6 operator contributions ( variables are in common block of anomcoup.f and set in mdata.f )       

        call ResetEWCouplings(sw2,gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest)

        vol2 = (vev/Lambda_BSM)**2
        vol4 = vol2**2        
        C33phiq3 = C_phiq_333 
        C33phiu = C_phiu_33        
        
c     END MARKUS      
c**********************************************************************************


      genfacs = Nc*mt**2/s**2*z**2/(1._dp-beta**2*z**2)
      genfacself = (2._dp - Nc**2*(1._dp-beta*z))/Nc/(1._dp-beta**2*z**2)
      genfacvert = (2._dp - Nc**2*(1._dp-beta*z))/Nc/(1._dp+beta*z)**2
      genfacbox = (Nc**2*(1._dp-beta*z)-2._dp)/Nc/(1._dp-beta**2*z**2)
      genfacschi = Nc*mt**4/s**2*z**2/(1._dp-beta**2*z**2)
      chifacs = 1d0/16._dp/sw2/cw2 ! gat_sq  MARKUS: replaced gat_sq coupl. by explicit expression
      chifac = (1d0/16._dp/sw2/cw2)/mz**2 !  MARKUS: replaced gat_sq coupl. by explicit expression
c --- sigma0 is jacobian additional to matrix element square
c --- sigma0 = gsq**2*beta/64._dp/pi/s/(Nc**2-1._dp)
      sigma0 = 1._dp
c      sigma0 = 1._dp/16._dp

      t = -s*(1._dp+beta*z)/2._dp+mt**2

C --- leading order differential cross section, use to normalize the result
C --- in the ref. (II.1)
      born = sigma0*(Nc**2*(1._dp+beta**2*z**2)-2._dp)/
     .     Nc/(1._dp-beta**2*z**2)**2*(1._dp-beta**4*z**4
     .     +2._dp*beta**2*(1._dp-beta**2)*(1._dp-z**2))
C --- the differential cross section is normalized when 
C --- sigma0 = (pi^2*alpha_s^2)/2 = gsq^2/32


C --- add triangle t,b contribution (II.16-17)
c MARKUS: replaced gw_sq coupl. by vev
      trih = 
     .     + 1d0/vev**2/pi**2*sigma0*mt**2*tevscale*
     .     beta**2/(1._dp - beta**2*z**2)/(s - mh**2)*
     .     ( tevscale*
     .     ( mt**2*(s - 4._dp*mt**2)*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep) - 2._dp*mt**2)
     .     + mb**2*(s - 4._dp*mb**2)*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep) - 2._dp*mb**2)

      
      trizx = 
     .     + 16._dp*alpha/pi*sigma0*mt**2/mz**2
     .     /(1._dp-beta**2*z**2)*(
     .     + gat_sq *mt**2*xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep)
     .     + gab*gat*mb**2*xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep)*0
     .     )

!        print *, "MCFM trizx",trizx
!        print *, "MAKU trizx",  -(alpha/Pi*
!      -      (-1 + beta**2)**2*s**2*
!      -      xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))
!      -      /(16.*CW2*MZ**2*SW2*(-1 + beta**2*z**2))
!        pause

C --- In the ref. (II.18-22), and Appendix B (B.1-15).
      vrts(1) = 
     .     (-2._dp*alpha*genfacs*sigma0*
     .     (beta**2*s*(2._dp*(gat_sq + gvt_sq)*mz**2 
     .     + (1._dp - beta**2)*(-3._dp*gat_sq + gvt_sq)*s)*
     .     db0(mt**2,mt**2,mz**2) + 2._dp*(2._dp*(gat_sq + gvt_sq)*mz**2 
     .     - beta**2*(-3._dp*gat_sq + gvt_sq)*s)*( - xI2(mt**2,
     .     mt**2,mz**2,musq,ep) + xI2(s,mt**2,mt**2,musq,ep)) + (4._dp*
     .     (gat_sq + gvt_sq)*mz**4 + 8._dp*beta**2*gat_sq*mz**2*s
     .     - beta**2*(gat_sq + gvt_sq + beta**2*(-3._dp*gat_sq + gvt_sq)
     .     )*s**2)*xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep)))/pi  



      vrts(2) = 
     .     (-8._dp*alpha*genfacs*gw_sq*sigma0*
     .     (-0.25_dp*beta**2*s*(4._dp*(mb**2 - mw**2) + (1._dp 
     .     - beta**2)*s)*db0(mt**2,mb**2,mw**2) + 0.5_dp*(4._dp*
     .     (-mb**2 + mw**2) + (1._dp + beta**2)*s)*
     .     (-xI2(mt**2,mb**2,mw**2,musq,ep) 
     .     + xI2(s,mb**2,mb**2,musq,ep)) 
     .     + 0.125_dp*(-4._dp*beta**2*s**2 + (4._dp*(-mb**2 + mw**2) 
     .     + (1._dp + beta**2)*s)**2)*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)))/pi


!     MARKUS: chi contribution, vertex, s-channel gluon
!       vrts(3) = 
!      .     (-8._dp*alpha*chifacs*genfacschi*sigma0*
!      .     (beta**2*s*db0(mt**2,mt**2,mz**2) + 2._dp*
!      .     (-xI2(mt**2,mt**2,mz**2,musq,ep) 
!      .     + xI2(s,mt**2,mt**2,musq,ep)) 
!      .     + (2._dp*mz**2 + beta**2*s)*
!      .     xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep)))/pi

           
          
        MARKUS_vrts3 = (alpha*MZ**2)/(16*MW**2*Pi*SW2)  /16d0  * ( 
     -     (-12*beta**2*(-1 + beta**2)*s*vol2*
     - (4*C33phiq3+4*C33phiq3**2*vol2+C33phiu*(2+C33phiu*vol2))*z**2)/MZ**2
     - +(24*beta**2*(-1+beta**2)**2*s*(1+2*C33phiq3*vol2 + C33phiu*vol2)**2*
     -     z**2*DB0(MT**2,MT**2,MZ**2))/(-1 + beta**2*z**2) + 
     -  (48*beta**2*(-1 + beta**2)*vol2*
     -     (4*C33phiq3**2*vol2 + C33phiu*(2 + C33phiu*vol2) + 
     -       C33phiq3*(4 + 8*C33phiu*vol2))*z**2*xI1(MT**2,musq,ep))/
     -   (MZ**2*(-1 + beta**2*z**2)) - 
     -  (48*beta**2*(-1 + beta**2)*vol2*
     -     (4*C33phiq3**2*vol2 + C33phiu*(2 + C33phiu*vol2) + 
     -       C33phiq3*(4 + 8*C33phiu*vol2))*z**2*xI1(MZ**2,musq,ep))/
     -   (MZ**2*(-1 + beta**2*z**2)) + 
     -  (24*(-1+beta**2)*(beta**2*(-1+beta**2)*(2*C33phiq3+C33phiu)*s*vol2*
     -        (2 + 2*C33phiq3*vol2 + C33phiu*vol2) + 
     -       2*MZ**2*((1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2 + 
     -          beta**2*(-1 + 4*C33phiq3*C33phiu*vol4)))*z**2*
     -     xI2(MT**2,MT**2,MZ**2,musq,ep))/(MZ**2*(-1 + beta**2*z**2)) + 
     -  (48*(-1 + beta**2)**2*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*z**2*
     -     xI2(s,MT**2,MT**2,musq,ep))/(-1 + beta**2*z**2) + 
     -  (24*(-1 + beta**2)**2*(2*MZ**2 + beta**2*s)*
     -     (1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*z**2*
     -     xI3(MT**2,MT**2,s,MT**2,MZ**2,MT**2,musq,ep))/(-1+beta**2*z**2))

!         print *, "MCFM vrts(3) ", vrts(3)
!         print *, "MAKU vrts(3) ", MARKUS_vrts3
!         pause             

              vrts(3) = MARKUS_vrts3
           
           

     
!     MARKUS: phi contribution vertex, s-channel gluon 
!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
!       vrts(4) = 
!      .     (-8._dp*alpha*genfacs*(0.125_dp/sw2)*sigma0*
!      .     (-0.03125_dp*(16._dp*beta**2*mb**4*s - 4._dp*beta**2*(1._dp 
!      .     - beta**2)*mw**2*s**2 + beta**2*(1._dp - beta**2 )**2*s**3
!      .     - 8._dp*mb**2*(2._dp*beta**2*mw**2*s + beta**2*(1._dp 
!      .     - beta**2)*s**2))*db0(mt**2,mb**2,mw**2) + 0.0625_dp*
!      .     (16._dp*mb**4 - 4._dp*(1._dp - beta**2)*mw**2*s 
!      .     - (1._dp - beta**4)*s**2 + mb**2*(-16._dp*mw**2 
!      .     + 8._dp*beta**2*s))*(xI2(mt**2,mb**2,mw**2,musq,ep) 
!      .     - xI2(s,mb**2,mb**2,musq,ep)) + 0.015625_dp*(64._dp*mb**6 
!      .     + 16._dp*(1._dp - beta**2)*mw**4*s + 8._dp*(1._dp 
!      .     - beta**4)*mw**2*s**2 + (1._dp - beta**2)**3*s**3 
!      .     - 16._dp*mb**4*(8._dp*mw**2 + (1._dp - beta**2)*s) 
!      .     + 4._dp*mb**2*(16._dp*mw**4 - (1._dp - beta**2)**2*
!      .     s**2))*xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)))
!      .     /(mw**2*pi)

    

        MARKUS_vrts4 = alpha/(64*MW**2*Pi*SW2)  /16d0  * (      
     -  -96*beta**2*(-1 + beta**2)*C33phiq3*s*vol2*(2 + C33phiq3*vol2)*z**2+
     - (24*beta**2*(-1 + beta**2)**2*s*(4*MW**2 + (-1 + beta**2)*s)*
     - (1 + C33phiq3*vol2)**2*z**2*DB0(MT**2,mb**2,MW**2))/(-1+beta**2*z**2)+ 
     - (192*beta**2*(-1+beta**2)*(1+C33phiq3*vol2*(4+C33phiq3*vol2))*z**2*
     -  0d0)/(-1 + beta**2*z**2) -
     - (384*beta**2*(-1 + beta**2)*C33phiq3*vol2*(2 + C33phiq3*vol2)*z**2*
     -  xI1(MW**2,musq,ep))/(-1 + beta**2*z**2) + 
     - (48*(-1 + beta**2)*(4*MW**2*
     -       ((1 + C33phiq3*vol2)**2 + 
     -  beta**2*(-1 + C33phiq3*vol2*(2 + C33phiq3*vol2))) + 
     - s*(1 -beta**4 + 2*(-1 + beta**2)**2*C33phiq3*vol2 + 
     - (-1+beta**2)**2*C33phiq3**2*vol4))*z**2*xI2(MT**2,0d0,MW**2,musq,ep))/
     - (-1+beta**2*z**2) + (48*(-1 + beta**2)**2*(4*MW**2 + s + beta**2*s)*
     - (1 + C33phiq3*vol2)**2*z**2*xI2(s,0d0,0d0,musq,ep))/(-1+beta**2*z**2)+
     - (12*(-1 + beta**2)**2*(4*MW**2 + (-1 + beta)**2*s)*
     - (4*MW**2 + (1 + beta)**2*s)*(1 + C33phiq3*vol2)**2*z**2*
     - xI3(MT**2,MT**2,s,0d0,MW**2,0d0,musq,ep))/(-1 + beta**2*z**2)       )

     
!         print *, "MCFM vrts(4) ", vrts(4)
!         print *, "MAKU vrts(4) ", MARKUS_vrts4
!         pause             

              vrts(4) = MARKUS_vrts4
     
     
     
!       vrts(5) = 
!      .     (4._dp*alpha*genfacs*gw_sq*mt**2*sigma0*
!      .     (beta**2*s*(-mh**2 + (1._dp - beta**2)*s)*db0(mt**2
!      .     ,mt**2,mh**2) + 2._dp*(mh**2 + beta**2*s)*(xI2(mt**2,mt**2,
!      .     mh**2,musq,ep) - xI2(s,mt**2,mt**2,musq,ep)) + 
!      .     (-2._dp*mh**4 - 3._dp*beta**2*mh**2*s + beta**2*(1._dp 
!      .     - beta**2)*s**2)*xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2
!      .     ,musq,ep)))/(mw**2*pi)
! 
!            print *, "vev1",1d0/dsqrt(8*3.14159265d0*alpha*gw_sq/mw**2)
!            print *, "vev2",vev

      vrts(5) =   !MARKUS: replaced gw_sq coupl. by vev
     .     (genfacs*mt**2/vev**2*sigma0*
     .     (beta**2*s*(-mh**2 + (1._dp - beta**2)*s)*db0(mt**2
     .     ,mt**2,mh**2) + 2._dp*(mh**2 + beta**2*s)*(xI2(mt**2,mt**2,
     .     mh**2,musq,ep) - xI2(s,mt**2,mt**2,musq,ep)) + 
     .     (-2._dp*mh**4 - 3._dp*beta**2*mh**2*s + beta**2*(1._dp 
     .     - beta**2)*s**2)*xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2
     .     ,musq,ep)))/(2*pi**2)
           
           
           
           

      slf(1) = 
     .     (0.125_dp*alpha*genfacself*sigma0*
     .     ((2._dp*(2._dp*(gat_sq + gvt_sq)*mz**2 + 
     .     (1._dp - beta**2)*(-3._dp*gat_sq + gvt_sq)*s)*
     .     (1._dp - beta**4*z**4 + (1._dp - beta**2)*
     .     (2._dp*beta**2 - beta*z - 3._dp*beta**2*z**2))*
     .     db0(mt**2,mt**2,mz**2))/(1._dp + beta*z) + 
     .     (16._dp*(gat_sq + gvt_sq)*(1._dp - beta**4*z**4 + 
     .     beta**2*(1._dp - beta**2)*(1._dp - 3._dp*z**2))*
     .     (-xI1(mt**2,musq,ep) + xI1(mz**2,musq,ep)))/ 
     .     ((1._dp - beta**2)*s*(1._dp + beta**2 + 2._dp*beta*z))  
     .     - (4._dp*((1._dp - beta**2)*(-3._dp*gat_sq + gvt_sq)*
     .     (1._dp + beta**2 - 2._dp*beta**4 -
     .     beta**2*(2._dp - 3._dp*beta**2)*z**2 - beta**4*z**4) 
     .     + (2._dp*(gat_sq + gvt_sq)*mz**2*(2._dp + 2._dp*beta**2 - 
     .     5._dp*beta**4 + 2._dp*beta**6 + beta**3*(3._dp - 2._dp*beta**2)*z 
     .     - 3._dp*beta**2*(2._dp - 3._dp*beta**2 + beta**4)*z**2 - 
     .     3._dp*beta**3*(1._dp - beta**2)*z**3 - beta**4*(2._dp 
     .     - beta**2)*z**4 - beta**5*z**5))/((1._dp 
     .     - beta**2)*s))*xI2(mt**2,mt**2,mz**2,musq,ep))/(1._dp 
     .     + beta*z)**2 + (4._dp*((2._dp*beta**2*(gat_sq + gvt_sq)*mz**2*
     .     (1._dp - z**2)*(2._dp + beta**2 - 2._dp*beta**4 + 
     .     (3._dp*beta - 2._dp*beta**3)*z + beta**4*z**2 + beta**3*z**3))
     .     /s + gvt_sq*(2._dp + 2._dp*beta**2 - 4._dp*beta**4 - 
     .     beta**6 + 2._dp*beta**8 + beta*(4._dp + 2._dp*beta**2 - 
     .     8._dp*beta**4 + 4._dp*beta**6)*z + beta**2*(-4._dp + 
     .     7._dp*beta**2 + beta**4 - 3._dp*beta**6)*z**2 - beta**3*
     .     (10._dp - 16._dp*beta**2 + 6._dp*beta**4)*z**3 + beta**4*(-5._dp 
     .     + 3._dp*beta**2 + beta**4)*z**4 - beta**5*(4._dp 
     .     - 2._dp*beta**2)*z**5 - beta**6*z**6) - gat_sq*
     .     (2._dp + 2._dp*beta**2 - 8._dp*beta**4 - 3._dp*beta**6 
     .     + 6._dp*beta**8 + 2._dp*beta*(2._dp - beta**2 
     .     - 8._dp*beta**4 + 6._dp*beta**6)*z - beta**2*(4._dp 
     .     - 5._dp*beta**2 - 7._dp*beta**4 + 9._dp*beta**6)*z**2 
     .     - 6._dp*beta**3*(1._dp - 4._dp*beta**2 + 3._dp*beta**4)*z**3 
     .     + beta**4*(1._dp - 3._dp*beta**2 + 3._dp*beta**4)*z**4 
     .     - beta**5*(4._dp - 6._dp*beta**2)*z**5 + beta**6*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta*z)**2*(1._dp 
     .     + beta**2 + 2._dp*beta*z))))/pi



      slf(2) = 
     .     (alpha*genfacself*gw_sq*sigma0*
     .     ((0.25_dp*(-4._dp*mb**2 + 4._dp*mw**2 - s + beta**2*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*
     .     db0(mt**2,mb**2,mw**2))/(1._dp + beta*z) + (4._dp*(1._dp - 
     .     beta**4*z**4 + beta**2*(1._dp - beta**2)*(1._dp - 
     .     3._dp*z**2))*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))
     .     /((1._dp - beta**2)*s*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.5_dp*((-4._dp*(mb**2 - mw**2)*(-2._dp - 2._dp*beta**2 + 
     .     5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp - 
     .     2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 + 
     .     beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 + 
     .     beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5))/(1._dp - 
     .     beta**2) - beta**2*s*(1._dp - z**2)*(2._dp + 
     .     beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z + 
     .     beta**3*z**2*(beta + z)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     (s*(1._dp + beta*z)**2) + (0.5_dp*beta**2*(1._dp - z**2)*
     .     (2._dp + beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z))*(4._dp*(-mb**2 + mw**2) + 
     .     s*(1._dp + beta**2 + 2._dp*beta*z))*xI2(t,mb**2,mw**2,musq,ep))
     .     /(s*(1._dp + beta*z)**2*(1._dp + beta**2 + 2._dp*beta*z))))/pi



!       slf(3) = 
!      .     (2._dp*alpha*chifac*genfacself*sigma0*
!      .     ((0.125_dp*(1._dp - beta**2)*mz**2*s*(1._dp 
!      .     - beta**4*z**4 + beta*(1._dp - beta**2)*
!      .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mt**2,
!      .     mz**2))/(1._dp + beta*z) + (0.5_dp*(1._dp - beta**4*z**4 
!      .     + (1._dp - beta**2)*(beta**2 - 3._dp*beta**2*z**2))*
!      .     (-xI1(mt**2,musq,ep) + xI1(mz**2,musq,ep)))/(1._dp 
!      .     + beta**2 + 2._dp*beta*z) - (0.25_dp*mz**2*(2._dp + 2._dp*beta**2 
!      .     - 5._dp*beta**4 + 2._dp*beta**6 + beta**3*(3._dp 
!      .     - 2._dp*beta**2)*z - 3._dp*beta**2*(2._dp - 3._dp*beta**2 
!      .     + beta**4)*z**2 - 3._dp*beta**3*(1._dp - beta**2)*z**3 
!      .     - beta**4*(2._dp - beta**2)*z**4 
!      .     - beta**5*z**5)*xI2(mt**2,mt**2,mz**2,musq,ep))/
!      .     (1._dp + beta*z)**2 + (0.125_dp*(1._dp - beta**2)*
!      .     (s*(1._dp + beta**2 - beta**4 - 3._dp*beta**2*(1._dp 
!      .     - beta**2)*z**2 - beta**4*z**4) + 
!      .     (2._dp*beta**2*mz**2*(1._dp - z**2)*(2._dp + beta**2 
!      .     - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
!      .     + beta**3*z**2*(beta + z)))/(1._dp + beta*z)**2)*
!      .     xI2(t,mt**2,mz**2,musq,ep))/
!      .     (1._dp + beta**2 + 2._dp*beta*z)))/pi


          MARKUS_slf3 = (alpha)/(64*MW**2*Pi*SW2)/16d0 *(
     -     (8*s*(-3 + beta*z)*(-7 + 9*beta*z)*
     -     (4*C33phiq3*(beta**2*C33phiu*vol2**2 - 
     -          beta**2*C33phiu*vol4 + 
     -          (-1 + beta**2)*vol2*(1 + beta*z)) + 
     -       4*beta**2*C33phiq3**2*
     -        (-(vol4*(-2 + beta**2 + z**2)) + 
     -          vol2**2*(-1 + beta*(beta + z - z**3))) + 
     -       C33phiu*(2*(-1 + beta**2)*vol2*(1 + beta*z) - 
     -          beta**2*C33phiu*vol4*(-2 + beta**2 + z**2) + 
     -          beta**2*C33phiu*vol2**2*(-1 + beta*(beta + z - z**3))
     -          )))/(3.*(-1 + beta*z)*(1 + beta*z)**2) - 
     -  (16*(-1 + beta)*(1 + beta)*s*
     -     (MZ + (2*C33phiq3 + C33phiu)*MZ*vol2)**2*(-7 + 9*beta*z)*
     -     (-1 + beta*(z + beta*
     -           (-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4))
     -          ))*DB0(MT**2,MT**2,MZ**2))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**2) + 
     -  (32*(-7 + 9*beta*z)*
     -     (2 + 8*C33phiq3*vol2 - 
     -       4*C33phiu*vol2*(-1 + C33phiq3*vol2) + 
     -       2*(4*C33phiq3**2 + 6*C33phiq3*C33phiu + C33phiu**2)*
     -        vol4 + 4*beta*
     -        (1 + 4*C33phiq3*vol2 + 2*C33phiu*vol2 - 
     -          2*C33phiq3*C33phiu*vol2**2 + 
     -          (4*C33phiq3**2 + 6*C33phiq3*C33phiu + C33phiu**2)*
     -           vol4)*z + beta**7*
     -        (2*(2*C33phiq3 + C33phiu)*vol2*
     -           (-1 + 2*C33phiq3*vol2 + C33phiu*vol2) - 
     -          (12*C33phiq3**2 + 16*C33phiq3*C33phiu + 
     -             3*C33phiu**2)*vol4)*z*(2 - 3*z**2 + z**4) + 
     -       beta**5*z*(-4 - 20*C33phiq3*vol2 - 10*C33phiu*vol2 + 
     -          24*C33phiq3*C33phiu*vol2**2 - 20*C33phiq3**2*vol4 - 
     -          48*C33phiq3*C33phiu*vol4 - 5*C33phiu**2*vol4 + 
     -          2*(6 + vol2*
     -              (36*C33phiq3**2*vol2 + 
     -                4*C33phiq3*(7 + 8*C33phiu*vol2) + 
     -                C33phiu*(14 + 9*C33phiu*vol2)) - 
     -             2*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 - 
     -          (4 + 20*C33phiq3**2*vol4 + 
     -             5*C33phiu*(2*vol2 + C33phiu*vol4) + 
     -             4*C33phiq3*
     -              (5*vol2 - 2*C33phiu*vol2**2 + 8*C33phiu*vol4))*
     -           z**4) + beta**4*
     -        (-2 - 4*C33phiq3*vol2 - 2*C33phiu*vol2 + 
     -          12*C33phiq3*C33phiu*vol2**2 - 4*C33phiq3**2*vol4 - 
     -          12*C33phiq3*C33phiu*vol4 - C33phiu**2*vol4 + 
     -          (8 + 4*C33phiq3**2*(9*vol2**2 + vol4) + 
     -             8*C33phiq3*
     -              (5*vol2*(1 + C33phiu*vol2) + C33phiu*vol4) + 
     -             C33phiu*(vol2*(20 + 9*C33phiu*vol2) + 
     -                C33phiu*vol4))*z**2 - 
     -          (8 + 44*C33phiq3**2*vol4 + 
     -             11*C33phiu*(2*vol2 + C33phiu*vol4) + 
     -             C33phiq3*
     -              (44*vol2 - 4*C33phiu*vol2**2 + 60*C33phiu*vol4))*
     -           z**4) + beta**6*
     -        (2*vol2*(4*C33phiq3**2*vol2 + 
     -             C33phiu*(-2 + C33phiu*vol2) + 
     -             C33phiq3*(-4 + 6*C33phiu*vol2)) - 
     -          4*(4*C33phiq3**2 + 7*C33phiq3*C33phiu + C33phiu**2)*
     -           vol4 + (-2 + 4*C33phiq3**2*(6*vol2**2 - 7*vol4) + 
     -             4*C33phiq3*
     -              (vol2*(-1 + 8*C33phiu*vol2) - 8*C33phiu*vol4) + 
     -             C33phiu*(-2*vol2 + 6*C33phiu*vol2**2 - 
     -                7*C33phiu*vol4))*z**2 + 
     -          (6 + (2*C33phiq3 + C33phiu)*vol2*
     -              (16 + 2*C33phiq3*vol2 + C33phiu*vol2) + 
     -             (28*C33phiq3**2 + 36*C33phiq3*C33phiu + 
     -                7*C33phiu**2)*vol4)*z**4 - 
     -          (2 + 6*(2*C33phiq3 + C33phiu)*vol2 + 
     -             (12*C33phiq3**2 + 16*C33phiq3*C33phiu + 
     -                3*C33phiu**2)*vol4)*z**6) + 
     -       beta**3*z*(4 - 12*z**2 + 
     -          4*C33phiq3**2*vol4*(9 - 17*z**2) + 
     -          C33phiu*(2*vol2 + C33phiu*vol4)*(9 - 17*z**2) + 
     -          4*C33phiq3*(9*vol2 + 14*C33phiu*vol4 + 
     -             (vol2*(-17 + 4*C33phiu*vol2) - 26*C33phiu*vol4)*
     -              z**2)) + 
     -       2*beta**2*(1 - 2*z**2 + 
     -          4*C33phiq3**2*vol4*(2 - 3*z**2) - 
     -          C33phiu*(2*vol2 + C33phiu*vol4)*(-2 + 3*z**2) + 
     -          2*C33phiq3*(4*vol2 - C33phiu*vol2**2 + 
     -             7*C33phiu*vol4 + 
     -             2*(vol2*(-3 + C33phiu*vol2) - 5*C33phiu*vol4)*z**2
     -             )))*xI1(MT**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z))
     -   + (32*(-7 + 9*beta*z)*
     -     (-2*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2 - 
     -       2*beta*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*z + 
     -       beta**3*z*(-2 - 20*C33phiq3*vol2 - 10*C33phiu*vol2 - 
     -          32*C33phiq3**2*vol2**2 - 
     -          32*C33phiq3*C33phiu*vol2**2 - 8*C33phiu**2*vol2**2 + 
     -          12*C33phiq3**2*vol4 + 3*C33phiu**2*vol4 + 
     -          3*(2 + 2*(2*C33phiq3 + C33phiu)*vol2*
     -              (3 + 4*C33phiq3*vol2 + 2*C33phiu*vol2) - 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4)*z**2) + 
     -       2*beta**2*(-1 - 8*C33phiq3*vol2 - 4*C33phiu*vol2 - 
     -          12*C33phiq3**2*vol2**2 - 
     -          12*C33phiq3*C33phiu*vol2**2 - 3*C33phiu**2*vol2**2 + 
     -          4*C33phiq3**2*vol4 + C33phiu**2*vol4 + 
     -          (3 + (2*C33phiq3 + C33phiu)*vol2*
     -              (8 + 5*(2*C33phiq3 + C33phiu)*vol2) - 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4)*z**2) + 
     -       beta**6*(2*(2*C33phiq3 + C33phiu)*vol2*
     -           (1 + 2*C33phiq3*vol2 + C33phiu*vol2) - 
     -          (4*C33phiq3**2 + C33phiu**2)*vol4)*
     -        (2 - 3*z**2 + z**4) + 
     -       beta**4*(2 + 4*C33phiq3*vol2 + 2*C33phiu*vol2 + 
     -          4*C33phiq3**2*vol4 + C33phiu**2*vol4 - 
     -          (6 + 2*(2*C33phiq3 + C33phiu)*vol2*
     -              (5 + 4*C33phiq3*vol2 + 2*C33phiu*vol2) + 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -          2*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*z**4) + 
     -       beta**5*z*(2 + 
     -          2*(2*C33phiq3 + C33phiu)*vol2*
     -           (4 + 6*C33phiq3*vol2 + 3*C33phiu*vol2) - 
     -          2*(4*C33phiq3**2 + C33phiu**2)*vol4 - 
     -          3*(2 + 2*(2*C33phiq3 + C33phiu)*vol2*
     -              (3 + 4*C33phiq3*vol2 + 2*C33phiu*vol2) - 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -          (2 + 2*(2*C33phiq3 + C33phiu)*vol2*
     -              (3 + 4*C33phiq3*vol2 + 2*C33phiu*vol2) - 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4)*z**4))*
     -     xI1(MZ**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z))
     -   + (8*(-7 + 9*beta*z)*
     -     (4*MZ**2*(2 + 5*(2*C33phiq3 + C33phiu)*vol2 + 
     -          2*(4*C33phiq3**2 + 6*C33phiq3*C33phiu + C33phiu**2)*
     -           vol4 - 3*beta**3*
     -           (1 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -             (4*C33phiq3**2 + 8*C33phiq3*C33phiu + C33phiu**2)*
     -              vol4)*z*(-1 + z**2) + 
     -          2*beta**2*(1 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -             4*C33phiq3**2*vol4 + 8*C33phiq3*C33phiu*vol4 + 
     -             C33phiu**2*vol4 - 
     -             (3 + 8*(2*C33phiq3 + C33phiu)*vol2 + 
     -                (6*C33phiq3 + C33phiu)*
     -                 (2*C33phiq3 + 3*C33phiu)*vol4)*z**2) + 
     -          beta**6*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*
     -           (2 - 3*z**2 + z**4) - 
     -          beta**5*(1 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -             (4*C33phiq3**2 + 8*C33phiq3*C33phiu + C33phiu**2)*
     -              vol4)*z*(2 - 3*z**2 + z**4) + 
     -          beta**4*(-5 - 24*C33phiq3*vol2 - 12*C33phiu*vol2 - 
     -             20*C33phiq3**2*vol4 - 28*C33phiq3*C33phiu*vol4 - 
     -             5*C33phiu**2*vol4 + 
     -             (9 + vol2*
     -                 (36*C33phiq3**2*vol2 + 
     -                   C33phiu*(22 + 9*C33phiu*vol2) + 
     -                   C33phiq3*(44 + 52*C33phiu*vol2)))*z**2 - 
     -             (2 + 5*(2*C33phiq3 + C33phiu)*vol2 + 
     -                2*(4*C33phiq3**2 + 6*C33phiq3*C33phiu + 
     -                   C33phiu**2)*vol4)*z**4)) + 
     -       2*(2*C33phiq3 + C33phiu)*s*
     -        (-((2*C33phiq3 + C33phiu)*vol4) - 
     -          beta**2*(2*C33phiq3 + C33phiu)*vol2**2*
     -           (1 - 4*beta**2 + 2*beta**4 + 
     -             beta*(3 - 5*beta**2 + 2*beta**4)*z - 
     -             4*(-1 + beta**2)**2*z**2 - 
     -             3*beta*(-1 + beta**2)**2*z**3 + 
     -             beta**2*(-1 + beta**2)*z**4 + 
     -             beta**3*(-1 + beta**2)*z**5) - 
     -          2*(-1 + beta**2)*vol2*(1 + beta*z)*
     -           (-1 + beta*
     -              (z + beta*
     -                 (-2 - beta*z + 3*z**2 + 
     -                   beta**2*(2 - 3*z**2 + z**4))))))*
     -     xI2(MT**2,MT**2,MZ**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**3) - 
     -  (16*(-7 + 9*beta*z)*
     -     (s + (2*C33phiq3 + C33phiu)*s*vol2 + 
     -       2*MZ**2*vol2*(C33phiu + 
     -          C33phiq3*(2 + 4*C33phiu*vol2)) + 
     -       beta*(s*(2 + 6*C33phiq3*vol2 + 3*C33phiu*vol2) + 
     -          4*MZ**2*vol2*
     -           (C33phiu + C33phiq3*(2 + 4*C33phiu*vol2)))*z + 
     -       beta**3*z*(2*MZ**2*
     -           (3 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 - 
     -             8*C33phiq3**2*vol4 - 2*C33phiu**2*vol4 + 
     -             (-3 - 7*C33phiu*vol2 - 
     -                2*C33phiq3*vol2*(7 + 8*C33phiu*vol2) + 
     -                8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4)*z**2)
     -           + s*(-4*(2*C33phiq3 + C33phiu)*vol2 - 
     -             7*(4*C33phiq3**2 + C33phiu**2)*vol4 + 
     -             (-6 - (2*C33phiq3 + C33phiu)*vol2 + 
     -                7*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2)) - 
     -       beta**5*z*(2*MZ**2*
     -           (5 + 22*C33phiq3*vol2 + 11*C33phiu*vol2 + 
     -             24*C33phiq3*C33phiu*vol2**2 + 
     -             16*C33phiq3**2*vol4 + 4*C33phiu**2*vol4 - 
     -             2*(3 + 8*C33phiu*vol2 + 
     -                4*C33phiq3*vol2*(4 + 5*C33phiu*vol2) + 
     -                8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4)*z**2 + 
     -             (1 + 3*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(3 + 4*C33phiu*vol2))*z**4) + 
     -          s*(4 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -             2*(-6 - 3*(2*C33phiq3 + C33phiu)*vol2 + 
     -                2*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -             2*(1 + 2*C33phiq3*vol2 - 8*C33phiq3**2*vol4 + 
     -                C33phiu*(vol2 - 2*C33phiu*vol4))*z**4)) + 
     -       beta**8*(2*(MZ + (2*C33phiq3 + C33phiu)*MZ*vol2)**2*
     -           (2 - 3*z**2 + z**4) + 
     -          s*z**2*(1 + 4*C33phiq3*vol2 + 2*C33phiu*vol2 + 
     -             8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4 - 
     -             3*(1 + 2*C33phiq3*vol2 + 4*C33phiq3**2*vol4 + 
     -                C33phiu*(vol2 + C33phiu*vol4))*z**2 + 
     -             (1 + 2*C33phiq3*vol2 + 4*C33phiq3**2*vol4 + 
     -                C33phiu*(vol2 + C33phiu*vol4))*z**4)) + 
     -       beta**7*z*(2*MZ**2*
     -           (1 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -             2*(2*C33phiq3 + C33phiu)**2*vol4)*
     -           (2 - 3*z**2 + z**4) + 
     -          s*(2 + 8*C33phiq3*vol2 + 16*C33phiq3**2*vol4 + 
     -             4*C33phiu*(vol2 + C33phiu*vol4) - 
     -             (6 + 5*(2*C33phiq3 + C33phiu)*vol2 + 
     -                4*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -             (2 + 2*(2*C33phiq3 + C33phiu)*vol2 - 
     -                (4*C33phiq3**2 + C33phiu**2)*vol4)*z**4 + 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4*z**6)) - 
     -       beta**2*(s*(2*C33phiq3*vol2 + C33phiu*vol2 + 
     -             8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4 - 
     -             2*(-1 + 2*C33phiq3*vol2 + 4*C33phiq3**2*vol4 + 
     -                C33phiu*(vol2 + C33phiu*vol4))*z**2) + 
     -          2*MZ**2*(-2 - 3*C33phiu*vol2 + 
     -             2*(1 + 2*C33phiu*vol2)*z**2 + 
     -             2*C33phiq3*vol2*
     -              (-3 - 2*C33phiu*vol2 + 4*(1 + C33phiu*vol2)*z**2)
     -             )) + beta**6*
     -        (s + s*(2*(2*C33phiq3*vol2 + 4*C33phiq3**2*vol4 + 
     -                C33phiu*(vol2 + C33phiu*vol4)) - 
     -             (5 + 6*C33phiq3*vol2 + 3*C33phiu*vol2)*z**2 + 
     -             (7 + 5*(2*C33phiq3 + C33phiu)*vol2 - 
     -                4*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**4 + 
     -             (-1 - (2*C33phiq3 + C33phiu)*vol2 + 
     -                2*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**6) + 
     -          2*MZ**2*(-3 + 4*z**2 - z**4 + 
     -             4*C33phiq3**2*vol4*(-1 + z**2)**3 + 
     -             C33phiu**2*vol4*(-1 + z**2)**3 - 
     -             C33phiu*vol2*(4 - 6*z**2 + z**4) + 
     -             2*C33phiq3*vol2*
     -              (-4 + 6*z**2 - z**4 + 
     -                2*C33phiu*vol2*(-1 + 2*z**2)))) + 
     -       beta**4*(s*(-2 + 6*z**2 - 4*z**4 - 
     -             2*C33phiq3*vol2*(2 + z**2 + 2*z**4) - 
     -             C33phiu*vol2*(2 + z**2 + 2*z**4) + 
     -             4*C33phiq3**2*vol4*(-1 - 7*z**2 + 8*z**4) + 
     -             C33phiu**2*vol4*(-1 - 7*z**2 + 8*z**4)) + 
     -          2*MZ**2*(-1 + z**2 - C33phiu*vol2*(-2 + z**2)**2 + 
     -             4*C33phiq3**2*vol4*(-2 - z**2 + 3*z**4) + 
     -             C33phiu**2*vol4*(-2 - z**2 + 3*z**4) - 
     -             2*C33phiq3*vol2*
     -              ((-2 + z**2)**2 + 
     -                2*C33phiu*vol2*(3 - 3*z**2 + z**4)))))*
     -     xI2(MT**2 - (s*(1d0+beta*z))/2d0,MT**2,MZ**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z))          
     -   )

!          print *, "MCFM slf(3)",slf(3)
!          print *, "MAKU slf(3)",MARKUS_slf3   
!          pause 
         
         slf(3) = MARKUS_slf3



!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
!       slf(4) = 
!      .     (alpha*genfacself*(0.125_dp/sw2)*sigma0*
!      .     ((-0.03125_dp*(16._dp*mb**2*(mb**2 - mw**2) 
!      .     - 4._dp*(1._dp - beta**2)*(2._dp*mb**2 + mw**2)*s 
!      .     + (1._dp - beta**2)**2*s**2)*(1._dp - beta**4*z**4 
!      .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
!      .     - 3._dp*beta*z**2))*db0(mt**2,mb**2,mw**2))/(1._dp + beta*z) 
!      .     + (0.5_dp*(4._dp*mb**2 + (1._dp - beta**2)*s)*(1._dp 
!      .     - beta**4*z**4 + beta**2*(1._dp - beta**2)*(1._dp 
!      .     - 3._dp*z**2))*(-xI1(mb**2,musq,ep) 
!      .     + xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
!      .     (1._dp + beta**2 + 2._dp*beta*z)) - (0.0625_dp*(8._dp*(1._dp 
!      .     - beta**2)**2*mb**2*s*(1._dp + beta**2 - 2._dp*beta**4 
!      .     - beta**2*(2._dp - 3._dp*beta**2)*z**2 
!      .     - beta**4*z**4) - 4._dp*(4._dp*mb**2*(-mb**2 + mw**2) 
!      .     + (1._dp - beta**2)*mw**2*s)*(-2._dp - 2._dp*beta**2 
!      .     + 5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp 
!      .     - 2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 
!      .     + beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 
!      .     + beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5) 
!      .     + beta**2*(1._dp - beta**2)**2*s**2*(1._dp - z**2)*
!      .     (2._dp + beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
!      .     + beta**3*z**2*(beta + z)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
!      .     ((1._dp - beta**2)*s*(1._dp + beta*z)**2) + 
!      .     (0.0625_dp*(-8._dp*mb**2*s*(-2._dp - 2._dp*beta**2 + 4._dp*beta**4 
!      .     + beta**6 - 2._dp*beta**8 - 2._dp*beta*(2._dp + beta**2 
!      .     - 4._dp*beta**4 + 2._dp*beta**6)*z + beta**2*(4._dp 
!      .     - 7._dp*beta**2 - beta**4 + 3._dp*beta**6)*z**2 
!      .     + 2._dp*beta**3*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**3 
!      .     + beta**4*(5._dp - 3._dp*beta**2 - beta**4)*z**4 
!      .     + 2._dp*beta**5*(2._dp - beta**2)*z**5 + beta**6*z**6) 
!      .     - beta**2*(1._dp - z**2)*(2._dp + beta**2 
!      .     - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
!      .     + beta**3*z**2*(beta + z))*(16._dp*mb**2*(mb**2 - mw**2) 
!      .     - 4._dp*(1._dp - beta**2)*mw**2*s - s**2*(1._dp 
!      .     - beta**4 + 2._dp*beta*(1._dp - beta**2)*z)))*
!      .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta*z)**2*
!      .     (1._dp + beta**2 + 2._dp*beta*z))))/(mw**2*pi)


         

          MARKUS_slf4=alpha/(64*MW**2*Pi*SW2)/16d0*(
     -   (-16*C33phiq3*s*(-3 + beta*z)*(-7 + 9*beta*z)*
     -     (-2*(-1 + beta**2)*vol2*(1 + beta*z) + 
     -       beta**2*C33phiq3*vol4*(-1 + z**2) + 
     -       beta**3*C33phiq3*vol2**2*z*(-1 + z**2)))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**2) - 
     -  (4*(-1 + beta**2)*s*(4*MW**2 + (-1 + beta**2)*s)*
     -     (1 + C33phiq3*vol2)**2*(-7 + 9*beta*z)*
     -     (-1 + beta*(z + beta*
     -           (-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4)))
     -       )*DB0(MT**2,mb**2,MW**2))/(3.*(-1 + beta*z)*(1 + beta*z)**2)
     -   + (64*(-7 + 9*beta*z)*
     -     (-1 - C33phiq3*(2*vol2 + C33phiq3*vol4) - 
     -       beta*(1 + 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z + 
     -       beta**6*C33phiq3*(2*vol2 + C33phiq3*vol4)*
     -        (2 - 3*z**2 + z**4) + 
     -       beta**4*(1 - (3 + 4*C33phiq3*vol2 + 2*C33phiq3**2*vol4)*
     -           z**2 + (1 + 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z**4)
     -         + beta**3*z*(-1 + 3*z**2 + 
     -          2*C33phiq3*(2*vol2 + C33phiq3*vol4)*(-2 + 3*z**2)) + 
     -       beta**2*(-1 + 3*z**2 + 
     -          C33phiq3*(2*vol2 + C33phiq3*vol4)*(-3 + 5*z**2)) + 
     -       beta**5*z*(1 - 3*z**2 + z**4 + 
     -          2*C33phiq3*vol2*(3 - 6*z**2 + 2*z**4) + 
     -          C33phiq3**2*vol4*(3 - 6*z**2 + 2*z**4)))*
     -     xI1(MW**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z))
     -   + (8*(-7 + 9*beta*z)*
     -     (4*MW**2*(2 + 6*C33phiq3*vol2 + 2*C33phiq3**2*vol4 - 
     -          3*beta**3*(1 + 4*C33phiq3*vol2 + C33phiq3**2*vol4)*z*
     -           (-1 + z**2) + 
     -          2*beta**2*(1 + 4*C33phiq3*vol2 + C33phiq3**2*vol4 - 
     -             (3 + 10*C33phiq3*vol2 + 3*C33phiq3**2*vol4)*z**2)
     -           + beta**6*(1 + C33phiq3*vol2)**2*
     -           (2 - 3*z**2 + z**4) - 
     -          beta**5*(1 + 4*C33phiq3*vol2 + C33phiq3**2*vol4)*z*
     -           (2 - 3*z**2 + z**4) + 
     -          beta**4*(-5 - 14*C33phiq3*vol2 - 5*C33phiq3**2*vol4 + 
     -             (9 + 26*C33phiq3*vol2 + 9*C33phiq3**2*vol4)*z**2 - 
     -             2*(1 + 3*C33phiq3*vol2 + C33phiq3**2*vol4)*z**4))
     -        + (-1 + beta**2)*s*
     -        (2*C33phiq3*vol2 + 
     -          beta**2*(-2 - 2*C33phiq3**2*vol4 + 
     -             2*(1 - 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z**2 + 
     -             3*beta*(1 + C33phiq3**2*vol4)*z*(-1 + z**2) + 
     -             beta**4*(1 + C33phiq3*vol2)**2*
     -              (2 - 3*z**2 + z**4) + 
     -             beta**3*(1 + C33phiq3**2*vol4)*z*
     -              (2 - 3*z**2 + z**4) + 
     -             beta**2*(-1 + z**2 + 
     -                C33phiq3**2*vol4*(-1 + z**2) - 
     -                2*C33phiq3*vol2*(3 - 5*z**2 + z**4)))))*
     -     xI2(MT**2,0d0,MW**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**3) - 
     -  (8*(-7 + 9*beta*z)*(4*MW**2 + s*(1 + beta**2 + 2*beta*z))*
     -     (2*C33phiq3*vol2 + 4*beta*C33phiq3*vol2*z + 
     -       beta**3*z*(3 - 7*C33phiq3**2*vol4 + 
     -          (-3 - 8*C33phiq3*vol2 + 7*C33phiq3**2*vol4)*z**2) + 
     -       beta**8*(1 + C33phiq3*vol2)**2*(2 - 3*z**2 + z**4) + 
     -       beta**7*(1 + 4*C33phiq3*vol2 + 3*C33phiq3**2*vol4)*z*
     -        (2 - 3*z**2 + z**4) + 
     -       beta**5*z*(-5 + 6*z**2 - z**4 - 
     -          4*C33phiq3*vol2*(3 - 5*z**2 + z**4) + 
     -          C33phiq3**2*vol4*(-3 + 2*z**2 + z**4)) + 
     -       beta**4*(-1 + z**2 - 
     -          2*C33phiq3*vol2*(3 - 3*z**2 + z**4) + 
     -          3*C33phiq3**2*vol4*(-1 - z**2 + 2*z**4)) + 
     -       beta**6*(-3 + 4*z**2 - z**4 + 
     -          2*C33phiq3*vol2*(-1 + 2*z**2) + 
     -          C33phiq3**2*vol4*(1 + 2*z**2 - 5*z**4 + 2*z**6)) + 
     -       2*beta**2*(1 - z**2 + 
     -          C33phiq3*(vol2 - 2*vol2*z**2 + 
     -             C33phiq3*vol4*(-1 + z**2))))*
     -     xI2(MT**2 - (s*(1d0+beta*z))/2d0,0d0,MW**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z)))



!          print *, "MCFM slf(4)",slf(4)
!          print *, "MAKU slf(4)",MARKUS_slf4
!          pause 
         
         slf(4) = MARKUS_slf4
         
         

!       slf(5) = 
!      .     (alpha*genfacself*gw_sq*sigma0*
!      .     ((-0.125_dp*(1._dp - beta**2)*s*(-mh**2 
!      .     + (1._dp - beta**2)*s)*(1._dp - beta**4*z**4 
!      .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
!      .     - 3._dp*beta*z**2))*db0(mt**2,mt**2,mh**2))/(1._dp + beta*z) 
!      .     + (0.5_dp*(1._dp - beta**4*z**4 + beta**2*(1._dp 
!      .     - beta**2)*(1._dp - 3._dp*z**2))*(xI1(mh**2,musq,ep) 
!      .     - xI1(mt**2,musq,ep)))/(1._dp + beta**2 + 2._dp*beta*z) 
!      .     + (0.25_dp*((1._dp - beta**2)**2*s*(1._dp + beta**2 
!      .     - 2._dp*beta**4 - beta**2*(2._dp - 3._dp*beta**2)*z**2 
!      .     - beta**4*z**4) + mh**2*(-2._dp - 2._dp*beta**2 
!      .     + 5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp 
!      .     - 2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 
!      .     + beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 
!      .     + beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5))*
!      .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp + beta*z)**2 - 
!      .     (0.125_dp*(1._dp - beta**2)*(s*(1._dp + beta**2 
!      .     - 5._dp*beta**4 - 2._dp*beta**6 + 4._dp*beta**8 
!      .     + 2._dp*beta*(1._dp - beta**2 - 5._dp*beta**4 
!      .     + 4._dp*beta**6)*z - beta**2*(2._dp - 2._dp*beta**2 
!      .     - 5._dp*beta**4 + 6._dp*beta**6)*z**2 - 2._dp*beta**3*(1._dp 
!      .     - 7._dp*beta**2 + 6._dp*beta**4)*z**3 + beta**4*(2._dp 
!      .     - 3._dp*beta**2 + 2._dp*beta**4)*z**4 - 2._dp*beta**5*(1._dp 
!      .     - 2._dp*beta**2)*z**5 + beta**6*z**6) - 2._dp*beta**2*mh**2*
!      .     (1._dp - z**2)*(2._dp + beta**2 - 2._dp*beta**4 + 
!      .     beta*(3._dp - 2._dp*beta**2)*z + beta**3*z**2*(beta + z)))*
!      .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta*z)**2*
!      .     (1._dp + beta**2 + 2._dp*beta*z))))/(mw**2*pi)


 !MARKUS: replaced gw_sq coupl. by vev     
      slf(5) = 
     .     (genfacself/vev**2*sigma0*
     .     ((-0.125_dp*(1._dp - beta**2)*s*(-mh**2 
     .     + (1._dp - beta**2)*s)*(1._dp - beta**4*z**4 
     .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mt**2,mh**2))/(1._dp + beta*z) 
     .     + (0.5_dp*(1._dp - beta**4*z**4 + beta**2*(1._dp 
     .     - beta**2)*(1._dp - 3._dp*z**2))*(xI1(mh**2,musq,ep) 
     .     - xI1(mt**2,musq,ep)))/(1._dp + beta**2 + 2._dp*beta*z) 
     .     + (0.25_dp*((1._dp - beta**2)**2*s*(1._dp + beta**2 
     .     - 2._dp*beta**4 - beta**2*(2._dp - 3._dp*beta**2)*z**2 
     .     - beta**4*z**4) + mh**2*(-2._dp - 2._dp*beta**2 
     .     + 5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp 
     .     - 2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 
     .     + beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5))*
     .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp + beta*z)**2 - 
     .     (0.125_dp*(1._dp - beta**2)*(s*(1._dp + beta**2 
     .     - 5._dp*beta**4 - 2._dp*beta**6 + 4._dp*beta**8 
     .     + 2._dp*beta*(1._dp - beta**2 - 5._dp*beta**4 
     .     + 4._dp*beta**6)*z - beta**2*(2._dp - 2._dp*beta**2 
     .     - 5._dp*beta**4 + 6._dp*beta**6)*z**2 - 2._dp*beta**3*(1._dp 
     .     - 7._dp*beta**2 + 6._dp*beta**4)*z**3 + beta**4*(2._dp 
     .     - 3._dp*beta**2 + 2._dp*beta**4)*z**4 - 2._dp*beta**5*(1._dp 
     .     - 2._dp*beta**2)*z**5 + beta**6*z**6) - 2._dp*beta**2*mh**2*
     .     (1._dp - z**2)*(2._dp + beta**2 - 2._dp*beta**4 + 
     .     beta*(3._dp - 2._dp*beta**2)*z + beta**3*z**2*(beta + z)))*
     .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta*z)**2*
     .     (1._dp + beta**2 + 2._dp*beta*z))))/(8d0*pi**2)

     
      vrt(1) = 
     .     (0.125_dp*alpha*genfacvert*sigma0*
     .     (4._dp*(1._dp - beta**2)*(gat_sq + gvt_sq) 
     .     - (4._dp*(2._dp*(gat_sq + gvt_sq)*mz**2 + 
     .     (1._dp - beta**2)*(-3._dp*gat_sq + gvt_sq)*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mt**2,mz**2))/
     .     (1._dp - beta*z) + (16._dp*(gat_sq + gvt_sq)*(1._dp 
     .     + 2._dp*beta**2 - beta**4 + beta*(1._dp + 4._dp*beta**2 
     .     - 3._dp*beta**4)*z - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(xI1(mt**2,musq,ep) 
     .     - xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + (4._dp*((2._dp*(gat_sq + gvt_sq)*mz**2*(1._dp + 8._dp*beta**2 
     .     - 11._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 
     .     - 3._dp*beta**4)*z - 2._dp*beta**2*(5._dp - 8._dp*beta**2 
     .     + 3._dp*beta**4)*z**2 - 6._dp*beta**3*(1._dp - beta**2)*
     .     z**3 - 2._dp*beta**4*(2._dp - beta**2)*z**4 
     .     - 2._dp*beta**5*z**5))/s + (1._dp - beta**2)**2*
     .     (gvt_sq*(1._dp + beta**2 - 4._dp*beta**4 
     .     - beta*(1._dp - beta**2)*z - 2._dp*beta**2*(1._dp 
     .     - 3._dp*beta**2)*z**2 - 2._dp*beta**4*z**4) + gat_sq*(1._dp 
     .     - 7._dp*beta**2 + 12._dp*beta**4 - beta*(1._dp 
     .     - beta**2)*z + 6._dp*beta**2*(1._dp - 3._dp*beta**2)*z**2 
     .     + 6._dp*beta**4*z**4)))*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     ((1._dp - beta**2)*(1._dp - beta**2*z**2)) + 
     .     (4._dp*((2._dp*(gat_sq + gvt_sq)*mz**2*(1._dp - 4._dp*beta**2 
     .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
     .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp - beta**2)*
     .     z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5))/s 
     .     - gvt_sq*(2._dp + 3._dp*beta**2 - 6._dp*beta**4 
     .     - beta**6 + 4._dp*beta**8 + beta*(3._dp + 7._dp*beta**2 
     .     - 13._dp*beta**4 + 7._dp*beta**6)*z - beta**2*(7._dp 
     .     - 18._dp*beta**2 + 3._dp*beta**4 + 6._dp*beta**6)*z**2 
     .     - beta**3*(14._dp - 26._dp*beta**2 + 12._dp*beta**4)*z**3 
     .     - beta**4*(10._dp - 6._dp*beta**2 - 2._dp*beta**4)*z**4 
     .     - beta**5*(8._dp - 4._dp*beta**2)*z**5 
     .     - 2._dp*beta**6*z**6) + gat_sq*(-2._dp + 5._dp*beta**2 
     .     - 10._dp*beta**4 - 7._dp*beta**6 + 12._dp*beta**8 
     .     - beta*(3._dp - 9._dp*beta**2 + 35._dp*beta**4 
     .     - 25._dp*beta**6)*z - beta**2*(1._dp - 6._dp*beta**2 
     .     - 11._dp*beta**4 + 18._dp*beta**6)*z**2 - 2._dp*beta**3*(1._dp 
     .     - 19._dp*beta**2 + 18._dp*beta**4)*z**3 + 2._dp*beta**4*(1._dp 
     .     - 3._dp*beta**2 + 3._dp*beta**4)*z**4 - 4._dp*beta**5*(2._dp 
     .     - 3._dp*beta**2)*z**5 + 2._dp*beta**6*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) + 2._dp*(1._dp - beta**2)**2*
     .     (gat_sq + gvt_sq)*s*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep)))/pi



      vrt(2) = 
     .     (alpha*genfacvert*gw_sq*sigma0*
     .     (1._dp - beta**2 + (0.5_dp*(4._dp*(mb**2 - mw**2) + 
     .     (1._dp - beta**2)*s)*(1._dp - beta**4*z**4 + 
     .     beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mb**2,mw**2))/
     .     (1._dp - beta*z) + (4._dp*(1._dp + 2._dp*beta**2 
     .     - beta**4 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(xI1(mb**2,musq,ep) 
     .     - xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.5_dp*((1._dp - beta**2)*s*(3._dp + 3._dp*beta**4 
     .     - 4._dp*beta**6 - beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z - 6._dp*beta**2*(1._dp - beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**6*z**4 
     .     - 2._dp*beta**5*z**5) + 4._dp*(-mb**2 + mw**2)*(1._dp 
     .     + 8._dp*beta**2 - 11._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**2*(5._dp 
     .     - 8._dp*beta**2 + 3._dp*beta**4)*z**2 - 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 - 2._dp*beta**4*(2._dp 
     .     - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mb**2,mw**2,musq,ep))/((1._dp - beta**2)*s*
     .     (1._dp - beta**2*z**2)) + (0.5_dp*(-s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)*(3._dp + 3._dp*beta**4 - 4._dp*beta**6 
     .     - beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     - 6._dp*beta**2*(1._dp - beta**4)*z**2 - 6._dp*beta**3*
     .     (1._dp - beta**2)*z**3 - 2._dp*beta**6*z**4 
     .     - 2._dp*beta**5*z**5) + 4._dp*(-mb**2 + mw**2)*
     .     (1._dp - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 
     .     + beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     + 2._dp*beta**2*(1._dp + 2._dp*beta**2 - 3._dp*beta**4)*z**2 
     .     + 6._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(1._dp - beta**2*z**2)) 
     .     + 2._dp*(1._dp - beta**2)*mb**2*
     .     xI3(0._dp,mt**2,t,mb**2,mb**2,mw**2,musq,ep)))/pi




!       vrt(3) = 
!      .     (2._dp*alpha*chifac*genfacvert*sigma0*
!      .     (0.125_dp*(1._dp - beta**2)**2*s - (0.25_dp*(1._dp 
!      .     - beta**2)*mz**2*s*(1._dp - beta**4*z**4 
!      .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
!      .     - 3._dp*beta*z**2))*db0(mt**2,mt**2,mz**2))/
!      .     (1._dp - beta*z) + (0.5_dp*(1._dp + 2._dp*beta**2 
!      .     - beta**4 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
!      .     - 2._dp*beta**4*z**4*(1._dp + beta*z) - 2._dp*(1._dp 
!      .     - beta**2)*(2._dp*beta**2*z**2 + 3._dp*beta**3*z**3))*
!      .     (xI1(mt**2,musq,ep) - xI1(mz**2,musq,ep)))/
!      .     ((1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) 
!      .     + (0.125_dp*(-(1._dp - beta**2)**2*s*(1._dp 
!      .     - beta*z)*(1._dp + beta**2 + 2._dp*beta*z) 
!      .     + 2._dp*mz**2*(1._dp + 8._dp*beta**2 - 11._dp*beta**4 
!      .     + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
!      .     - 2._dp*beta**2*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 
!      .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**4*
!      .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5))*
!      .     xI2(mt**2,mt**2,mz**2,musq,ep))/(1._dp - beta**2*z**2) 
!      .     + (0.125_dp*(1._dp - beta**2)*(2._dp*mz**2*(1._dp 
!      .     - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
!      .     - 8._dp*beta**2 + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp 
!      .     + 2._dp*beta**2 - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
!      .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
!      .     + beta*s*(-beta*(1._dp + beta**4) + (1._dp - 7._dp*beta**2 
!      .     + beta**4 + beta**6)*z + beta*(3._dp - 12._dp*beta**2 
!      .     + 7._dp*beta**4)*z**2 + 6._dp*beta**2*(1._dp - beta**2)*
!      .     z**3 + 2._dp*beta**3*(4._dp - 3._dp*beta**2)*z**4 
!      .     + 4._dp*beta**4*z**5 + 2._dp*beta**5*z**6))*
!      .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
!      .     (1._dp - beta**2*z**2)) 
!      .     - 0.0625_dp*(1._dp - beta**2)**2*s**2*
!      .     (1._dp + beta**2 + 2._dp*beta*z)*
!      .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep)))/pi


     
     
       MARKUS_vrt3 = (alpha)/(64*MW**2*Pi*SW2)/16d0 * ( 
     -    (16*s*(-7 + 9*beta*z)*
     -     (-1 - (2*C33phiq3 + C33phiu)*
     -        (8*vol2 + (2*C33phiq3 + C33phiu)*vol4) + 
     -       beta**5*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)*z + 
     -       beta*(1 - 3*(2*C33phiq3 + C33phiu)*vol2)*z + 
     -       beta**4*(-1 - (2*C33phiq3 + C33phiu)*vol2 - 
     -          (3*(2*C33phiq3 + C33phiu)*vol2 + 
     -             2*(4*C33phiq3**2 + 2*C33phiq3*C33phiu + 
     -                C33phiu**2)*vol4)*z**2 + 
     -          (4*C33phiq3**2 + C33phiu**2)*vol4*z**4) + 
     -       2*beta**3*z*(-1 + 2*C33phiq3*vol2 - 
     -          4*C33phiq3**2*vol4*(-1 + z**2) + 
     -          C33phiu*(vol2 + C33phiu*vol4 - C33phiu*vol4*z**2)) + 
     -       beta**2*(2 - 8*C33phiq3**2*vol4*(-2 + z**2) + 
     -          2*C33phiq3*vol2*
     -           (9 + 2*C33phiu*vol2 + (3 + 2*C33phiu*vol2)*z**2) + 
     -          C33phiu*(-2*C33phiu*vol4*(-2 + z**2) + 
     -             3*vol2*(3 + z**2)))))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**2) + 
     -  (32*(-1 + beta**2)*s*(MZ + (2*C33phiq3 + C33phiu)*MZ*vol2)**2*
     -     (-7 + 9*beta*z)*(-1 + 
     -       beta*(z + beta*(-2 - beta*z + 3*z**2 + 
     -             beta**2*(2 - 3*z**2 + z**4))))*
     -     DB0(MT**2,MT**2,MZ**2))/(3.*(-1 + beta*z)*(1 + beta*z)**2)
     -   + (64*(-7 + 9*beta*z)*
     -     ((-((4*C33phiq3**2 + C33phiu**2)*vol4) - 
     -          beta*(4*C33phiq3**2 + C33phiu**2)*vol4*z + 
     -          beta**3*(4*C33phiq3**2 + C33phiu**2)*vol4*z*
     -           (-1 + 3*z**2) + 
     -          beta**2*(1 + 2*C33phiq3*vol2 + C33phiu*vol2 - 
     -             4*C33phiq3**2*vol4 - C33phiu**2*vol4 + 
     -             (-1 - (2*C33phiq3 + C33phiu)*vol2 + 
     -                3*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2) + 
     -          beta**5*(4*C33phiq3**2 + C33phiu**2)*vol4*z*
     -           (1 - 3*z**2 + z**4) - 
     -          beta**6*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)*
     -           (2 - 3*z**2 + z**4) + 
     -          beta**4*((-1 + z**2)**2 + 
     -             2*C33phiq3*vol2*(-1 + z**2)**2 + 
     -             C33phiu*vol2*(-1 + z**2)**2 + 
     -             4*C33phiq3**2*vol4*(1 - 3*z**2 + z**4) + 
     -             C33phiu**2*vol4*(1 - 3*z**2 + z**4)))/
     -        (1 + beta**2 + 2*beta*z) + 
     -       (1 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -          (4*C33phiq3**2 + 8*C33phiq3*C33phiu + C33phiu**2)*vol4
     -          )*(-1 + beta*
     -           (z + beta*(-2 - beta*z + 3*z**2 + 
     -                beta**2*(2 - 3*z**2 + z**4)))))*
     -     xI1(MT**2,musq,ep))/(3.*(-1 + beta*z)*(1 + beta*z)**2) - 
     -  (64*(-7 + 9*beta*z)*(-1 - 
     -       vol2*(3*C33phiu + C33phiq3*(6 + 8*C33phiu*vol2)) - 
     -       2*(4*C33phiq3**2 + C33phiu**2)*vol4 - 
     -       beta*(1 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -          2*(2*C33phiq3 + C33phiu)**2*vol4)*z + 
     -       2*beta**2*(-1 - 
     -          4*vol2*(C33phiu + C33phiq3*(2 + 3*C33phiu*vol2)) - 
     -          2*(4*C33phiq3**2 + C33phiu**2)*vol4 + 
     -          (2 + 7*C33phiu*vol2 + 
     -             2*C33phiq3*vol2*(7 + 10*C33phiu*vol2) + 
     -             16*C33phiq3**2*vol4 + 4*C33phiu**2*vol4)*z**2) + 
     -       beta**6*vol2*(4*C33phiq3**2*vol2 + 
     -          C33phiu*(2 + C33phiu*vol2) + 
     -          C33phiq3*(4 + 8*C33phiu*vol2))*(2 - 3*z**2 + z**4) + 
     -       beta**4*(1 + 2*C33phiq3*vol2 + C33phiu*vol2 + 
     -          4*C33phiq3**2*vol4 + C33phiu**2*vol4 - 
     -          (4 + 8*vol2*
     -              (C33phiu + 2*C33phiq3*(1 + C33phiu*vol2)) + 
     -             5*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -          2*(1 + 2*C33phiu*vol2 + 
     -             4*C33phiq3*vol2*(1 + C33phiu*vol2) + 
     -             4*C33phiq3**2*vol4 + C33phiu**2*vol4)*z**4) + 
     -       beta**3*z*(-4 - 12*C33phiu*vol2 - 5*C33phiu**2*vol4 + 
     -          3*(2 + 6*C33phiu*vol2 + 3*C33phiu**2*vol4)*z**2 + 
     -          4*C33phiq3*vol2*(3 + 4*C33phiu*vol2)*(-2 + 3*z**2) + 
     -          4*C33phiq3**2*vol4*(-5 + 9*z**2)) + 
     -       beta**5*z*(3 + 9*C33phiu*vol2 + 4*C33phiu**2*vol4 - 
     -          3*(2 + 6*C33phiu*vol2 + 3*C33phiu**2*vol4)*z**2 + 
     -          (2 + 6*C33phiu*vol2 + 3*C33phiu**2*vol4)*z**4 + 
     -          2*C33phiq3*vol2*(3 + 4*C33phiu*vol2)*
     -           (3 - 6*z**2 + 2*z**4) + 
     -          4*C33phiq3**2*vol4*(4 - 9*z**2 + 3*z**4)))*
     -     xI1(MZ**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z))
     -   + (16*(-7 + 9*beta*z)*
     -     (-2*MZ**2*(1 + 5*(2*C33phiq3 + C33phiu)*vol2 + 
     -          2*(4*C33phiq3**2 + 8*C33phiq3*C33phiu + C33phiu**2)*
     -           vol4 + beta*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)*z + 
     -          2*beta**3*z*(2 + 16*C33phiq3*vol2 + 8*C33phiu*vol2 + 
     -             24*C33phiq3*C33phiu*vol2**2 + 
     -             12*C33phiq3**2*vol4 + 3*C33phiu**2*vol4 - 
     -             3*(1 + 3*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(3 + 4*C33phiu*vol2) + 
     -                4*C33phiq3**2*vol4 + C33phiu**2*vol4)*z**2) + 
     -          2*beta**6*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*
     -           (2 - 3*z**2 + z**4) + 
     -          beta**4*(-11 - 50*C33phiq3*vol2 - 25*C33phiu*vol2 - 
     -             56*C33phiq3*C33phiu*vol2**2 - 
     -             40*C33phiq3**2*vol4 - 10*C33phiu**2*vol4 + 
     -             8*(2 + 5*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(5 + 6*C33phiu*vol2) + 
     -                8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4)*z**2 - 
     -             2*(2 + 5*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(5 + 6*C33phiu*vol2) + 
     -                8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4)*z**4) + 
     -          2*beta**2*(4 + 9*C33phiu*vol2 + 3*C33phiu**2*vol4 - 
     -             (5 + 14*C33phiu*vol2 + 5*C33phiu**2*vol4)*z**2 + 
     -             4*C33phiq3**2*vol4*(3 - 5*z**2) + 
     -             2*C33phiq3*vol2*
     -              (9 + 10*C33phiu*vol2 - 
     -                2*(7 + 9*C33phiu*vol2)*z**2)) + 
     -          beta**5*z*(-3 - 11*C33phiu*vol2 - 4*C33phiu**2*vol4 + 
     -             6*(1 + 3*C33phiu*vol2 + C33phiu**2*vol4)*z**2 - 
     -             2*(1 + 3*C33phiu*vol2 + C33phiu**2*vol4)*z**4 - 
     -             8*C33phiq3**2*vol4*(2 - 3*z**2 + z**4) - 
     -             2*C33phiq3*vol2*
     -              (11 + 6*z**2*(-3 + z**2) + 
     -                8*C33phiu*vol2*(2 - 3*z**2 + z**4)))) + 
     -       (-1 + beta**2)*s*
     -        (-1 - (2*C33phiq3 + C33phiu)*
     -           (5*vol2 + 2*(2*C33phiq3 + C33phiu)*vol4) - 
     -          beta*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)*z + 
     -          2*beta**2*(-4*vol2*
     -              (C33phiu + 2*C33phiq3*(1 + C33phiu*vol2)) - 
     -             2*(4*C33phiq3**2 + C33phiu**2)*vol4 + 
     -             (1 + 9*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(9 + 8*C33phiu*vol2) + 
     -                16*C33phiq3**2*vol4 + 4*C33phiu**2*vol4)*z**2)
     -           + 2*beta**3*z*
     -           (1 + 12*C33phiq3**2*vol4*(-1 + z**2) + 
     -             2*C33phiq3*vol2*
     -              (-5 + 6*z**2 + 6*C33phiu*vol2*(-1 + z**2)) + 
     -             C33phiu*(3*C33phiu*vol4*(-1 + z**2) + 
     -                vol2*(-5 + 6*z**2))) + 
     -          beta**4*(1 - 2*z**2 + 
     -             8*C33phiq3**2*vol4*(2 - 4*z**2 + z**4) + 
     -             2*C33phiq3*vol2*
     -              (9 + 8*C33phiu*vol2 - 
     -                2*(9 + 8*C33phiu*vol2)*z**2 + 
     -                4*(1 + C33phiu*vol2)*z**4) + 
     -             C33phiu*(2*C33phiu*vol4*(2 - 4*z**2 + z**4) + 
     -                vol2*(9 - 18*z**2 + 4*z**4))) + 
     -          beta**5*z*(-1 + 
     -             8*C33phiq3**2*vol4*(2 - 3*z**2 + z**4) + 
     -             2*C33phiq3*vol2*
     -              (7 + 4*z**2*(-3 + z**2) + 
     -                4*C33phiu*vol2*(2 - 3*z**2 + z**4)) + 
     -             C33phiu*(2*C33phiu*vol4*(2 - 3*z**2 + z**4) + 
     -                vol2*(7 + 4*z**2*(-3 + z**2))))))*
     -     xI2(MT**2,MT**2,MZ**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**3) + 
     -  (16*(-7 + 9*beta*z)*(beta*s*
     -        (-((-1 + beta)*beta*(1 + beta)*(1 + beta**4)) - 
     -          (beta + beta**3 - 3*beta**5 + beta**7)*
     -           (2*C33phiq3 + C33phiu)*vol2 + 
     -          2*beta*(-2 - beta**2 + 2*beta**4)*
     -           (4*C33phiq3**2 + C33phiu**2)*vol4 + 
     -          (-1 + 8*beta**2 - 8*beta**4 + beta**8 + 
     -             (1 - 6*beta**4 + 4*beta**6 + beta**8)*
     -              (2*C33phiq3 + C33phiu)*vol2 + 
     -             2*beta**2*(-7 + 4*beta**4)*
     -              (4*C33phiq3**2 + C33phiu**2)*vol4)*z + 
     -          beta*(-3 + 15*beta**2 - 19*beta**4 + 7*beta**6 + 
     -             (5 + beta**2 - 15*beta**4 + 9*beta**6)*
     -              (2*C33phiq3 + C33phiu)*vol2 + 
     -             2*(2 - 7*beta**2 + 2*beta**6)*
     -              (4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 - 
     -          2*beta**2*(3*(-1 + beta**2)**2 + 
     -             2*(-1 + beta**4)*(2*C33phiq3 + C33phiu)*vol2 + 
     -             (-7 + 4*(beta**2 + beta**4))*
     -              (4*C33phiq3**2 + C33phiu**2)*vol4)*z**3 - 
     -          2*beta**3*(4 - 7*beta**2 + 3*beta**4 + 
     -             (2 - 5*beta**2 + 3*beta**4)*(2*C33phiq3 + C33phiu)*
     -              vol2 + (-8 + 4*beta**2 + 3*beta**4)*
     -              (4*C33phiq3**2 + C33phiu**2)*vol4)*z**4 - 
     -          2*beta**4*(-2*(-1 + beta**2)*
     -              (1 + 2*C33phiq3*vol2 + C33phiu*vol2) + 
     -             (-4 + beta**2)*(4*C33phiq3**2 + C33phiu**2)*vol4)*
     -           z**5 + 2*beta**5*
     -           ((-1 + beta**2)*
     -              (1 + 2*C33phiq3*vol2 + C33phiu*vol2) + 
     -             (2 + beta**2)*(4*C33phiq3**2 + C33phiu**2)*vol4)*
     -           z**6 + 2*beta**6*(4*C33phiq3**2 + C33phiu**2)*vol4*
     -           z**7) + 2*MZ**2*
     -        (-1 - (2*C33phiq3 + C33phiu)*vol2 - 
     -          2*(4*C33phiq3**2 + C33phiu**2)*vol4 - 
     -          beta*(1 + 2*C33phiq3*vol2 + 16*C33phiq3**2*vol4 + 
     -             C33phiu*(vol2 + 4*C33phiu*vol4))*z + 
     -          beta**3*z*(9 + 26*C33phiq3*vol2 + 13*C33phiu*vol2 + 
     -             16*C33phiq3*C33phiu*vol2**2 + 
     -             2*(-3 - 5*C33phiu*vol2 - 
     -                2*C33phiq3*vol2*(5 + 4*C33phiu*vol2) + 
     -                16*C33phiq3**2*vol4 + 4*C33phiu**2*vol4)*z**2)
     -           + beta**2*(5 + 14*C33phiq3*vol2 + 7*C33phiu*vol2 + 
     -             8*C33phiq3*C33phiu*vol2**2 + 
     -             2*(-1 - 2*vol2*
     -                 (C33phiu + 2*C33phiq3*(1 + C33phiu*vol2)) + 
     -                (4*C33phiq3**2 + C33phiu**2)*vol4)*z**2) + 
     -          2*beta**8*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*
     -           (2 - 3*z**2 + z**4) + 
     -          beta**6*(-7 + 10*z**2 - 2*z**4 + 
     -             C33phiu*vol2*(-9 + 12*z**2 - 2*z**4) + 
     -             2*C33phiq3*vol2*
     -              (-9 - 4*C33phiu*vol2 + 
     -                4*(3 + C33phiu*vol2)*z**2 - 2*z**4) + 
     -             8*C33phiq3**2*vol4*(-1 + 2*z**2 - 3*z**4 + z**6) + 
     -             2*C33phiu**2*vol4*(-1 + 2*z**2 - 3*z**4 + z**6)) + 
     -          beta**5*z*(-13 - 25*C33phiu*vol2 - 
     -             8*C33phiu**2*vol4 + 
     -             4*(3 + 7*C33phiu*vol2 + C33phiu**2*vol4)*z**2 - 
     -             2*(1 + 3*C33phiu*vol2)*z**4 + 
     -             16*C33phiq3**2*vol4*(-2 + z**2) - 
     -             2*C33phiq3*vol2*
     -              (25 - 28*z**2 + 6*z**4 + 
     -                8*C33phiu*vol2*(3 - 4*z**2 + z**4))) + 
     -          beta**4*(-1 - 2*z**2 + 
     -             C33phiu*vol2*(-5 + 4*z**2 - 2*z**4) + 
     -             8*C33phiq3**2*vol4*(-1 - z**2 + 3*z**4) + 
     -             2*C33phiu**2*vol4*(-1 - z**2 + 3*z**4) - 
     -             2*C33phiq3*vol2*
     -              (5 - 4*z**2 + 2*z**4 + 
     -                4*C33phiu*vol2*(2 - 3*z**2 + z**4))) + 
     -          beta**7*z*(5 - 6*z**2 + 2*z**4 + 
     -             16*C33phiq3**2*vol4*(2 - 3*z**2 + z**4) + 
     -             4*C33phiu**2*vol4*(2 - 3*z**2 + z**4) + 
     -             C33phiu*vol2*(13 + 6*z**2*(-3 + z**2)) + 
     -             2*C33phiq3*vol2*
     -              (13 + 6*z**2*(-3 + z**2) + 
     -                8*C33phiu*vol2*(2 - 3*z**2 + z**4)))))*
     -     xI2(MT**2 - (s*(1d0 + beta*z))/2d0,MT**2,MZ**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z))
     -   - (8*(-1 + beta**2)**2*s**2*
     -     (1 + 2*C33phiq3*vol2 + C33phiu*vol2)*
     -     (1 + beta**2 + 2*beta*z)*(-7 + 9*beta*z)*
     -  xI3(0d0,MT**2,MT**2 - (s*(1d0 + beta*z))/2d0,MT**2,MT**2,MZ**2,
     -      musq,ep))/(3.*(1 + beta*z)**2)   )

!          print *, "MCFM vrt(3)",vrt(3)
!          print *, "MAKU vrt(3)",MARKUS_vrt3   
!          pause 
         
         vrt(3) = MARKUS_vrt3
     
     


!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
!       vrt(4) = 
!      .     (alpha*genfacvert*(0.125_dp/sw2)*sigma0*
!      .     (0.125_dp*(1._dp - beta**2)*(4._dp*mb**2 + (1._dp 
!      .     - beta**2)*s) + (0.0625_dp*(16._dp*mb**2*(mb**2 
!      .     - mw**2) - 4._dp*(1._dp - beta**2)*(2._dp*mb**2 
!      .     + mw**2)*s + (1._dp - beta**2)**2*s**2)*(1._dp 
!      .     - beta**4*z**4 + beta*(1._dp - beta**2)*
!      .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mb**2,
!      .     mw**2))/(1._dp - beta*z) + (0.5_dp*(4._dp*mb**2 + 
!      .     (1._dp - beta**2)*s)*(1._dp + 2._dp*beta**2 - beta**4 
!      .     + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**4
!      .     *z**4*(1._dp + beta*z) - 2._dp*(1._dp - beta**2)*
!      .     (2._dp*beta**2*z**2 + 3._dp*beta**3*z**3))*(xI1(mb**2,musq,ep) 
!      .     - xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
!      .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + 
!      .     (0.0625_dp*(-8._dp*(1._dp - beta**2)**2*mb**2*s*(1._dp 
!      .     - beta**2 + 4._dp*beta**4 + beta*(1._dp - beta**2)*z 
!      .     - 6._dp*beta**4*z**2 + 2._dp*beta**4*z**4) + (-16._dp*mb**4 
!      .     + 4._dp*mw**2*(4._dp*mb**2 + (1._dp - beta**2)*s))*
!      .     (1._dp + 8._dp*beta**2 - 11._dp*beta**4 + 4._dp*beta**6 + 
!      .     beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**2*
!      .     (5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 - 6._dp*beta**3*
!      .     (1._dp - beta**2)*z**3 - 2._dp*beta**4*
!      .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5) 
!      .     - (1._dp - beta**2)**2*s**2*(1._dp - 4._dp*beta**2 
!      .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
!      .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
!      .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
!      .     - beta**2)*z**3 + 2._dp*beta**6*z**4 
!      .     + 2._dp*beta**5*z**5))*xI2(mt**2,mb**2,mw**2,musq,ep))/
!      .     ((1._dp - beta**2)*s*(1._dp - beta**2*z**2)) 
!      .     + (0.0625_dp*((-16._dp*mb**4 + 4._dp*mw**2*(4._dp*mb**2 
!      .     + (1._dp - beta**2)*s))*(1._dp - 4._dp*beta**2 
!      .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
!      .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
!      .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
!      .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
!      .     + (1._dp - beta**2)*s**2*(1._dp + beta**2 + 2._dp*beta*z)*
!      .     (1._dp - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 
!      .     + beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
!      .     + 2._dp*beta**2*(1._dp + 2._dp*beta**2 - 3._dp*beta**4)*z**2 
!      .     + 6._dp*beta**3*(1._dp - beta**2)*z**3 + 2._dp*beta**6*z**4 
!      .     + 2._dp*beta**5*z**5) + 8._dp*beta*mb**2*s*(-beta*(3._dp 
!      .     - 4._dp*beta**2 - beta**4 + 4._dp*beta**6) + (1._dp 
!      .     - 11._dp*beta**2 + 13._dp*beta**4 - 7._dp*beta**6)*z 
!      .     + beta*(5._dp - 18._dp*beta**2 + 5._dp*beta**4 
!      .     + 6._dp*beta**6)*z**2 + 2._dp*beta**2*(5._dp - 11._dp*beta**2 
!      .     + 6._dp*beta**4)*z**3 + 2._dp*beta**3*(5._dp - 3._dp*beta**2 
!      .     - beta**4)*z**4 + 4._dp*beta**4*(2._dp 
!      .     - beta**2)*z**5 + 2._dp*beta**5*z**6))*
!      .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta**2 
!      .     + 2._dp*beta*z)*(1._dp - beta**2*z**2)) 
!      .     - 0.25_dp*(1._dp - beta**2)*mb**2*(-4._dp*mb**2 + 3._dp*s 
!      .     + beta**2*s + 4._dp*beta*s*z)*
!      .     xI3(0._dp,mt**2,t,mb**2,mb**2,mw**2,musq,ep)))/(mw**2*pi)

           

     
       MARKUS_vrt4 = (alpha)/(64*MW**2*Pi*SW2)/16d0 * ( 
     -   (-16*s*(-7 + 9*beta*z)*(1 + 14*C33phiq3*vol2 + C33phiq3**2*vol4 + 
     -beta**5*(-1 + C33phiq3**2*vol4)*z + 
     -beta*(-1 + 8*C33phiq3*vol2 + C33phiq3**2*vol4)*z + 
     -2*beta**3*z*(1 - 4*C33phiq3*vol2 + C33phiq3**2*vol4*(-3 + 2*z**2)) + 
     -2*beta**2*(-1 + C33phiq3**2*vol4*(-3 + 2*z**2) - 
     -C33phiq3*vol2*(7 + 3*z**2)) + 
     -beta**4*(1 + 6*C33phiq3*vol2*z**2 + 
     -C33phiq3**2*vol4*(-1 + 4*z**2 - 2*z**4))))/
     -(3.*(-1 + beta*z)*(1 + beta*z)**2) + 
     -(8*(-1 + beta**2)*s*(4*MW**2 + (-1+beta**2)*s)*(1 + C33phiq3*vol2)**2*
     -(-7 + 9*beta*z)*(-1 + beta*
     -(z + beta*(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4))))*
     -DB0(MT**2,MB**2,MW**2))/(3.*(-1 + beta*z)*(1 + beta*z)**2) + 
     -(64*(1 + 4*C33phiq3*vol2 + C33phiq3**2*vol4)*(-7 + 9*beta*z)*
     -(-1 + beta*(z + beta*(-2 - beta*z + 3*z**2 + 
     -beta**2*(2 - 3*z**2 + z**4))))*0d0)/
     -(3.*(-1 + beta*z)*(1 + beta*z)**2) + 
     -(64*(-7 + 9*beta*z)*((beta**2*(-1 + z**2)*
     -(1 + 3*C33phiq3**2*vol4 + 6*beta*C33phiq3**2*vol4*z + 
     -beta**4*(1 + C33phiq3**2*vol4)*(-2 + z**2) + 
     -2*beta**3*C33phiq3**2*vol4*z*(-2 + z**2) + 
     -beta**2*(1 - z**2 + C33phiq3**2*vol4*(1 + z**2))))/
     -(1 + beta**2 + 2*beta*z) - 
     -(1 + 4*C33phiq3*vol2 + 3*C33phiq3**2*vol4)*
     -(-1 + beta*(z + beta*
     -(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4)))))*
     -xI1(MW**2,musq,ep))/(3.*(-1 + beta*z)*(1 + beta*z)**2) + 
     -(8*(-7 + 9*beta*z)*(2*(-((-1 + beta**2)*s*(1 + C33phiq3**2*vol4)) +
     -4*MW**2*(1 + 4*C33phiq3*vol2 + C33phiq3**2*vol4))*(1 + beta*z)*
     -(-1 + beta*(z + beta*
     -(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4)))) + 
     -(1 + C33phiq3*vol2)*(4*MW**2*
     -(1 - C33phiq3*vol2 + 
     -beta*(-4*beta - 4*beta**5*(1 + C33phiq3*vol2) + 
     - beta**3*(7 + 5*C33phiq3*vol2) + 
     - (-1 + beta**2)**2*(-1 + C33phiq3*vol2)*z + 
     - 2*beta*(1 - 4*beta**2 + 3*beta**4)*(1 + C33phiq3*vol2)*z**2 - 
     - 2*beta**3*(-1 + beta**2)*(1 + C33phiq3*vol2)*z**4)) + 
     -(-1 + beta**2)**2*s*
     -(3 + C33phiq3*vol2 + 
     -beta*(z - C33phiq3*vol2*z + beta**2*(-1 + C33phiq3*vol2)*z - 
     - beta*(3 + C33phiq3*vol2)*(-1 + 2*z**2) - 
     - 2*beta**3*(1 + C33phiq3*vol2)*(2 - 3*z**2 + z**4)))))*
     -xI2(MT**2,0d0,MW**2,musq,ep))/(3.*(-1 + beta*z)*(1 + beta*z)**3) + 
     -(8*(-7 + 9*beta*z)*(4*MW**2*
     -(-1 - 3*C33phiq3**2*vol4 - beta*(z + 7*C33phiq3**2*vol4*z) + 
     -2*beta**8*(1 + C33phiq3*vol2)**2*(2 - 3*z**2 + z**4) + 
     -beta**2*(5 - 2*z**2 - 4*C33phiq3*vol2*(-1 + z**2) + 
     -C33phiq3**2*vol4*(-5 + 6*z**2)) + 
     -beta**3*z*(9 - 6*z**2 - 8*C33phiq3*vol2*(-1 + z**2) + 
     -C33phiq3**2*vol4*(-9 + 22*z**2)) + 
     -beta**5*z*(-13 + 12*z**2 - 2*z**4 - 
     -8*C33phiq3*vol2*(3 - 4*z**2 + z**4) + 
     -C33phiq3**2*vol4*(-3 - 4*z**2 + 2*z**4)) + 
     -beta**4*(-1 - 2*z**2 - 4*C33phiq3*vol2*(2 - 3*z**2 + z**4) + 
     -C33phiq3**2*vol4*(-3 - 2*z**2 + 12*z**4)) + 
     -beta**6*(-7 + 10*z**2 - 2*z**4 + 4*C33phiq3*vol2*(-1 + z**2) + 
     -C33phiq3**2*vol4*(3 - 2*z**2 - 10*z**4 + 4*z**6)) + 
     -beta**7*z*(5 - 6*z**2 + 2*z**4 + 
     -8*C33phiq3*vol2*(2 - 3*z**2 + z**4) + 
     -C33phiq3**2*vol4*(11 + 6*z**2*(-3 + z**2)))) + 
     -s*(1 + beta**2 + 2*beta*z)*
     -(-1 + 4*C33phiq3*vol2 + C33phiq3**2*vol4 + 
     -beta*(-1 + 8*C33phiq3*vol2 + C33phiq3**2*vol4)*z + 
     -beta**2*(5 + 4*C33phiq3*vol2 - 5*C33phiq3**2*vol4 + 
     -2*(-1 - 4*C33phiq3*vol2 + C33phiq3**2*vol4)*z**2) + 
     -beta**3*z*(9 - 17*C33phiq3**2*vol4 + 
     -2*(-3 - 8*C33phiq3*vol2 + 7*C33phiq3**2*vol4)*z**2) + 
     -2*beta**8*(1 + C33phiq3*vol2)**2*(2 - 3*z**2 + z**4) + 
     -beta**5*z*(-13 + 12*z**2 - 2*z**4 - 
     -8*C33phiq3*vol2*(3 - 5*z**2 + z**4) + 
     -C33phiq3**2*vol4*(-3 + 4*z**2 + 2*z**4)) + 
     -beta**4*(-1 - 2*z**2 - 4*C33phiq3*vol2*(3 - 3*z**2 + z**4) + 
     -C33phiq3**2*vol4*(-7 - 2*z**2 + 12*z**4)) + 
     -beta**6*(-7 + 10*z**2 - 2*z**4 + 4*C33phiq3*vol2*(-1+2*z**2)+
     -C33phiq3**2*vol4*(3 + 2*z**2 - 10*z**4 + 4*z**6)) + 
     -beta**7*z*(5 - 6*z**2 + 2*z**4 + 
     -8*C33phiq3*vol2*(2 - 3*z**2 + z**4) + 
     -C33phiq3**2*vol4*(11 + 6*z**2*(-3 + z**2)))))*
     -xI2(MT**2 - (s*(1d0 + beta*z))/2d0,0d0,MW**2,musq,ep))/
     -(3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z)) )

!          print *, "MCFM vrt(4)",vrt(4)
!          print *, "MAKU vrt(4)",MARKUS_vrt4
!          pause 
         
         vrt(4) = MARKUS_vrt4
     
     



!       vrt(5) = 
!      .     (alpha*genfacvert*gw_sq*sigma0*
!      .     (0.125_dp*(1._dp - beta**2)**2*s + (0.25_dp*(1._dp 
!      .     - beta**2)*s*(-mh**2 + (1._dp - beta**2)*s)*
!      .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
!      .     (2._dp*beta - z - 3._dp*beta*z**2))*
!      .     db0(mt**2,mt**2,mh**2))/(1._dp - beta*z) 
!      .     + (0.5_dp*(1._dp + 2._dp*beta**2 - beta**4 
!      .     + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
!      .     - 2._dp*beta**4*z**4*(1._dp + beta*z) 
!      .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
!      .     + 3._dp*beta**3*z**3))*(-xI1(mh**2,musq,ep) 
!      .     + xI1(mt**2,musq,ep)))/((1._dp - beta*z)*
!      .     (1._dp + beta**2 + 2._dp*beta*z)) + (0.125_dp*(-(1._dp 
!      .     - beta**2)**2*s*(1._dp + 5._dp*beta**2 - 8._dp*beta**4 
!      .     + beta*(1._dp - beta**2)*z - 6._dp*beta**2*(1._dp 
!      .     - 2._dp*beta**2)*z**2 - 4._dp*beta**4*z**4) 
!      .     + 2._dp*mh**2*(1._dp + 8._dp*beta**2 - 11._dp*beta**4 
!      .     + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
!      .     - 2._dp*beta**2*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 
!      .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**4*
!      .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5))*
!      .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp - beta**2*z**2) 
!      .     + (0.125_dp*(1._dp - beta**2)*(2._dp*mh**2*(1._dp 
!      .     - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
!      .     - 8._dp*beta**2 + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp 
!      .     + 2._dp*beta**2 - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
!      .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
!      .     + beta*s*(beta*(3._dp - 8._dp*beta**2 - 5._dp*beta**4 
!      .     + 8._dp*beta**6) + (1._dp + beta**2 - 23._dp*beta**4 
!      .     + 17._dp*beta**6)*z - beta*(1._dp - 11._dp*beta**4 
!      .     + 12._dp*beta**6)*z**2 - 2._dp*beta**2*(1._dp - 13._dp*beta**2 
!      .     + 12._dp*beta**4)*z**3 + 2._dp*beta**3*(2._dp - 3._dp*beta**2 
!      .     + 2._dp*beta**4)*z**4 - 4._dp*beta**4*(1._dp 
!      .     - 2._dp*beta**2)*z**5 + 2._dp*beta**5*z**6))*
!      .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
!      .     (1._dp - beta**2*z**2)) 
!      .     + 0.0625_dp*(1._dp - beta**2)**2*s**2*(3._dp 
!      .     - beta**2 + 2._dp*beta*z)*
!      .     xI3(0._dp,mt**2,t,mt**2,mt**2,mh**2,musq,ep)))/(mw**2*pi)

     
 !MARKUS: replaced gw_sq coupl. by vev     
      vrt(5) = 
     .     (genfacvert/vev**2*sigma0*
     .     (0.125_dp*(1._dp - beta**2)**2*s + (0.25_dp*(1._dp 
     .     - beta**2)*s*(-mh**2 + (1._dp - beta**2)*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*
     .     db0(mt**2,mt**2,mh**2))/(1._dp - beta*z) 
     .     + (0.5_dp*(1._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(-xI1(mh**2,musq,ep) 
     .     + xI1(mt**2,musq,ep)))/((1._dp - beta*z)*
     .     (1._dp + beta**2 + 2._dp*beta*z)) + (0.125_dp*(-(1._dp 
     .     - beta**2)**2*s*(1._dp + 5._dp*beta**2 - 8._dp*beta**4 
     .     + beta*(1._dp - beta**2)*z - 6._dp*beta**2*(1._dp 
     .     - 2._dp*beta**2)*z**2 - 4._dp*beta**4*z**4) 
     .     + 2._dp*mh**2*(1._dp + 8._dp*beta**2 - 11._dp*beta**4 
     .     + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**2*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**4*
     .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp - beta**2*z**2) 
     .     + (0.125_dp*(1._dp - beta**2)*(2._dp*mh**2*(1._dp 
     .     - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     - 8._dp*beta**2 + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp 
     .     + 2._dp*beta**2 - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
     .     + beta*s*(beta*(3._dp - 8._dp*beta**2 - 5._dp*beta**4 
     .     + 8._dp*beta**6) + (1._dp + beta**2 - 23._dp*beta**4 
     .     + 17._dp*beta**6)*z - beta*(1._dp - 11._dp*beta**4 
     .     + 12._dp*beta**6)*z**2 - 2._dp*beta**2*(1._dp - 13._dp*beta**2 
     .     + 12._dp*beta**4)*z**3 + 2._dp*beta**3*(2._dp - 3._dp*beta**2 
     .     + 2._dp*beta**4)*z**4 - 4._dp*beta**4*(1._dp 
     .     - 2._dp*beta**2)*z**5 + 2._dp*beta**5*z**6))*
     .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) 
     .     + 0.0625_dp*(1._dp - beta**2)**2*s**2*(3._dp 
     .     - beta**2 + 2._dp*beta*z)*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mh**2,musq,ep)))/(8*pi**2)     
     
     
     
      bx(1) = 
     .     (0.125_dp*alpha*genfacbox*sigma0*
     .     (-2._dp*beta*(gat_sq + gvt_sq)*z*(1._dp - z**2) + 
     .     (8._dp*beta*(gat_sq + gvt_sq)*(1._dp - z**2)*
     .     (4._dp*beta - z - 2._dp*beta*z**2)*(xI1(mt**2,musq,ep) 
     .     - xI1(mz**2,musq,ep)))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)) + (4._dp*((-(gat_sq + gvt_sq)*mz**2*
     .     (-6._dp*beta**3 - (5._dp - 8._dp*beta**2 + 4._dp*beta**4)*z 
     .     - beta*(5._dp - 12._dp*beta**2)*z**2 + (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 + beta*(3._dp 
     .     - 4._dp*beta**2)*z**4))/(beta*s) - gat_sq*(3._dp*beta**2 
     .     - 5._dp*beta**4 - 3._dp*beta*(3._dp - 5._dp*beta**2 
     .     + 2._dp*beta**4)*z - (1._dp + 8._dp*beta**2 
     .     - 9._dp*beta**4)*z**2 + beta*(4._dp - 7._dp*beta**2 
     .     + 3._dp*beta**4)*z**3 + beta**2*(5._dp - 3._dp*beta**2)*z**4) 
     .     + gvt_sq*(5._dp*beta**2 - 3._dp*beta**4 + beta*(1._dp + beta**2 
     .     - 2._dp*beta**4)*z + (1._dp - 4._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - beta**3*(1._dp - beta**2)*z**3 - beta**2*
     .     (1._dp + beta**2)*z**4))*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     (1._dp + beta*z) - 2._dp*((2._dp*(gat_sq + gvt_sq)*mz**2*z*(5._dp 
     .     - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2))/(beta*s) 
     .     - gat_sq*(2._dp - 3._dp*beta*(5._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 + 3._dp*beta*(3._dp - 2._dp*beta**2)*z**3) 
     .     - gvt_sq*(2._dp + beta*(1._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 + beta*(1._dp + 2._dp*beta**2)*z**3))*
     .     xI2(s,mt**2,mt**2,musq,ep) + (4._dp*((2._dp*beta*(gat_sq 
     .     + gvt_sq)*mz**2*(beta + z)*(1._dp - 3._dp*beta**2 + beta*(1._dp 
     .     - 2._dp*beta**2)*z + 2._dp*beta**2*z**2 + beta**3*z**3))/s 
     .     - gvt_sq*(1._dp + 8._dp*beta**2 - 3._dp*beta**6 
     .     + beta*(4._dp + 14._dp*beta**2 - 11._dp*beta**4 
     .     - 2._dp*beta**6)*z - beta**2*(2._dp + beta**2 
     .     + 3._dp*beta**4)*z**2 - beta**3*(11._dp - 6._dp*beta**2 
     .     - beta**4)*z**3 - 2._dp*beta**4*(1._dp 
     .     - beta**2)*z**4 - beta**5*z**5) - gat_sq*
     .     (1._dp + 5._dp*beta**6 + beta*(4._dp - 10._dp*beta**2 
     .     + 5._dp*beta**4 + 6._dp*beta**6)*z + beta**2*(2._dp 
     .     - 17._dp*beta**2 + 9._dp*beta**4)*z**2 + beta**3*(1._dp 
     .     - 2._dp*beta**2 - 3._dp*beta**4)*z**3 + 6._dp*beta**4*(1._dp 
     .     - beta**2)*z**4 - beta**5*z**5))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta*z)*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + 2._dp*((8._dp*(gat_sq 
     .     + gvt_sq)*mz**4)/s + 4._dp*mz**2*(-gat_sq*(1._dp 
     .     - 4._dp*beta**2 - beta*z) + gvt_sq*(3._dp + beta*z)) 
     .     + s*(gat_sq*(4._dp - 6._dp*beta**2 + 7._dp*beta**4 
     .     - 2._dp*beta*(2._dp - 3._dp*beta**2)*z + beta**2*z**2) 
     .     + gvt_sq*(4._dp + 2._dp*beta**2 - beta**4 
     .     + 2._dp*beta*(2._dp - beta**2)*z + beta**2*z**2)))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep) - 2._dp*((8._dp* 
     .     (gat_sq + gvt_sq)*mz**4*(1._dp + beta*z))/s + 4._dp*mz**2*(1._dp 
     .     + beta*z)*(-gat_sq + 4._dp*beta**2*gat_sq + 3._dp*gvt_sq 
     .     + beta*(gat_sq + gvt_sq)*z) + s*(gat_sq*(4._dp - 6._dp*beta**2 
     .     + 7._dp*beta**4 - beta*(1._dp - 2._dp*beta**2 
     .     - 6._dp*beta**4)*z - 3._dp*beta**2*(1._dp - 2._dp*beta**2)*z**2 
     .     + beta**3*z**3) + gvt_sq*(4._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(7._dp + 2._dp*beta**2 - 2._dp*beta**4)*z + beta**2*(5._dp 
     .     - 2._dp*beta**2)*z**2 + beta**3*z**3)))*xI3(0._dp,mt**2,t,mt**2
     .     ,mt**2,mz**2,musq,ep) + 2._dp*((-2._dp*(gat_sq + gvt_sq)*
     .     mz**4*z*(5._dp - 8._dp*beta**2 - 3._dp*z**2 
     .     + 2._dp*beta**2*z**2))/(beta*s) + beta*s*(gvt_sq*(beta*(4._dp 
     .     - beta**2) + 4._dp*(1._dp + beta**2 - beta**4)*z 
     .     - beta**3*z**2 + beta**2*z**3 + beta**4*z**3) 
     .     - gat_sq*(-3._dp*beta**3 - 4._dp*(1._dp - 3._dp*beta**2 
     .     + 3._dp*beta**4)*z + beta*(4._dp - 3._dp*beta**2)*z**2 
     .     - beta**2*(5._dp - 3._dp*beta**2)*z**3)) + 2._dp*mz**2*
     .     (gvt_sq*(1._dp + beta**2 + 4._dp*beta*z - (1._dp 
     .     - beta**2)*z**2 + 2._dp*beta*z**3) + gat_sq*(1._dp 
     .     + beta**2 - 4._dp*beta*(3._dp - 4._dp*beta**2)*z - (1._dp 
     .     - beta**2)*z**2 + 2._dp*beta*(3._dp - 2._dp*beta**2)*z**3))
     .     )*xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep) + ((16._dp*
     .     (gat_sq + gvt_sq)*mz**6)/s + 8._dp*mz**4*(-gat_sq*(1._dp 
     .     - 5._dp*beta**2 - 2._dp*beta*z) + gvt_sq*(3._dp + beta**2 
     .     + 2._dp*beta*z)) + 2._dp*mz**2*s*(gvt_sq*(4._dp + 9._dp*beta**2 
     .     - 2._dp*beta**4 + 10._dp*beta*z + beta**2*(2._dp + beta**2)*z**2
     .     ) + gat_sq*(4._dp - 7._dp*beta**2 + 14._dp*beta**4 - 2._dp*beta*
     .     (3._dp - 8._dp*beta**2)*z + beta**2*(2._dp + beta**2)*z**2)) 
     .     + s**2*(-gat_sq*(3._dp - 7._dp*beta**2 + 5._dp*beta**4 
     .     - 6._dp*beta**6 - beta*(3._dp - 8._dp*beta**2 
     .     + 12._dp*beta**4)*z + 3._dp*beta**2*(1._dp - beta**2 
     .     - beta**4)*z**2 - beta**3*z**3) + gvt_sq*(1._dp 
     .     + 7._dp*beta**2 - beta**4 - 2._dp*beta**6 + beta*(3._dp 
     .     + 8._dp*beta**2 - 4._dp*beta**4)*z + beta**2*(1._dp 
     .     + 3._dp*beta**2 - beta**4)*z**2 + beta**3*z**3)))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,mt**2,mz**2,
     .     musq,ep)))/pi



      bx(2) = 
     .     (0.5_dp*alpha*genfacbox*gw_sq*sigma0*
     .     (-beta*z*(1._dp - z**2) + (4._dp*beta*(1._dp 
     .     - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
     .     (xI1(mb**2,musq,ep) - xI1(mw**2,musq,ep)))/(s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + (0.5_dp*(4._dp*(-mb**2 
     .     + mw**2)*(6._dp*beta**3 + (5._dp - 8._dp*beta**2 + 4._dp*beta**4
     .     )*z + beta*(5._dp - 12._dp*beta**2)*z**2 - (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 - beta*(3._dp 
     .     - 4._dp*beta**2)*z**4) + s*(2._dp*beta**3*(5._dp - beta**2) 
     .     - beta*(3._dp + 5._dp*beta**2)*z**4 + (1._dp 
     .     - beta**2)*((5._dp + 12._dp*beta**2 - 4._dp*beta**4)*z 
     .     + 9._dp*beta*z**2 - (3._dp + 4._dp*beta**2 
     .     - 2._dp*beta**4)*z**3)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     (beta*s*(1._dp + beta*z)) - (0.5_dp*(4._dp*(-mb**2 
     .     + mw**2)*z*(5._dp - 4._dp*beta**2 - (3._dp - 2._dp*beta**2
     .     )*z**2) - s*(4._dp*beta - (5._dp + 5._dp*beta**2 
     .     - 4._dp*beta**4)*z - 4._dp*beta*z**2 + (3._dp + 5._dp*beta**2 
     .     - 2._dp*beta**4)*z**3))*xI2(s,mb**2,mb**2,musq,ep))/(beta*s) + 
     .     ((4._dp*beta*(-mb**2 + mw**2)*(beta + z)*(1._dp 
     .     - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2)*z 
     .     + 2._dp*beta**2*z**2 + beta**3*z**3) - s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(2._dp + 5._dp*beta**2 - beta**4 + beta*
     .     (3._dp - 6._dp*beta**2 + 2._dp*beta**4)*z - beta**2*(7._dp 
     .     - 2._dp*beta**2)*z**2 + beta**3*(2._dp - beta**2)*z**3 
     .     - beta**4*z**4))*xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp 
     .     + beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + (0.5_dp*(16._dp*
     .     (-mb**2 + mw**2)**2 - 8._dp*beta*mb**2*s*(2._dp*beta + z) 
     .     + 8._dp*mw**2*s*(2._dp + beta**2 + beta*z) + s**2*(7._dp 
     .     + 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp + beta**2)*z 
     .     + 2._dp*beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep))
     .     /s - (0.5_dp*(16._dp*(-mb**2 + mw**2)**2*(1._dp + beta*z) 
     .     - 8._dp*beta*mb**2*s*(beta + z)*(2._dp + beta*z) + 8._dp*mw**2*
     .     s*(1._dp + beta*z)*(2._dp + beta**2 + beta*z) + s**2*(1._dp 
     .     + beta*z)*(7._dp + 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp 
     .     + beta**2)*z + 2._dp*beta**2*z**2))*xI3(0._dp,mt**2,t,mb**2
     .     ,mb**2,mw**2,musq,ep))/s + (0.125_dp*(-16._dp*(-mb**2 
     .     + mw**2)**2*z*(5._dp - 8._dp*beta**2 - (3._dp 
     .     - 2._dp*beta**2)*z**2) - 8._dp*mb**2*s*(2._dp*beta*(1._dp 
     .     + beta**2) - (5._dp - 8._dp*beta**2)*(1._dp + beta**2)*z 
     .     - 2._dp*beta*(1._dp - beta**2)*z**2 + (3._dp 
     .     - 2._dp*beta**2)*(1._dp + beta**2)*z**3) + 8._dp*mw**2*s*
     .     (2._dp*beta*(1._dp + beta**2) - (5._dp - 5._dp*beta**2 
     .     - 8._dp*beta**4)*z - 2._dp*beta*(1._dp - beta**2)*z**2 + 
     .     (3._dp + 3._dp*beta**2 - 2._dp*beta**4)*z**3) - s**2*
     .     (-4._dp*beta*(1._dp + 4._dp*beta**2 + beta**4) + (5._dp 
     .     - 30._dp*beta**2 + beta**4 - 8._dp*beta**6)*z + 4._dp*beta*
     .     (1._dp + 2._dp*beta**2 - beta**4)*z**2 - (3._dp 
     .     + 4._dp*beta**2 + 11._dp*beta**4 - 2._dp*beta**6)*z**3))*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep))/(beta*s) 
     .     + (0.125_dp*(64._dp*(-mb**2 + mw**2)**3 + s**3*(1._dp 
     .     + beta**2 + 2._dp*beta*z)*(7._dp + 2._dp*beta**2 + beta**4 
     .     + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2) 
     .     + 8._dp*s*(10._dp*mw**4 + 2._dp*(-mb**2 + mw**2)**2*
     .     (3._dp*beta**2 + 4._dp*beta*z) - 4._dp*mb**2*mw**2*(3._dp 
     .     + beta**2*z**2) + 2._dp*mb**4*(1._dp + 2._dp*beta**2*z**2)) 
     .     + 4._dp*s**2*(mw**2*(11._dp + 8._dp*beta**2 + 3._dp*beta**4 
     .     + 4._dp*beta*(3._dp + 2._dp*beta**2)*z + 6._dp*beta**2*z**2) 
     .     - mb**2*(11._dp - 4._dp*beta**2 + 3._dp*beta**4 
     .     + 8._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*
     .     (6._dp + beta**2)*z**2)))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mb**2,mb**2,mb**2,mw**2,
     .     musq,ep))/s))/pi


!       bx(3) = 
!      .     (0.5_dp*alpha*chifac*genfacbox*mt**2*sigma0*
!      .     (-beta*z*(1._dp - z**2) + (4._dp*beta*(1._dp 
!      .     - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
!      .     (xI1(mt**2,musq,ep) - xI1(mz**2,musq,ep)))/(s*(1._dp 
!      .     + beta**2 + 2._dp*beta*z)) + ((-2._dp*beta*(1._dp 
!      .     - beta**2)*s*(2._dp - beta**2 + beta*z - z**2 
!      .     - beta*z**3) + 2._dp*mz**2*(6._dp*beta**3 + (5._dp 
!      .     - 8._dp*beta**2 + 4._dp*beta**4)*z + beta*(5._dp 
!      .     - 12._dp*beta**2)*z**2 - (3._dp - 4._dp*beta**2 
!      .     + 2._dp*beta**4)*z**3 - beta*(3._dp - 4._dp*beta**2)*z**4)
!      .     )*xI2(mt**2,mt**2,mz**2,musq,ep))/(beta*s*(1._dp + beta*z)) 
!      .     - (1._dp*(-beta*s*(2._dp + beta*z)*(1._dp - z**2) 
!      .     + 2._dp*mz**2*z*(5._dp - 4._dp*beta**2 - (3._dp 
!      .     - 2._dp*beta**2)*z**2))*xI2(s,mt**2,mt**2,musq,ep))/(beta*s)
!      .     + ((4._dp*beta*mz**2*(beta + z)*(1._dp - 3._dp*beta**2 
!      .     + beta*(1._dp - 2._dp*beta**2)*z + 2._dp*beta**2*z**2 
!      .     + beta**3*z**3) + 2._dp*s*(1._dp - 4._dp*beta**2 + beta**6 
!      .     + beta*(2._dp - 10._dp*beta**2 + 3._dp*beta**4)*z + beta**2*
!      .     (3._dp - 5._dp*beta**2)*z**2 + 2._dp*beta**3*(3._dp 
!      .     - beta**2)*z**3 + 4._dp*beta**4*z**4 + beta**5*z**5))*
!      .     xI2(t,mt**2,mz**2,musq,ep))/(s*(1._dp + beta*z)*(1._dp 
!      .     + beta**2 + 2._dp*beta*z)) + ((8._dp*mz**4 + 4._dp*beta*mz**2*s*
!      .     (beta + z) + s**2*(2._dp - 2._dp*beta**2 + beta**4 
!      .     + 2._dp*beta*z + beta**2*z**2))*
!      .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/s
!      .     - (1._dp*(s**2*(2._dp - 2._dp*beta**2 + beta**4 
!      .     + 3._dp*beta*z + 3._dp*beta**2*z**2 + beta**3*z**3) + 4._dp*
!      .     (1._dp + beta*z)*(2._dp*mz**4 + beta*mz**2*s*(beta + z)))*
!      .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep))/s - 
!      .     (1._dp*(-beta**2*s**2*(beta + z)*(1._dp + beta*z) 
!      .     + 2._dp*mz**4*z*(5._dp - 8._dp*beta**2 - (3._dp 
!      .     - 2._dp*beta**2)*z**2) - 2._dp*beta*mz**2*s*(1._dp + beta**2 
!      .     - 2._dp*beta*(1._dp - 2._dp*beta**2)*z - (1._dp 
!      .     - beta**2)*z**2 + beta*(1._dp - beta**2)*z**3))*
!      .     xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep))/(beta*s) + 
!      .     (0.5_dp*(16._dp*mz**6 + 16._dp*beta*mz**4*s*(beta + z) 
!      .     + beta*s**3*(beta + z)*(1._dp + beta*z)**2 + 2._dp*mz**2*s**2*
!      .     (2._dp - beta**2 + 2._dp*beta**4 + 2._dp*beta*(1._dp 
!      .     + 2._dp*beta**2)*z + beta**2*(2._dp + beta**2)*z**2))*
!      .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,mt**2,mz**2,musq,
!      .     ep))/s))/pi

         MARKUS_bx3= alpha/(32d0*MZ**2*Pi*SW2)/8d0*
     -    (    (-2*MZ**2*s*(-7 + 9*beta*z)*
     -     (3*(4*C33phiq3**2 + C33phiu**2)*vol4 + 
     -       beta**5*vol2*z*(-4*C33phiq3 - 2*C33phiu + 
     -          2*(-2*C33phiq3 + C33phiu)**2*vol2 + 
     -          (-4*C33phiq3 + 12*C33phiq3**2*vol2 + 
     -             C33phiu*(-2 + 3*C33phiu*vol2))*z**2) + 
     -       beta**3*z*(-1 + z**2 + 4*C33phiu*vol2*(2 + z**2) - 
     -          8*C33phiq3**2*vol4*(4 + 5*z**2) - 
     -          2*C33phiu**2*vol4*(4 + 5*z**2) + 
     -          4*C33phiq3*vol2*
     -           (4 + 3*C33phiu*vol2 + (2 + C33phiu*vol2)*z**2)) - 
     -       beta*z*(-1 + z**2 + 4*C33phiq3**2*vol4*(-14 + z**2) + 
     -          C33phiu**2*vol4*(-14 + z**2) + 
     -          2*C33phiu*vol2*(3 + z**2) + 
     -          4*C33phiq3*vol2*(3 + z**2 + C33phiu*vol2*(1 + z**2)))
     -        - beta**2*(-z**2 + z**4 + 2*C33phiu*vol2*(1 + z**4) + 
     -          4*C33phiq3**2*vol4*(4 - 7*z**2 + z**4) + 
     -          C33phiu**2*vol4*(4 - 7*z**2 + z**4) + 
     -          4*C33phiq3*vol2*
     -           (1 + z**4 + C33phiu*vol2*(2 - 3*z**2 + z**4))) + 
     -       beta**4*(-z**2 + z**4 + 2*C33phiu*vol2*(1 + z**4) + 
     -          4*C33phiq3**2*vol4*(2 - 9*z**2 + 2*z**4) + 
     -          C33phiu**2*vol4*(2 - 9*z**2 + 2*z**4) + 
     -          4*C33phiq3*vol2*
     -           (1 + z**4 + C33phiu*vol2*(2 - 3*z**2 + z**4)))))/
     -   (3.*(-1 + beta*z)*(MW + beta*MW*z)**2) - 
     -  (8*MZ**2*(-7 + 9*beta*z)*
     -     (4*beta**2 - 8*beta**4 + 4*beta**6 + 
     -       12*beta**2*C33phiq3*vol2 - 12*beta**4*C33phiq3*vol2 + 
     -       6*beta**2*C33phiu*vol2 - 6*beta**4*C33phiu*vol2 + 
     -       8*beta**2*C33phiq3*C33phiu*vol2**2 + 
     -       8*beta**4*C33phiq3*C33phiu*vol2**2 + 
     -       4*C33phiq3**2*vol4 - 72*beta**4*C33phiq3**2*vol4 + 
     -       36*beta**6*C33phiq3**2*vol4 - 
     -       16*beta**6*C33phiq3*C33phiu*vol4 + C33phiu**2*vol4 - 
     -       18*beta**4*C33phiu**2*vol4 + 9*beta**6*C33phiu**2*vol4 + 
     -       beta*((-1 + beta**2)**2*(-1 + 4*beta**2) - 
     -          2*(1 - 3*beta**2 + 2*beta**6)*(2*C33phiq3 + C33phiu)*
     -           vol2 - 4*(1 - 9*beta**4 + 8*beta**6)*C33phiq3*
     -           C33phiu*vol2**2 + 
     -          beta**2*(-6 - 11*beta**2 + 9*beta**4)*
     -           (4*C33phiq3**2 + C33phiu**2)*vol4)*z - 
     -       beta**2*(7 + 4*C33phiq3*vol2*(6 + 5*C33phiu*vol2) + 
     -          24*C33phiq3**2*vol4 + 
     -          6*C33phiu*(2*vol2 + C33phiu*vol4) - 
     -          2*beta**2*(7 + 6*C33phiu*vol2 - 
     -             4*C33phiq3*vol2*(-3 + C33phiu*vol2) + 
     -             44*C33phiq3**2*vol4 + 11*C33phiu**2*vol4) + 
     -          beta**4*(7 + 
     -             4*(16*C33phiq3**2 - 7*C33phiq3*C33phiu + 
     -                4*C33phiu**2)*vol4))*z**2 - 
     -       beta*(-1 + beta**2)*
     -        (1 + 2*C33phiu*vol2 + 
     -          4*C33phiq3*vol2*(1 + C33phiu*vol2) + 
     -          4*C33phiq3**2*vol4 + C33phiu**2*vol4 + 
     -          beta**2*(-7 - 8*C33phiu*vol2 - 
     -             4*C33phiq3*vol2*(4 + C33phiu*vol2) + 
     -             8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4) + 
     -          3*beta**4*(2 - 
     -             2*vol2*(C33phiu + C33phiq3*(2 + 8*C33phiu*vol2)) + 
     -             5*(4*C33phiq3**2 + C33phiu**2)*vol4))*z**3 + 
     -       beta**2*(3*(1 + 2*C33phiu*vol2 + 
     -             4*C33phiq3*vol2*(1 + C33phiu*vol2) + 
     -             4*C33phiq3**2*vol4 + C33phiu**2*vol4) - 
     -          beta**2*(6 + 6*(2*C33phiq3 + C33phiu)*vol2 + 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4) + 
     -          3*beta**4*(1 + 
     -             2*(4*C33phiq3**2 - 2*C33phiq3*C33phiu + 
     -                C33phiu**2)*vol4))*z**4 + 
     -       beta**3*(2*(-1 + beta**2)**2 - 
     -          2*(-2 + beta**2 + beta**4)*(2*C33phiq3 + C33phiu)*
     -           vol2 + 8*(1 + beta**2 - 2*beta**4)*C33phiq3*C33phiu*
     -           vol2**2 + (2 + beta**2 + 5*beta**4)*
     -           (4*C33phiq3**2 + C33phiu**2)*vol4)*z**5)*
     -     xI1(MT**2,musq,ep))/
     -   (3.*(-1 + beta**2)*(-1 + beta*z)*(1 + beta**2 + 2*beta*z)*
     -     (MW + beta*MW*z)**2) + 
     -  (8*MZ**2*(-7 + 9*beta*z)*
     -     (4*beta**2 - 8*beta**4 + 4*beta**6 + 
     -       12*beta**2*C33phiq3*vol2 - 12*beta**4*C33phiq3*vol2 + 
     -       6*beta**2*C33phiu*vol2 - 6*beta**4*C33phiu*vol2 + 
     -       8*beta**2*C33phiq3*C33phiu*vol2**2 + 
     -       8*beta**4*C33phiq3*C33phiu*vol2**2 + 
     -       4*C33phiq3**2*vol4 - 72*beta**4*C33phiq3**2*vol4 + 
     -       36*beta**6*C33phiq3**2*vol4 - 
     -       16*beta**6*C33phiq3*C33phiu*vol4 + C33phiu**2*vol4 - 
     -       18*beta**4*C33phiu**2*vol4 + 9*beta**6*C33phiu**2*vol4 + 
     -       beta*((-1 + beta**2)**2*(-1 + 4*beta**2) - 
     -          2*(1 - 3*beta**2 + 2*beta**6)*(2*C33phiq3 + C33phiu)*
     -           vol2 - 4*(1 - 9*beta**4 + 8*beta**6)*C33phiq3*
     -           C33phiu*vol2**2 + 
     -          beta**2*(-6 - 11*beta**2 + 9*beta**4)*
     -           (4*C33phiq3**2 + C33phiu**2)*vol4)*z - 
     -       beta**2*(7 + 4*C33phiq3*vol2*(6 + 5*C33phiu*vol2) + 
     -          24*C33phiq3**2*vol4 + 
     -          6*C33phiu*(2*vol2 + C33phiu*vol4) - 
     -          2*beta**2*(7 + 6*C33phiu*vol2 - 
     -             4*C33phiq3*vol2*(-3 + C33phiu*vol2) + 
     -             44*C33phiq3**2*vol4 + 11*C33phiu**2*vol4) + 
     -          beta**4*(7 + 
     -             4*(16*C33phiq3**2 - 7*C33phiq3*C33phiu + 
     -                4*C33phiu**2)*vol4))*z**2 - 
     -       beta*(-1 + beta**2)*
     -        (1 + 2*C33phiu*vol2 + 
     -          4*C33phiq3*vol2*(1 + C33phiu*vol2) + 
     -          4*C33phiq3**2*vol4 + C33phiu**2*vol4 + 
     -          beta**2*(-7 - 8*C33phiu*vol2 - 
     -             4*C33phiq3*vol2*(4 + C33phiu*vol2) + 
     -             8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4) + 
     -          3*beta**4*(2 - 
     -             2*vol2*(C33phiu + C33phiq3*(2 + 8*C33phiu*vol2)) + 
     -             5*(4*C33phiq3**2 + C33phiu**2)*vol4))*z**3 + 
     -       beta**2*(3*(1 + 2*C33phiu*vol2 + 
     -             4*C33phiq3*vol2*(1 + C33phiu*vol2) + 
     -             4*C33phiq3**2*vol4 + C33phiu**2*vol4) - 
     -          beta**2*(6 + 6*(2*C33phiq3 + C33phiu)*vol2 + 
     -             (4*C33phiq3**2 + C33phiu**2)*vol4) + 
     -          3*beta**4*(1 + 
     -             2*(4*C33phiq3**2 - 2*C33phiq3*C33phiu + 
     -                C33phiu**2)*vol4))*z**4 + 
     -       beta**3*(2*(-1 + beta**2)**2 - 
     -          2*(-2 + beta**2 + beta**4)*(2*C33phiq3 + C33phiu)*
     -           vol2 + 8*(1 + beta**2 - 2*beta**4)*C33phiq3*C33phiu*
     -           vol2**2 + (2 + beta**2 + 5*beta**4)*
     -           (4*C33phiq3**2 + C33phiu**2)*vol4)*z**5)*
     -     xI1(MZ**2,musq,ep))/
     -   (3.*(-1 + beta**2)*(-1 + beta*z)*(1 + beta**2 + 2*beta*z)*
     -     (MW + beta*MW*z)**2) + 
     -  (4*MZ**2*(-7 + 9*beta*z)*
     -     (2*beta**11*(-2*C33phiq3 + C33phiu)**2*s*vol4*z**2*
     -        (-2 + z**2) + (MZ + (2*C33phiq3 + C33phiu)*MZ*vol2)**2*
     -        z*(-5 + 3*z**2) + 
     -       beta**10*s*z*(1 - 2*C33phiq3*vol2 - C33phiu*vol2 - 
     -          40*C33phiq3**2*vol4 + 40*C33phiq3*C33phiu*vol4 - 
     -          10*C33phiu**2*vol4 + 
     -          14*(-2*C33phiq3 + C33phiu)**2*vol4*z**2 - 
     -          4*(-2*C33phiq3 + C33phiu)**2*vol4*z**4) + 
     -       beta**7*(-5*s + 
     -          2*(-2*s*vol2*
     -              (2*C33phiq3 + C33phiu + 15*C33phiq3*C33phiu*vol2)
     -              + 10*(4*C33phiq3**2 + C33phiu**2)*s*vol4 + 
     -             MZ**2*(-3 + 
     -                3*vol2*
     -                 (2*C33phiq3 + C33phiu + 
     -                   8*C33phiq3*C33phiu*vol2) - 
     -                4*(4*C33phiq3**2 + C33phiu**2)*vol4)) + 
     -          (28*MZ**2*(1 - 4*C33phiq3*C33phiu*vol2**2) + 
     -             2*s*(2 - 5*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(-5 + 22*C33phiu*vol2)) + 
     -             3*(4*C33phiq3**2 + C33phiu**2)*(10*MZ**2 - 17*s)*
     -              vol4)*z**2 - 
     -          (12*MZ**2 + 3*s + 
     -             2*(2*C33phiq3 + C33phiu)*(MZ**2 - 5*s)*vol2 + 
     -             4*C33phiq3*C33phiu*(-10*MZ**2 + s)*vol2**2 + 
     -             (4*C33phiq3**2 + C33phiu**2)*(12*MZ**2 - 13*s)*vol4
     -             )*z**4) + 
     -       beta**3*(2*MZ**2*
     -           (-3 + 12*C33phiq3*C33phiu*vol2**2 + 
     -             (20 + 26*C33phiu*vol2 + 
     -                4*C33phiq3*vol2*(13 + 6*C33phiu*vol2) + 
     -                60*C33phiq3**2*vol4 + 15*C33phiu**2*vol4)*z**2
     -              - 2*(5 + 8*C33phiu*vol2 + 
     -                4*C33phiq3*vol2*(4 + 3*C33phiu*vol2) + 
     -                16*C33phiq3**2*vol4 + 4*C33phiu**2*vol4)*z**4)
     -           + s*(-7 - 4*vol2*
     -              (2*C33phiu + C33phiq3*(4 + C33phiu*vol2)) - 
     -             2*(4*C33phiq3**2 + C33phiu**2)*vol4 - 
     -             (-4 + 2*vol2*
     -                 (C33phiu + 2*C33phiq3*(1 + 6*C33phiu*vol2)) + 
     -                9*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -             (-1 + 6*C33phiu*vol2 + 
     -                4*C33phiq3*vol2*(3 + 7*C33phiu*vol2) + 
     -                28*C33phiq3**2*vol4 + 7*C33phiu**2*vol4)*z**4))
     -        + beta*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)*
     -        (2*MZ**2*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)*z**2*
     -           (-5 + 3*z**2) + 
     -          s*(2 - z**2 - (2*C33phiq3 + C33phiu)*vol2*(-1 + z**2))
     -          ) + beta**2*z*
     -        (MZ**2*(18 + 6*vol2*
     -              (5*C33phiu + 2*C33phiq3*(5 + 4*C33phiu*vol2)) + 
     -             14*(4*C33phiq3**2 + C33phiu**2)*vol4 - 
     -             (15 + 4*C33phiq3*vol2*(13 + 11*C33phiu*vol2) + 
     -                52*C33phiq3**2*vol4 + 
     -                13*C33phiu*(2*vol2 + C33phiu*vol4))*z**2 + 
     -             3*(1 + 2*C33phiu*vol2 + 
     -                4*C33phiq3*vol2*(1 + C33phiu*vol2) + 
     -                4*C33phiq3**2*vol4 + C33phiu**2*vol4)*z**4) + 
     -          s*(3 - 2*z**2 + 8*C33phiq3**2*vol4*(-2 + z**2) + 
     -             2*C33phiq3*vol2*
     -              (-1 + 4*C33phiu*vol2*(-2 + z**2)) + 
     -             C33phiu*(-vol2 + 2*C33phiu*vol4*(-2 + z**2)))) + 
     -       beta**5*(s*(9 + 8*C33phiu*vol2 - 13*C33phiu**2*vol4 + 
     -             3*(-2 + 4*C33phiu*vol2 + 13*C33phiu**2*vol4)*
     -              z**2 + (3 - 2*C33phiu*(7*vol2 + 6*C33phiu*vol4))*
     -              z**4 - 4*C33phiq3**2*vol4*
     -              (13 - 39*z**2 + 12*z**4) + 
     -             4*C33phiq3*vol2*
     -              (4 + 7*C33phiu*vol2 + 2*(3 + C33phiu*vol2)*z**2 - 
     -                (7 + 11*C33phiu*vol2)*z**4)) + 
     -          2*MZ**2*(6 - 3*C33phiu*vol2 + 8*C33phiu**2*vol4 - 
     -             (27 + 16*C33phiu*vol2 + 23*C33phiu**2*vol4)*z**2 + 
     -             (12 + 11*C33phiu*vol2 + 6*C33phiu**2*vol4)*z**4 + 
     -             4*C33phiq3**2*vol4*(8 - 23*z**2 + 6*z**4) - 
     -             2*C33phiq3*vol2*
     -              (3 + 16*z**2 - 11*z**4 + 
     -                2*C33phiu*vol2*(9 - 11*z**2 + z**4)))) + 
     -       beta**4*z*(MZ**2*
     -           (-31 - 34*C33phiu*vol2 - 25*C33phiu**2*vol4 + 
     -             (35 + 40*C33phiu*vol2 + 33*C33phiu**2*vol4)*z**2 - 
     -             2*(5 + 8*C33phiu*vol2 + 4*C33phiu**2*vol4)*z**4 - 
     -             4*C33phiq3**2*vol4*(25 - 33*z**2 + 8*z**4) - 
     -             4*C33phiq3*vol2*
     -              (17 + 3*C33phiu*vol2 - 
     -                5*(4 + C33phiu*vol2)*z**2 + 
     -                2*(4 + 3*C33phiu*vol2)*z**4)) + 
     -          2*s*(-5 + 3*z**2 + 
     -             4*C33phiq3**2*vol4*(9 - 8*z**2 + 2*z**4) + 
     -             2*C33phiq3*vol2*
     -              (3 - 3*z**2 + 2*z**4 + 
     -                2*C33phiu*vol2*(6 - 5*z**2 + 2*z**4)) + 
     -             C33phiu*(C33phiu*vol4*(9 - 8*z**2 + 2*z**4) + 
     -                vol2*(3 - 3*z**2 + 2*z**4)))) + 
     -       2*beta**8*z*(s*(-3 + z**2 + 
     -             4*C33phiq3**2*vol4*(19 - 20*z**2 + 2*z**4) + 
     -             2*C33phiq3*vol2*
     -              (3 - 16*C33phiu*vol2 + 
     -                (-3 + 22*C33phiu*vol2)*z**2 + 
     -                2*(1 - 2*C33phiu*vol2)*z**4) + 
     -             C33phiu*(3*vol2 + 19*C33phiu*vol4 - 
     -                (3*vol2 + 20*C33phiu*vol4)*z**2 + 
     -                2*(vol2 + C33phiu*vol4)*z**4)) + 
     -          MZ**2*(-5 + C33phiu*vol2 - 6*C33phiu**2*vol4 + 
     -             (7 - 3*C33phiu*vol2 + 10*C33phiu**2*vol4)*z**2 + 
     -             (-2 + C33phiu*(vol2 - 3*C33phiu*vol4))*z**4 - 
     -             4*C33phiq3**2*vol4*(6 - 10*z**2 + 3*z**4) + 
     -             2*C33phiq3*vol2*
     -              (1 - 3*z**2 + z**4 + 
     -                2*C33phiu*vol2*(6 - 10*z**2 + 3*z**4)))) + 
     -       beta**6*z*(MZ**2*
     -           (4*(7 + 2*C33phiq3*vol2*(3 - 8*C33phiu*vol2) + 
     -                36*C33phiq3**2*vol4 + 
     -                3*C33phiu*(vol2 + 3*C33phiu*vol4)) - 
     -             (37 + 2*vol2*
     -                 (7*C33phiu + C33phiq3*(14 - 46*C33phiu*vol2))
     -                 + 43*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -             (11 + 8*C33phiu*vol2 + 
     -                4*C33phiq3*vol2*(4 - 3*C33phiu*vol2) + 
     -                12*C33phiq3**2*vol4 + 3*C33phiu**2*vol4)*z**4)
     -           + 2*s*(6 - 3*z**2 - 
     -             4*C33phiq3**2*vol4*(21 - 20*z**2 + 2*z**4) - 
     -             2*C33phiq3*vol2*
     -              (5 + 2*C33phiu*vol2 - 6*z**2 + 
     -                4*(1 + C33phiu*vol2)*z**4) - 
     -             C33phiu*(C33phiu*vol4*(21 - 20*z**2 + 2*z**4) + 
     -                vol2*(5 - 6*z**2 + 4*z**4)))) + 
     -       beta**9*(2*MZ**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4)*
     -           z**2*(-2 + z**2) + 
     -          s*(1 + C33phiu*vol2 - 6*C33phiu**2*vol4 + 
     -             (-1 + 2*C33phiu*(vol2 + 13*C33phiu*vol4))*z**2 + 
     -             (1 - 2*C33phiu*(vol2 + 5*C33phiu*vol4))*z**4 - 
     -             8*C33phiq3**2*vol4*(3 - 13*z**2 + 5*z**4) + 
     -             2*C33phiq3*vol2*
     -              (1 + 2*z**2 - 2*z**4 + 
     -                2*C33phiu*vol2*(8 + 7*z**2*(-3 + z**2))))))*
     -     xI2(MT**2,MT**2,MZ**2,musq,ep))/
     -   (3.*beta*(-1 + beta**2)*MW**2*(-1 + beta*z)*(1 + beta*z)**3)
     -   + (2*MZ**2*(-7 + 9*beta*z)*
     -     (4*beta**6*(-2*C33phiq3 + C33phiu)**2*s*vol4*z*
     -        (-2 + z**2) - 2*beta*s*
     -        (1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*(-1 + z**2) - 
     -       2*beta**3*s*(-1 - 2*(2*C33phiq3 + C33phiu)*vol2 + 
     -          (-2*C33phiq3 + C33phiu)**2*vol4)*(-1 + z**2) + 
     -       2*(MZ + (2*C33phiq3 + C33phiu)*MZ*vol2)**2*z*
     -        (-5 + 3*z**2) + 
     -       beta**4*z*(s*(-1 + 10*C33phiu*vol2 + 
     -             4*C33phiq3*vol2*(5 + 3*C33phiu*vol2)) + 
     -          13*(4*C33phiq3**2 + C33phiu**2)*s*vol4 - 
     -          s*(-1 + 6*vol2*
     -              (C33phiu + 2*C33phiq3*(1 + C33phiu*vol2)) + 
     -             5*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -          4*MZ**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4)*
     -           (-2 + z**2)) + 
     -       beta**2*z*(2*MZ**2*
     -           (9 + 10*C33phiu*vol2 + 7*C33phiu**2*vol4 - 
     -             (5 + 6*C33phiu*vol2 + 3*C33phiu**2*vol4)*z**2 + 
     -             4*C33phiq3**2*vol4*(7 - 3*z**2) + 
     -             4*C33phiq3*vol2*
     -              (5 + C33phiu*vol2 - (3 + C33phiu*vol2)*z**2)) + 
     -          s*(4*C33phiq3**2*vol4*(-11 + 7*z**2) + 
     -             (1 + C33phiu*vol2)*
     -              (1 - 11*C33phiu*vol2 + (-1 + 7*C33phiu*vol2)*z**2)
     -               + 4*C33phiq3*vol2*
     -              (-5 - 11*C33phiu*vol2 + (3 + 7*C33phiu*vol2)*z**2)
     -             )))*xI2(s,MT**2,MT**2,musq,ep))/
     -   (3.*beta*MW**2*(-1 + beta**2*z**2)) - 
     -  (4*MZ**2*(-7 + 9*beta*z)*
     -     (-(s*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)) - 
     -       2*(4*C33phiq3**2 + C33phiu**2)*MZ**2*vol4 - 
     -       beta*(2*MZ**2*(1 + 6*C33phiq3*vol2 + 3*C33phiu*vol2 + 
     -             2*(2*C33phiq3 + C33phiu)**2*vol4) + 
     -          s*(3 + 4*C33phiq3*vol2 + 8*C33phiq3**2*vol4 + 
     -             2*C33phiu*(vol2 + C33phiu*vol4)))*z + 
     -       2*beta**10*(-2*C33phiq3 + C33phiu)**2*s*vol4*z**2*
     -        (-2 + z**2) + beta**5*z*
     -        (2*MZ**2*(2*(2 + 5*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(5 + 6*C33phiu*vol2) + 
     -                28*C33phiq3**2*vol4 + 7*C33phiu**2*vol4) - 
     -             2*(-1 + vol2*
     -                 (C33phiu + C33phiq3*(2 + 8*C33phiu*vol2)) + 
     -                6*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 - 
     -             (1 + 2*C33phiq3*vol2 + C33phiu*vol2)*z**4) + 
     -          s*(-17 - 14*C33phiq3*vol2 - 7*C33phiu*vol2 - 
     -             76*C33phiq3*C33phiu*vol2**2 + 
     -             128*C33phiq3**2*vol4 + 32*C33phiu**2*vol4 + 
     -             2*(8 + vol2*
     -                 (C33phiu + C33phiq3*(2 + 36*C33phiu*vol2)) - 
     -                12*(4*C33phiq3**2 + C33phiu**2)*vol4)*z**2 + 
     -             (-5 + 4*vol2*
     -                 (C33phiu + C33phiq3*(2 + 3*C33phiu*vol2)))*z**4
     -             )) + beta**2*
     -        (s*(5 + 8*C33phiq3*vol2 + 4*C33phiu*vol2 + 
     -             8*C33phiq3**2*vol4 + 2*C33phiu**2*vol4 - 
     -             (5 + 5*C33phiu*vol2 + 
     -                2*C33phiq3*vol2*(5 + 8*C33phiu*vol2) + 
     -                32*C33phiq3**2*vol4 + 8*C33phiu**2*vol4)*z**2)
     -           + 2*MZ**2*(-1 - 2*z**2 + 
     -             4*C33phiq3**2*vol4*(-1 + z**2) + 
     -             C33phiu**2*vol4*(-1 + z**2) - 
     -             2*C33phiu*vol2*(3 + z**2) - 
     -             4*C33phiq3*vol2*(3 + 5*C33phiu*vol2 + z**2))) + 
     -       beta**9*s*z*(1 - C33phiu*vol2 + 
     -          8*C33phiq3**2*vol4*(-5 + z**2 + z**4) + 
     -          2*C33phiu**2*vol4*(-5 + z**2 + z**4) - 
     -          2*C33phiq3*(vol2 + 4*C33phiu*vol4*(-5 + z**2 + z**4)))
     -         + beta**7*z*(2*MZ**2*
     -           (-5 + z**2 + z**4 + 
     -             8*C33phiq3**2*vol4*(-3 + 2*z**2) + 
     -             2*C33phiu**2*vol4*(-3 + 2*z**2) + 
     -             C33phiu*vol2*(1 - 3*z**2 + z**4) + 
     -             2*C33phiq3*vol2*
     -              (1 + 12*C33phiu*vol2 - 
     -                (3 + 8*C33phiu*vol2)*z**2 + z**4)) + 
     -          s*(2 - 7*z**2 + 5*z**4 + 
     -             C33phiu*vol2*(6 + 7*z**2 - 4*z**4) + 
     -             2*C33phiq3*vol2*
     -              (6 + 36*C33phiu*vol2 + 
     -                (7 - 22*C33phiu*vol2)*z**2 - 
     -                2*(2 + C33phiu*vol2)*z**4) + 
     -             4*C33phiq3**2*vol4*(10 + 2*z**2 - 9*z**4 + z**6) + 
     -             C33phiu**2*vol4*(10 + 2*z**2 - 9*z**4 + z**6))) - 
     -       beta**3*z*(2*MZ**2*
     -           (-2 + 8*vol2*
     -              (C33phiu + C33phiq3*(2 + 5*C33phiu*vol2)) - 
     -             10*(4*C33phiq3**2 + C33phiu**2)*vol4 + 
     -             (3 - 5*C33phiu*vol2 - 
     -                2*C33phiq3*vol2*(5 + 16*C33phiu*vol2) + 
     -                32*C33phiq3**2*vol4 + 8*C33phiu**2*vol4)*z**2)
     -           + s*(-17 - 4*C33phiu*vol2 + 3*C33phiu**2*vol4 + 
     -             (9 + 9*C33phiu*vol2 + C33phiu**2*vol4)*z**2 + 
     -             4*C33phiq3**2*vol4*(3 + z**2) + 
     -             2*C33phiq3*vol2*
     -              (-4 + 9*z**2 + 2*C33phiu*vol2*(9 + 5*z**2)))) + 
     -       beta**6*(-(s*(1 + 18*z**2 - 12*z**4 + z**6 + 
     -               C33phiu*vol2*(1 - 14*z**2 + z**4 - z**6) + 
     -               8*C33phiq3**2*vol4*
     -                (-6 - 16*z**2 + 17*z**4 + z**6) + 
     -               2*C33phiu**2*vol4*
     -                (-6 - 16*z**2 + 17*z**4 + z**6) + 
     -               2*C33phiq3*vol2*
     -                (1 + 2*C33phiu*vol2 - 
     -                  2*(7 + 2*C33phiu*vol2)*z**2 + 
     -                  (1 - 22*C33phiu*vol2)*z**4 + 
     -                  (-1 + 4*C33phiu*vol2)*z**6))) + 
     -          2*MZ**2*(-3 + 2*z**4 - 
     -             C33phiu*vol2*(-3 + 2*z**2 + z**4) + 
     -             4*C33phiq3**2*vol4*(-4 + 4*z**2 - z**4 + z**6) + 
     -             C33phiu**2*vol4*(-4 + 4*z**2 - z**4 + z**6) - 
     -             2*C33phiq3*vol2*
     -              (-3 + 2*z**2 + z**4 + 
     -                2*C33phiu*vol2*(-6 + 2*z**2 + 3*z**4)))) + 
     -       beta**4*(2*MZ**2*
     -           (4 + 3*C33phiu*vol2 + 9*C33phiu**2*vol4 + 
     -             4*(1 + C33phiu*(vol2 + 2*C33phiu*vol4))*z**2 + 
     -             (-3 + C33phiu*(vol2 - 15*C33phiu*vol4))*z**4 + 
     -             4*C33phiq3**2*vol4*(9 + 8*z**2 - 15*z**4) + 
     -             2*C33phiq3*vol2*
     -              (3 - 2*C33phiu*vol2 + 4*z**2 + z**4 + 
     -                8*C33phiu*vol2*z**4)) + 
     -          s*(-4 - 3*C33phiu*vol2 + 7*C33phiu**2*vol4 + 
     -             (20 - 8*C33phiu*vol2 - 15*C33phiu**2*vol4)*z**2 + 
     -             2*(-5 + 6*C33phiu**2*vol4)*z**4 + 
     -             4*C33phiq3**2*vol4*(7 - 15*z**2 + 12*z**4) + 
     -             2*C33phiq3*vol2*
     -              (-3 - 8*z**2 + 
     -                2*C33phiu*vol2*(-7 - 13*z**2 + 6*z**4)))) + 
     -       beta**8*(2*MZ**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4)*
     -           z**2*(-2 + z**2) + 
     -          s*(1 + 3*z**2 - 2*z**4 + z**6 - 
     -             C33phiu*vol2*(-1 + z)*(1 + z)*(1 + z**4) - 
     -             4*C33phiq3**2*vol4*(6 + 6*z**2 - 9*z**4 + z**6) - 
     -             C33phiu**2*vol4*(6 + 6*z**2 - 9*z**4 + z**6) + 
     -             2*C33phiq3*vol2*
     -              (-((-1 + z)*(1 + z)*(1 + z**4)) + 
     -                2*C33phiu*vol2*(8 + 11*z**2 - 15*z**4 + 2*z**6))
     -             )))*xI2(MT**2 - (s*(1d0 + beta*z))/2d0,MT**2,MZ**2,
     -      musq,ep))/
     -   (3.*MW**2*(-1 + beta*z)*(1 + beta*z)**3*
     -     (1 + beta**2 + 2*beta*z)) - 
     -  (2*MZ**2*(-7 + 9*beta*z)*
     -     (8*MZ**4*(-1 - 12*C33phiq3**2*vol4 + 
     -          (4*C33phiq3 - 3*C33phiu)*C33phiu*vol4 + 
     -          beta**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4)) + 
     -       4*beta*MZ**2*s*(beta*
     -           (-1 - 36*C33phiq3**2*vol4 + 
     -             3*(4*C33phiq3 - 3*C33phiu)*C33phiu*vol4 + 
     -             beta**2*(1 + 3*(-2*C33phiq3 + C33phiu)**2*vol4)) + 
     -          (-1 - 12*C33phiq3**2*vol4 + 
     -             (4*C33phiq3 - 3*C33phiu)*C33phiu*vol4 + 
     -             beta**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4))*z)
     -        + s**2*(-2*(1 + 2*C33phiq3*vol2 + C33phiu*vol2) + 
     -          beta**6*(1 + 5*(-2*C33phiq3 + C33phiu)**2*vol4) + 
     -          4*beta**5*(-2*C33phiq3 + C33phiu)**2*vol4*z + 
     -          2*beta**3*(1 - 20*C33phiq3**2*vol4 + 
     -             4*C33phiq3*C33phiu*vol4 - 5*C33phiu**2*vol4)*z + 
     -          2*beta*(-1 + (2*C33phiq3 + C33phiu)**2*vol4)*z + 
     -          beta**4*(-3 - 52*C33phiq3**2*vol4 + 
     -             12*C33phiq3*C33phiu*vol4 - 13*C33phiu**2*vol4 - 
     -             (-1 + 4*C33phiq3*vol2 + 2*C33phiu*vol2 + 
     -                (-2*C33phiq3 + C33phiu)**2*vol4)*z**2) + 
     -          beta**2*(4 - z**2 - 4*C33phiq3**2*vol4*(2 + z**2) + 
     -             4*C33phiq3*vol2*
     -              (1 + z**2 - C33phiu*vol2*(-2 + z**2)) + 
     -             C33phiu*(2*vol2*(1 + z**2) - 
     -                C33phiu*vol4*(2 + z**2)))))*
     -     xI3(0d0,0d0,s,MT**2,MT**2,MT**2,musq,ep))/
     -   (3.*MW**2*(-1 + beta**2*z**2)) + 
     -  (2*MZ**2*(-7 + 9*beta*z)*
     -     (8*MZ**4*(-1 - 12*C33phiq3**2*vol4 + 
     -          (4*C33phiq3 - 3*C33phiu)*C33phiu*vol4 + 
     -          beta**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4))*
     -        (1 + beta*z)**2 + 
     -       4*beta*s*(MZ + beta*MZ*z)**2*
     -        (beta*(-1 - 36*C33phiq3**2*vol4 + 
     -             3*(4*C33phiq3 - 3*C33phiu)*C33phiu*vol4 + 
     -             beta**2*(1 + 3*(-2*C33phiq3 + C33phiu)**2*vol4)) + 
     -          (-1 - 12*C33phiq3**2*vol4 + 
     -             (4*C33phiq3 - 3*C33phiu)*C33phiu*vol4 + 
     -             beta**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4))*z)
     -        + s**2*(-2 + (2*C33phiq3 + C33phiu)*
     -           (-vol2 + (2*C33phiq3 + C33phiu)*vol4) + 
     -          beta*(-5 + (2*C33phiq3 + C33phiu)*
     -              (-vol2 + 2*(2*C33phiq3 + C33phiu)*vol4))*z + 
     -          4*beta**8*(-2*C33phiq3 + C33phiu)**2*vol4*z**2 + 
     -          beta**7*z*(1 - 2*C33phiq3*vol2 - C33phiu*vol2 + 
     -             40*C33phiq3**2*vol4 - 40*C33phiq3*C33phiu*vol4 + 
     -             10*C33phiu**2*vol4 + 
     -             4*(-2*C33phiq3 + C33phiu)**2*vol4*z**2) + 
     -          beta**3*z*(7 - 4*z**2 + 
     -             8*C33phiq3**2*vol4*(-8 + z**2) + 
     -             2*C33phiu**2*vol4*(-8 + z**2) + 
     -             C33phiu*vol2*(-3 + 4*z**2) + 
     -             2*C33phiq3*vol2*(-3 + 12*C33phiu*vol2 + 4*z**2)) + 
     -          beta**5*z*(-3 + 4*z**2 + C33phiu*vol2*(5 - 4*z**2) - 
     -             8*C33phiq3**2*vol4*(10 + 7*z**2) - 
     -             2*C33phiu**2*vol4*(10 + 7*z**2) + 
     -             2*C33phiq3*vol2*
     -              (5 + 4*C33phiu*vol2 + 
     -                4*(-1 + 2*C33phiu*vol2)*z**2)) + 
     -          beta**2*(4 + 2*C33phiq3*(vol2 + 10*C33phiu*vol4) - 
     -             6*z**2 + 4*C33phiq3**2*vol4*(-2 + z**2) + 
     -             C33phiu*(vol2 + C33phiu*vol4*(-2 + z**2))) + 
     -          beta**4*(-3 + 6*z**2 - z**4 + 
     -             C33phiu*vol2*(-1 + 2*z**4) - 
     -             4*C33phiq3**2*vol4*(15 + 20*z**2 + 2*z**4) - 
     -             C33phiu**2*vol4*(15 + 20*z**2 + 2*z**4) - 
     -             2*C33phiq3*vol2*
     -              (1 - 2*z**4 + 2*C33phiu*vol2*(2 - 10*z**2 + z**4))
     -             ) + beta**6*
     -           (1 + z**4 + 4*C33phiq3**2*vol4*(6 - 5*z**2) + 
     -             C33phiu**2*vol4*(6 - 5*z**2) + 
     -             C33phiu*(vol2 - 2*vol2*z**4) + 
     -             2*C33phiq3*vol2*
     -              (1 - 2*z**4 + 2*C33phiu*vol2*(-4 - 6*z**2 + z**4))
     -             )))*xI3(0d0,MT**2,MT**2 - (s*(1d0 + beta*z))/2d0,MT**2,
     -      MT**2,MZ**2,musq,ep))/
     -   (3.*(-1 + beta*z)*(MW + beta*MW*z)**2) + 
     -  (2*MZ**2*(-7 + 9*beta*z)*
     -     (2*beta**8*(-2*C33phiq3 + C33phiu)**2*s**2*vol4*z*
     -        (-4 + z**2) - 2*beta*s*
     -        (MZ + (2*C33phiq3 + C33phiu)*MZ*vol2)**2*(-1 + z**2) - 
     -       2*beta**7*(-2*C33phiq3 + C33phiu)**2*s**2*vol4*
     -        (1 + z**2) + 2*MZ**4*
     -        (1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*z*(-5 + 3*z**2)
     -        + 4*beta**4*z*(MZ**4*
     -           (1 + (-2*C33phiq3 + C33phiu)**2*vol4)*(-4 + z**2) + 
     -          s**2*vol2*(-2*C33phiq3 - C33phiu - 
     -             (4*C33phiq3**2 + 6*C33phiq3*C33phiu + C33phiu**2)*
     -              vol2 + (2*C33phiq3 + C33phiu)**2*vol2*z**2) + 
     -          MZ**2*s*(3 + 5*C33phiu*vol2 + 14*C33phiu**2*vol4 - 
     -             (1 + 3*C33phiu*vol2 + 2*C33phiu**2*vol4)*z**2 - 
     -             8*C33phiq3**2*vol4*(-7 + z**2) - 
     -             2*C33phiq3*vol2*
     -              (-5 + 4*C33phiu*vol2 + (3 + 2*C33phiu*vol2)*z**2))
     -          ) + beta**5*s*
     -        (-2*MZ**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4)*
     -           (1 + z**2) + 
     -          s*(-1 - 2*C33phiu*vol2 + 7*C33phiu**2*vol4 + 
     -             (-1 + 2*C33phiu*vol2 + 3*C33phiu**2*vol4)*z**2 + 
     -             4*C33phiq3**2*vol4*(7 + 3*z**2) + 
     -             4*C33phiq3*vol2*
     -              (-1 + z**2 + C33phiu*vol2*(-3 + z**2)))) + 
     -       beta**3*s*(s*(1 + 2*C33phiq3*vol2 + C33phiu*vol2)*
     -           (1 + z**2 + (2*C33phiq3 + C33phiu)*vol2*(1 - 3*z**2))
     -            + 4*MZ**2*(z**2 + 4*C33phiq3**2*vol4*(2 + z**2) + 
     -             2*C33phiq3*vol2*(-1 - 2*C33phiu*vol2 + z**2) + 
     -             C33phiu*(vol2*(-1 + z**2) + 
     -                C33phiu*vol4*(2 + z**2)))) + 
     -       beta**6*s*z*(2*MZ**2*
     -           (1 + 3*(-2*C33phiq3 + C33phiu)**2*vol4)*(-4 + z**2)
     -           + s*(-1 + 4*C33phiq3**2*vol4*(19 - 2*z**2) + 
     -             C33phiu*(2*vol2 + C33phiu*vol4*(19 - 2*z**2)) - 
     -             4*C33phiq3*vol2*(-1 + C33phiu*vol2*(3 + 2*z**2))))
     -        + beta**2*z*((s + (2*C33phiq3 + C33phiu)*s*vol2)**2 + 
     -          2*MZ**4*(13 + 10*C33phiu*vol2 + 19*C33phiu**2*vol4 - 
     -             (5 + 6*C33phiu*vol2 + 3*C33phiu**2*vol4)*z**2 + 
     -             4*C33phiq3**2*vol4*(19 - 3*z**2) - 
     -             4*C33phiq3*vol2*
     -              (-5 + 3*C33phiu*vol2 + (3 + C33phiu*vol2)*z**2))
     -           + 2*MZ**2*s*
     -           (4*C33phiq3**2*vol4*(-8 + 5*z**2) + 
     -             4*C33phiq3*vol2*
     -              (-5 - 8*C33phiu*vol2 + (3 + 5*C33phiu*vol2)*z**2)
     -              + (1 + C33phiu*vol2)*
     -              (-2 + z**2 + C33phiu*vol2*(-8 + 5*z**2)))))*
     -     xI3(MT**2,MT**2,s,MT**2,MZ**2,MT**2,musq,ep))/
     -   (3.*beta*MW**2*(-1 + beta**2*z**2)) - 
     -  (MZ**2*(-7 + 9*beta*z)*
     -     (16*MZ**6*(-1 - 12*C33phiq3**2*vol4 + 
     -          (4*C33phiq3 - 3*C33phiu)*C33phiu*vol4 + 
     -          beta**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4)) + 
     -       16*beta*MZ**4*s*
     -        (-(beta*(1 + 24*C33phiq3**2*vol4 - 
     -               8*C33phiq3*C33phiu*vol4 + 6*C33phiu**2*vol4)) + 
     -          beta**3*(1 + 2*(-2*C33phiq3 + C33phiu)**2*vol4) - 
     -          (1 + 12*C33phiq3**2*vol4 - 4*C33phiq3*C33phiu*vol4 + 
     -             3*C33phiu**2*vol4)*z + 
     -          beta**2*(1 + (-2*C33phiq3 + C33phiu)**2*vol4)*z) + 
     -       2*MZ**2*s**2*(-2*(1 + 2*C33phiq3*vol2 + C33phiu*vol2) - 
     -          2*beta**3*(1 + 68*C33phiq3**2*vol4 - 
     -             20*C33phiq3*C33phiu*vol4 + 17*C33phiu**2*vol4)*z + 
     -          4*beta**5*(1 + 3*(-2*C33phiq3 + C33phiu)**2*vol4)*z + 
     -          2*beta*(-1 + (2*C33phiq3 + C33phiu)**2*vol4)*z + 
     -          beta**6*(2 + z**2 + 
     -             (-2*C33phiq3 + C33phiu)**2*vol4*(10 + z**2)) + 
     -          beta**2*(3 + 2*C33phiu*vol2 - 5*C33phiu**2*vol4 + 
     -             2*(-1 + C33phiu*(vol2 - 2*C33phiu*vol4))*z**2 + 
     -             4*C33phiq3*vol2*(1 + 3*C33phiu*vol2 + z**2) - 
     -             4*C33phiq3**2*vol4*(5 + 4*z**2)) + 
     -          beta**4*(-3 + z**2 - 12*C33phiq3**2*vol4*(9 + z**2) + 
     -             4*C33phiq3*
     -              (-(vol2*z**2) + C33phiu*vol4*(7 + z**2)) - 
     -             C33phiu*(2*vol2*z**2 + 3*C33phiu*vol4*(9 + z**2))))
     -         + beta*s**3*(-((1 + 2*C33phiq3*vol2 + C33phiu*vol2)**2*
     -             z) + 8*beta**6*(-2*C33phiq3 + C33phiu)**2*vol4*z + 
     -          2*beta**7*(-2*C33phiq3 + C33phiu)**2*vol4*
     -           (2 + z**2) + 
     -          beta**4*z*(2 - 88*C33phiq3**2*vol4 + 
     -             24*C33phiq3*C33phiu*vol4 - 22*C33phiu**2*vol4 - 
     -             (-1 + 4*C33phiq3*vol2 + 2*C33phiu*vol2 + 
     -                (-2*C33phiq3 + C33phiu)**2*vol4)*z**2) + 
     -          beta**5*(-2*(12*C33phiq3**2 + 4*C33phiq3*C33phiu + 
     -                3*C33phiu**2)*vol4 - 
     -             (-1 + 2*C33phiu*vol2 + 20*C33phiq3**2*vol4 + 
     -                5*C33phiu**2*vol4 + 
     -                4*C33phiq3*(vol2 - C33phiu*vol4))*z**2) + 
     -          beta*(-1 - 2*C33phiu*vol2 + C33phiu**2*vol4 + 
     -             2*(-1 + C33phiu**2*vol4)*z**2 + 
     -             4*C33phiq3*vol2*
     -              (-1 + C33phiu*vol2 + 2*C33phiu*vol2*z**2) + 
     -             4*C33phiq3**2*(vol4 + 2*vol4*z**2)) + 
     -          beta**2*z*(-1 + 2*C33phiu*vol2 + C33phiu**2*vol4 - 
     -             (1 - 2*C33phiu*vol2 + C33phiu**2*vol4)*z**2 - 
     -             4*C33phiq3**2*vol4*(-1 + z**2) + 
     -             4*C33phiq3*vol2*
     -              (1 + z**2 - C33phiu*vol2*(-3 + z**2))) + 
     -          beta**3*(1 + z**2 - 4*C33phiq3**2*vol4*(9 + 5*z**2) + 
     -             4*C33phiq3*vol2*
     -              (1 + z**2 - C33phiu*vol2*(-5 + z**2)) + 
     -             C33phiu*(2*vol2*(1 + z**2) - 
     -                C33phiu*vol4*(9 + 5*z**2)))))*
     -     xI4(0d0,0d0,MT**2,MT**2,s,MT**2 - (s*(1d0 + beta*z))/2d0,MT**2,
     -      MT**2,MT**2,MZ**2,musq,ep))/(3.*MW**2*(-1 + beta**2*z**2)) )     
    
!         print *, "MCFM bx(3) ", bx(3)      
!         print *, "MAKU bx(3) ", MARKUS_bx3
!         pause             

              bx(3) = MARKUS_bx3


!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
!       bx(4) = 
!      .     (0.25_dp*alpha*genfacbox*(0.125_dp/sw2)*sigma0*
!      .     (-0.25_dp*beta*(4._dp*mb**2 + (1._dp - beta**2)*s)*z*(1._dp 
!      .     - z**2)*1 + (beta*(4._dp*mb**2 + (1._dp - beta**2)*s)*
!      .     (1._dp - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
!      .     (xI1(mb**2,musq,ep)*1 - xI1(mw**2,musq,ep)*1))/(s*(1._dp 
!      .     + beta**2 + 2._dp*beta*z)) - (0.125_dp*((1._dp - beta**2
!      .     )**2*s**2*(8._dp*beta - 2._dp*beta**3 - (5._dp 
!      .     - 4._dp*beta**2 - 4._dp*beta**4)*z - 9._dp*beta*z**2 + (3._dp 
!      .     - 4._dp*beta**2 - 2._dp*beta**4)*z**3 + 3._dp*beta*z**4) 
!      .     + 16._dp*beta*(1._dp - beta**2)*mb**2*s*(2._dp 
!      .     - 3._dp*beta**2 + beta*(3._dp - 2._dp*beta**2)*z 
!      .     - (1._dp - 3._dp*beta**2)*z**2 - beta*(2._dp 
!      .     - beta**2)*z**3 - beta**2*z**4) + (16._dp*mb**2*
!      .     (mb**2 - mw**2) - 4._dp*(1._dp - beta**2)*mw**2*s)*
!      .     (6._dp*beta**3 + (5._dp - 8._dp*beta**2 + 4._dp*beta**4)*z 
!      .     + beta*(5._dp - 12._dp*beta**2)*z**2 - (3._dp 
!      .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 - beta*(3._dp 
!      .     - 4._dp*beta**2)*z**4))*xI2(mt**2,mb**2,mw**2,musq,ep)*1)/
!      .     (beta*s*(1._dp + beta*z)) + (0.125_dp*((16._dp*mb**2*(mb**2 
!      .     - mw**2) - 4._dp*(1._dp - beta**2)*mw**2*s)*z*(5._dp 
!      .     - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
!      .     + 8._dp*beta*mb**2*s*(2._dp + beta*(5._dp - 4._dp*beta**2)*z 
!      .     - 2._dp*z**2 - beta*(3._dp - 2._dp*beta**2)*z**3) + (1._dp 
!      .     - beta**2)*s**2*(4._dp*beta - (5._dp - 3._dp*beta**2 
!      .     - 4._dp*beta**4)*z - 4._dp*beta*z**2 + (3._dp - 3._dp*beta**2 
!      .     - 2._dp*beta**4)*z**3))*xI2(s,mb**2,mb**2,musq,ep)*1)/(beta*s) 
!      .     - (0.25_dp*(beta*(16._dp*mb**2*(mb**2 - mw**2) 
!      .     - 4._dp*(1._dp - beta**2)*mw**2*s)*(beta + z)*(1._dp 
!      .     - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2)*z 
!      .     + 2._dp*beta**2*z**2 + beta**3*z**3) - (1._dp 
!      .     - beta**2)*s**2*(1._dp + beta**2 + 2._dp*beta*z)*
!      .     (2._dp - 5._dp*beta**2 + beta**4 + beta*(1._dp - 2._dp*beta**2 
!      .     - 2._dp*beta**4)*z + beta**2*(3._dp - 2._dp*beta**2)*z**2 
!      .     + beta**3*(2._dp + beta**2)*z**3 + beta**4*z**4) 
!      .     - 8._dp*mb**2*s*(1._dp - 6._dp*beta**2 + 3._dp*beta**6 
!      .     + beta*(2._dp - 16._dp*beta**2 + 7._dp*beta**4 + 2._dp*beta**6)*z 
!      .     + beta**2*(4._dp - 9._dp*beta**2 + 3._dp*beta**4)*z**2 
!      .     + beta**3*(9._dp - 4._dp*beta**2 - beta**4)*z**3 
!      .     + 2._dp*beta**4*(3._dp - beta**2)*z**4 + beta**5*z**5))*
!      .     xI2(t,mb**2,mw**2,musq,ep)*1)/(s*(1._dp + beta*z)*(1._dp 
!      .     + beta**2 + 2._dp*beta*z)) + (0.125_dp*(64._dp*mb**6 
!      .     + 16._dp*mw**4*(4._dp*mb**2 + (1._dp - beta**2)*s) 
!      .     + 8._dp*beta*(1._dp - beta**2)*mw**2*s**2*(beta + z) 
!      .     + 32._dp*mb**2*mw**2*s*(1._dp + beta*z) + 16._dp*mb**4*
!      .     (-8._dp*mw**2 + (1._dp - beta**2)*s - 2._dp*beta*s*z) 
!      .     + 4._dp*mb**2*s**2*(3._dp - 2._dp*beta**2 + beta**4 
!      .     + 4._dp*beta*(2._dp - beta**2)*z + 2._dp*beta**2*z**2) 
!      .     + (1._dp - beta**2)*s**3*(3._dp - 2._dp*beta**2 + beta**4 
!      .     + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2))*
!      .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep)*1)/s 
!      .     - (0.125_dp*(64._dp*mb**6*(1._dp + beta*z) 
!      .     + 16._dp*(1._dp - beta**2)*mw**4*s*(1._dp + beta*z) 
!      .     + 8._dp*beta*(1._dp - beta**2)*mw**2*s**2*(beta + z)*
!      .     (1._dp + beta*z) + (1._dp - beta**2)*s**3*(1._dp + beta*z)*
!      .     (3._dp - 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp + beta**2)*z 
!      .     + 2._dp*beta**2*z**2) - 16._dp*mb**4*(8._dp*mw**2*(1._dp 
!      .     + beta*z) - s*(1._dp - beta**2 - beta*(3._dp 
!      .     - beta**2)*z - 2._dp*beta**2*z**2)) + mb**2*(4._dp*s**2*
!      .     (3._dp - 2._dp*beta**2 + beta**4 + 9._dp*beta*z - 2._dp*beta**3*z 
!      .     - beta**5*z + 10._dp*beta**2*z**2 - 4._dp*beta**4*z**2 
!      .     + 2._dp*beta**3*z**3) + 32._dp*(1._dp + beta*z)*(2._dp*mw**4 
!      .     + mw**2*s*(1._dp + beta*z))))*xI3(0._dp,mt**2,t,mb**2,mb**2
!      .     ,mw**2,musq,ep)*1)/s + (0.03125_dp*(-16._dp*(4._dp*mb**6 
!      .     - 8._dp*mb**4*mw**2 + mw**4*(4._dp*mb**2 + (1._dp - beta**2
!      .     )*s))*z*(5._dp - 8._dp*beta**2 - (3._dp - 2._dp*beta**2
!      .     )*z**2) + 64._dp*beta*mb**2*mw**2*s*(1._dp + beta**2 
!      .     + 2._dp*beta*z - (1._dp - beta**2)*z**2) 
!      .     + 8._dp*(1._dp - beta**2)*mw**2*s**2*(2._dp*beta*(1._dp 
!      .     + beta**2) - (5._dp - beta**2 - 8._dp*beta**4)*z 
!      .     - 2._dp*beta*(1._dp - beta**2)*z**2 + (3._dp - beta**2 
!      .     - 2._dp*beta**4)*z**3) - 16._dp*mb**4*s*(4._dp*beta*(1._dp 
!      .     + beta**2) - (5._dp - 17._dp*beta**2 + 8._dp*beta**4)*z 
!      .     - 4._dp*beta*(1._dp - beta**2)*z**2 + (3._dp - 9._dp*beta**2 
!      .     + 2._dp*beta**4)*z**3) + 4._dp*mb**2*s**2*(8._dp*beta**3*(2._dp 
!      .     - beta**2) + (5._dp - 2._dp*beta**2 + 21._dp*beta**4 
!      .     - 8._dp*beta**6)*z + 8._dp*beta**3*(2._dp - beta**2)*z**2 
!      .     - (3._dp - beta**4 - 2._dp*beta**6)*z**3) + (1._dp 
!      .     - beta**2)*s**3*(4._dp*beta*(1._dp + beta**4) - (5._dp 
!      .     - 22._dp*beta**2 + 9._dp*beta**4 - 8._dp*beta**6)*z 
!      .     - 4._dp*beta*(1._dp - 2._dp*beta**2 - beta**4)*z**2 
!      .     + (3._dp - 4._dp*beta**2 + 3._dp*beta**4 - 2._dp*beta**6)*z**3))*
!      .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)*1)/(beta*s) 
!      .     - (0.03125_dp*(256._dp*mb**8 - 64._dp*(1._dp - beta**2)*
!      .     mw**6*s - 16._dp*(1._dp - beta**2)*mw**4*s**2*(1._dp 
!      .     + 3._dp*beta**2 + 4._dp*beta*z) - (1._dp - beta**2)*
!      .     s**4*(1._dp + beta**2 + 2._dp*beta*z)*(3._dp - 2._dp*beta**2 
!      .     + beta**4 + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2) 
!      .     - 4._dp*(1._dp - beta**2)*mw**2*s**3*(3._dp + 3._dp*beta**4 
!      .     + 4._dp*beta*(1._dp + 2._dp*beta**2)*z + 6._dp*beta**2*z**2) 
!      .     - 128._dp*mb**6*(6._dp*mw**2 + beta*s*z*(2._dp + beta*z)) 
!      .     + 32._dp*mb**4*(24._dp*mw**4 + 2._dp*mw**2*s*(3._dp + beta**2 
!      .     + 8._dp*beta*z + 2._dp*beta**2*z**2) - s**2*(1._dp 
!      .     - 3._dp*beta**2 + beta**4 - 6._dp*beta*z + 2._dp*beta**3*z 
!      .     - 5._dp*beta**2*z**2 + 2._dp*beta**4*z**2)) - 8._dp*mb**2*
!      .     (32._dp*mw**6 + 16._dp*mw**4*s*(1._dp + beta**2 + 2._dp*beta*z) 
!      .     + 2._dp*mw**2*s**2*(5._dp + beta**4 + 12._dp*beta*z 
!      .     + 2._dp*beta**2*(2._dp + beta**2)*z**2) + s**3*(-2._dp 
!      .     + 6._dp*beta**2 - 2._dp*beta**4 + 4._dp*beta*z + 4._dp*beta**3*z 
!      .     - 2._dp*beta**5*z + 9._dp*beta**2*z**2 - 4._dp*beta**4*z**2 
!      .     + beta**6*z**2 + 2._dp*beta**3*z**3)))*xI4(0._dp,0._dp,mt**2,
!      .     mt**2,s,t,mb**2,mb**2,mb**2,mw**2,musq,ep)*1)/s))/(mw**2*pi)



         MARKUS_bx4= alpha/(64d0*Pi*SW2*MW**2)/8d0*  (
     -(-4*s*(-7 + 9*beta*z)*
     -(6*C33phiq3**2*vol4 - 
     -beta*z*(-1 + z**2 + C33phiq3**2*vol4*(-27 + z**2) + 
     -  2*C33phiq3*vol2*(7 + z**2)) + 
     -2*beta**5*C33phiq3*vol2*z*
     -(-2*(1 + z**2) + C33phiq3*vol2*(2 + 3*z**2)) + 
     -beta**3*z*(-1 + z**2 + 6*C33phiq3*vol2*(3 + z**2) - 
     -  3*C33phiq3**2*vol4*(5 + 7*z**2)) - 
     -beta**2*(-z**2 + z**4 + 
     -  C33phiq3**2*vol4*(8 - 13*z**2 + z**4) + 
     -  2*C33phiq3*vol2*(2 + z**2 + z**4)) + 
     -beta**4*(-z**2 + z**4 + 
     -  2*C33phiq3*vol2*(2 + z**2 + z**4) + 
     -  C33phiq3**2*vol4*(4 - 17*z**2 + 3*z**4))))/
     -(3.*(-1 + beta*z)*(1 + beta*z)**2) + 
     -(16*(-7 + 9*beta*z)*(-4*beta**2*(-1 + beta**2)*
     -(1 + C33phiq3*vol2 + beta**2*(-1 + 2*C33phiq3*vol2)) + 
     -2*(1 - 2*beta**2 - 14*beta**4 + 7*beta**6)*C33phiq3**2*
     -vol4 + beta*((-1 + beta**2)**2*(-1 + 4*beta**2) - 
     -  2*(1 - 9*beta**4 + 8*beta**6)*C33phiq3*vol2 + 
     -  (1 - 18*beta**2 - 13*beta**4 + 14*beta**6)*
     -   C33phiq3**2*vol4)*z - 
     -beta**2*(-1 + beta**2)*
     -(-7 - 5*C33phiq3*(2*vol2 + C33phiq3*vol4) + 
     -  beta**2*(7 - 14*C33phiq3*vol2 + 25*C33phiq3**2*vol4))*
     -z**2 - beta*(-1 + beta**2)*
     -(1 + 2*C33phiq3*vol2 + C33phiq3**2*vol4 + 
     -  6*beta**4*(1 - 4*C33phiq3*vol2 + 
     -     4*C33phiq3**2*vol4) + 
     -  beta**2*(-7 - 2*C33phiq3*vol2 + 11*C33phiq3**2*vol4))*
     -z**3 + beta**2*
     -(3*(-1 + beta**2)**2 - 6*(-1 + beta**4)*C33phiq3*vol2 + 
     -  (3 + 4*beta**2 + 9*beta**4)*C33phiq3**2*vol4)*z**4 + 
     -2*beta**3*(-((-1 + beta**2)*
     -     (1 + 2*C33phiq3*vol2 + 
     -       beta**2*(-1 + 4*C33phiq3*vol2))) + 
     -  (1 + 3*beta**2 + 4*beta**4)*C33phiq3**2*vol4)*z**5)*
     -xI1(MW**2,musq,ep))/
     -(3.*(-1 + beta**2)*(-1 + beta*z)*(1 + beta*z)**2*
     -(1 + beta**2 + 2*beta*z)) + 
     -(2*(-7 + 9*beta*z)*(2*beta**11*s*(-1 + C33phiq3*vol2)**2*z**2*
     -(-2 + z**2) + (4*MW**2 + s)*(1 + C33phiq3*vol2)**2*z*
     -(-5 + 3*z**2) + 
     -2*beta*(1 + C33phiq3*vol2)*
     -(4*s - (20*MW**2 + 7*s)*(1 + C33phiq3*vol2)*z**2 + 
     -  3*(4*MW**2 + s)*(1 + C33phiq3*vol2)*z**4) + 
     -beta**2*z*(72*MW**2 + 27*s + 96*C33phiq3*MW**2*vol2 + 
     -  18*C33phiq3*s*vol2 + 40*C33phiq3**2*MW**2*vol4 - 
     -  5*C33phiq3**2*s*vol4 - 
     -  2*(s*(11 + 10*C33phiq3*vol2 + C33phiq3**2*vol4) + 
     -     MW**2*(30 + 22*C33phiq3*(2*vol2 + C33phiq3*vol4)))*
     -   z**2 + 3*(4*MW**2 + s)*
     -   (1 + 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z**4) + 
     -2*beta**10*s*z*(-1 + z**2 + 
     -  C33phiq3**2*vol4*(-3 + 7*z**2 - 2*z**4) + 
     -  2*C33phiq3*vol2*(2 - 4*z**2 + z**4)) + 
     -2*beta**9*(4*MW**2*(-1 + C33phiq3*vol2)**2*z**2*
     -   (-2 + z**2) + 
     -  s*(1 + 4*C33phiq3*vol2 - C33phiq3**2*vol4 + 
     -     2*(2 - 6*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*
     -      z**2 - (1 + C33phiq3**2*vol4)*z**4)) + 
     -beta**4*z*(4*MW**2*
     -   (-31 - 6*C33phiq3*vol2 - 19*C33phiq3**2*vol4 + 
     -     (35 + 10*C33phiq3*vol2 + 31*C33phiq3**2*vol4)*
     -      z**2 - 2*
     -      (5 + 6*C33phiq3*vol2 + 3*C33phiq3**2*vol4)*z**4)
     -   + s*(-49 - 6*C33phiq3*vol2 + 11*C33phiq3**2*vol4 + 
     -     2*(23 + 8*C33phiq3*vol2 - C33phiq3**2*vol4)*z**2 + 
     -     (-9 + 6*C33phiq3*vol2 + 11*C33phiq3**2*vol4)*z**4))
     - + beta**6*z*(4*MW**2*
     -   (4*(7 - 8*C33phiq3*vol2 + 11*C33phiq3**2*vol4) + 
     -     (-37 + 46*C33phiq3*vol2 - 49*C33phiq3**2*vol4)*
     -       z**2 + (11 - 6*C33phiq3*vol2 - 
     -         5*C33phiq3**2*vol4)*z**4) + 
     -   s*(35 + 6*C33phiq3*vol2 - 41*C33phiq3**2*vol4 + 
     -      12*(-3 - C33phiq3*vol2 + C33phiq3**2*vol4)*z**2 + 
     -      (9 + 13*C33phiq3*(-2*vol2 + C33phiq3*vol4))*z**4))
     - + 2*beta**5*(s*
     -    (15 + 28*C33phiq3*vol2 + 5*C33phiq3**2*vol4 + 
     -      (-25 - 22*C33phiq3*vol2 + C33phiq3**2*vol4)*z**2 - 
     -      2*(-7 + 13*C33phiq3*vol2 + C33phiq3**2*vol4)*z**4)
     -    + 4*MW**2*(6 - 27*z**2 + 12*z**4 + 
     -      C33phiq3**2*vol4*(10 - 19*z**2) - 
     -      2*C33phiq3*vol2*(9 - 11*z**2 + z**4))) + 
     -beta**8*z*(-8*MW**2*
     -    (5 - 7*z**2 + 2*z**4 - 
     -      2*C33phiq3*vol2*(6 - 10*z**2 + 3*z**4) + 
     -      C33phiq3**2*vol4*(7 - 13*z**2 + 4*z**4)) + 
     -   s*(-6 + 7*z**2 - 3*z**4 + 
     -      C33phiq3**2*vol4*(46 - 25*z**2 - 23*z**4) + 
     -      2*C33phiq3*vol2*(-8 + 13*z**2 + 5*z**4))) + 
     -2*beta**3*(4*MW**2*
     -    (-3 + 6*C33phiq3*vol2 + 3*C33phiq3**2*vol4 + 
     -      2*(10 + 6*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*
     -       z**2 - 2*
     -       (5 + 6*C33phiq3*vol2 + 3*C33phiq3**2*vol4)*z**4)
     -    + s*(-13 + 23*z**2 - 11*z**4 + 
     -      2*C33phiq3*vol2*(-9 + 15*z**2 + z**4) + 
     -      C33phiq3**2*vol4*(-15 + 13*z**2 + 9*z**4))) + 
     -2*beta**7*(4*MW**2*
     -    (-3 + 14*z**2 - 6*z**4 + 
     -      C33phiq3**2*vol4*(-5 + 16*z**2 - 6*z**4) + 
     -      2*C33phiq3*vol2*(6 - 14*z**2 + 5*z**4)) + 
     -   s*(-7 + 7*z**2 - 6*z**4 + 
     -      2*C33phiq3*vol2*(-9 + 7*z**2 + 10*z**4) + 
     -      C33phiq3**2*vol4*(11 - 5*z**2*(3 + 2*z**2)))))*
     -xI2(MT**2,0d0,MW**2,musq,ep))/
     -(3.*beta*(-1 + beta**2)*(-1 + beta*z)*(1 + beta*z)**3) + 
     -(2*(-7 + 9*beta*z)*(2*beta**6*s*(-1 + C33phiq3*vol2)**2*z*
     - (-2 + z**2) - 4*beta*s*(1 + C33phiq3*vol2)**2*
     - (-1 + z**2) - 4*beta**3*s*
     - (-1 - 2*C33phiq3*vol2 + 3*C33phiq3**2*vol4)*(-1 + z**2)
     - + (4*MW**2 + s)*(1 + C33phiq3*vol2)**2*z*
     - (-5 + 3*z**2) + 
     -2*beta**2*z*(4*s - 
     -   2*C33phiq3*s*(4*vol2 + 3*C33phiq3*vol4) + 
     -   s*(-3 + 6*C33phiq3*vol2 + 7*C33phiq3**2*vol4)*z**2 + 
     -   2*MW**2*(9 + 2*C33phiq3*vol2 + 5*C33phiq3**2*vol4 - 
     -      (5 + 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z**2)) + 
     -beta**4*z*(8*MW**2*(-1 + C33phiq3*vol2)**2*(-2 + z**2) + 
     -   s*(1 + z**2 + 2*C33phiq3*vol2*(9 - 7*z**2) + 
     -      C33phiq3**2*vol4*(-3 + 5*z**2))))*
     -xI2(s,0d0,0d0,musq,ep))/(3.*beta*(-1 + beta**2*z**2)) - 
     -(4*(-7 + 9*beta*z)*(-2*(s + 4*C33phiq3**2*MW**2*vol4) - 
     -beta*(4*MW**2*(1 + 4*C33phiq3*vol2 + 
     -      3*C33phiq3**2*vol4) + s*(7 + 5*C33phiq3**2*vol4))*z
     -  + beta**10*s*(-1 + C33phiq3*vol2)**2*z**2*
     - (-2 + z**2) + beta**9*s*z*
     - (-1 - 5*z**2 + 3*z**4 + 
     -   4*C33phiq3*vol2*(1 + z**2 - z**4) + 
     -   C33phiq3**2*vol4*(-3 + z**2 + z**4)) + 
     -beta**3*z*(s*(21 + 59*C33phiq3**2*vol4 - 
     -      (13 + 36*C33phiq3*vol2 + 59*C33phiq3**2*vol4)*z**2)
     -     -4*MW**2*(-1 + C33phiq3*vol2)*
     -    (2 - 3*z**2 + C33phiq3*vol2*(-18 + 13*z**2))) + 
     -beta**4*(s*(1 + 24*z**2 - 13*z**4 - 
     -      4*C33phiq3*(vol2 + 8*vol2*z**2) + 
     -      C33phiq3**2*vol4*(31 + 10*z**2 - 17*z**4)) + 
     -   4*MW**2*(4 + 4*z**2 - 3*z**4 + 
     -      C33phiq3**2*vol4*(14 + 12*z**2 - 27*z**4) + 
     -      2*C33phiq3*vol2*(-1 + 4*z**4))) + 
     -beta**5*z*(4*MW**2*
     -    (4 + 2*z**2 - z**4 + 4*C33phiq3*vol2*(3 - 2*z**2) + 
     -      C33phiq3**2*vol4*(24 - 26*z**2 + z**4)) + 
     -   s*(-8 + 17*z**2 - 7*z**4 + 
     -      C33phiq3**2*vol4*(48 - 21*z**2 - 13*z**4) + 
     -      4*C33phiq3*vol2*(-6 + 5*z**2 + 4*z**4))) + 
     -beta**8*(4*MW**2*(-1 + C33phiq3*vol2)**2*z**2*
     -    (-2 + z**2) + 
     -   s*(1 - 6*z**2 + z**4 + 2*z**6 + 
     -      4*C33phiq3*vol2*(1 + 4*z**2 - 4*z**4) + 
     -      C33phiq3**2*vol4*(-1 - 6*z**2 + 7*z**4))) + 
     -beta**6*(s*(-5 - 6*z**2 + 11*z**4 - 2*z**6 + 
     -      2*C33phiq3*vol2*(-3 + 14*z**2 + 9*z**4) + 
     -      C33phiq3**2*vol4*(7 + 30*z**2 - 35*z**4 - 12*z**6))
     -     + 4*MW**2*(-3 + 2*z**4 - 
     -      2*C33phiq3*vol2*(-6 + 2*z**2 + 3*z**4) + 
     -      C33phiq3**2*vol4*(-5 + 8*z**2 - 4*z**4 + 2*z**6)))
     - + beta**7*z*(-4*MW**2*
     -    (5 - z**2*(1 + z**2) + 
     -      4*C33phiq3*vol2*(-3 + 2*z**2) + 
     -      C33phiq3**2*vol4*(7 - 7*z**2 + z**4)) + 
     -   s*(-5 + z**2 + 4*z**4 + 
     -      4*C33phiq3*vol2*(5 + 3*z**2 - 3*z**4) + 
     -      C33phiq3**2*vol4*(9 - 5*z**2 - 16*z**4 + 4*z**6)))
     - + beta**2*(s*(5 + 6*C33phiq3*vol2 + 
     -      23*C33phiq3**2*vol4 - 
     -      2*(5 + 8*C33phiq3*vol2 + 18*C33phiq3**2*vol4)*z**2)
     -     + 4*MW**2*(-1 - 2*z**2 + 
     -      C33phiq3*(-10*vol2 + C33phiq3*vol4*(-1 + 4*z**2))))
     -)*xI2(MT**2 - (s*(1d0+beta*z))/2d0,0d0,MW**2,musq,ep))/
     -(3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z))
     -- (2*(-7 + 9*beta*z)*
     -(16*MW**4*(-1 + 2*C33phiq3*vol2 + 
     -   beta**2*(-1 + C33phiq3*vol2)**2 - 5*C33phiq3**2*vol4)
     - + 8*MW**2*s*(beta**4*(-1 + C33phiq3*vol2)**2 - 
     -   8*C33phiq3**2*vol4 + 
     -   beta**2*(-1 + 2*C33phiq3*vol2 - 5*C33phiq3**2*vol4) + 
     -   beta**3*(-1 + C33phiq3*vol2)**2*z + 
     -   beta*(-1 + 2*C33phiq3*vol2 - 5*C33phiq3**2*vol4)*z) + 
     -s**2*(-3 + beta**6*(-1 + C33phiq3*vol2)**2 - 
     -   5*C33phiq3*(2*vol2 + 3*C33phiq3*vol4) + 
     -   2*beta**5*(-1 + C33phiq3*vol2)**2*z - 
     -   8*beta**3*C33phiq3**2*vol4*z - 
     -   2*beta*(1 - 2*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*z + 
     -   beta**4*(-3 + 6*C33phiq3*vol2 - 7*C33phiq3**2*vol4 + 
     -      2*(1 - 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z**2) + 
     -   beta**2*(5 + 6*C33phiq3*vol2 - 19*C33phiq3**2*vol4 - 
     -      2*(1 - 2*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*z**2))
     -)*xI3(0d0,0d0,s,0d0,0d0,0d0,musq,ep))/(-3 + 3*beta**2*z**2)+ 
     -(2*(-7 + 9*beta*z)*(16*MW**4*
     - (-1 + 2*C33phiq3*vol2 + 
     -   beta**2*(-1 + C33phiq3*vol2)**2 - 5*C33phiq3**2*vol4)
     - + 8*MW**2*s*(beta**4*(-1 + C33phiq3*vol2)**2 - 
     -   8*C33phiq3**2*vol4 + 
     -   beta**2*(-1 + 2*C33phiq3*vol2 - 5*C33phiq3**2*vol4) + 
     -   beta**3*(-1 + C33phiq3*vol2)**2*z + 
     -   beta*(-1 + 2*C33phiq3*vol2 - 5*C33phiq3**2*vol4)*z) + 
     -s**2*(-3 + beta**6*(-1 + C33phiq3*vol2)**2 - 
     -   5*C33phiq3*(2*vol2 + 3*C33phiq3*vol4) + 
     -   2*beta**5*(-1 + C33phiq3*vol2)**2*z - 
     -   8*beta**3*C33phiq3**2*vol4*z - 
     -   2*beta*(1 - 2*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*z + 
     -   beta**4*(-3 + 6*C33phiq3*vol2 - 7*C33phiq3**2*vol4 + 
     -      2*(1 - 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z**2) + 
     -   beta**2*(5 + 6*C33phiq3*vol2 - 19*C33phiq3**2*vol4 - 
     -      2*(1 - 2*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*z**2))
     -)*xI3(0d0,MT**2,MT**2-(s*(1d0+beta*z))/2d0,0d0,0d0,MW**2,musq,
     -ep))/(-3 + 3*beta*z) + 
     -((-7 + 9*beta*z)*(2*beta**8*s**2*(-1 + C33phiq3*vol2)**2*z*
     - (-4 + z**2) - 4*beta*s*(4*MW**2 + s)*
     - (1 + C33phiq3*vol2)**2*(-1 + z**2) - 
     -4*beta**7*s**2*(-1 + C33phiq3*vol2)**2*(1 + z**2) + 
     -(4*MW**2 + s)**2*(1 + C33phiq3*vol2)**2*z*
     - (-5 + 3*z**2) - 
     -4*beta**5*s*(4*MW**2*(-1 + C33phiq3*vol2)**2*
     -    (1 + z**2) + 
     -   s*(-1 + z**2 + C33phiq3**2*vol4*(-9 + z**2) - 
     -      6*C33phiq3*vol2*(-1 + z**2))) + 
     -beta**2*z*(s**2*
     -    (27 + 38*C33phiq3*vol2 + 55*C33phiq3**2*vol4 + 
     -      (-7 + 10*C33phiq3*vol2 + 13*C33phiq3**2*vol4)*z**2)
     -     + 16*MW**2*s*
     -    (3 - 6*C33phiq3*vol2 + 13*C33phiq3**2*vol4 + 
     -      2*(-1 + 2*C33phiq3*(vol2 + C33phiq3*vol4))*z**2) - 
     -   16*MW**4*(-13 + 5*z**2 + 
     -      C33phiq3**2*vol4*(-25 + z**2) + 
     -      2*C33phiq3*vol2*(3 + z**2))) - 
     -beta**6*s*z*(-16*MW**2*(-1 + C33phiq3*vol2)**2*
     -    (-4 + z**2) + 
     -   s*(-17 + 5*z**2 + C33phiq3**2*vol4*(-29 + z**2) + 
     -      2*C33phiq3*vol2*(7 + z**2))) + 
     -beta**4*z*(32*MW**4*(-1 + C33phiq3*vol2)**2*
     -    (-4 + z**2) + 
     -   8*MW**2*s*(7 + 6*C33phiq3*vol2 + 
     -      19*C33phiq3**2*vol4 + 
     -      (-1 - 10*C33phiq3*vol2 + 3*C33phiq3**2*vol4)*z**2)
     -    + s**2*(-31 + 7*z**2 - 10*C33phiq3*vol2*(3 + z**2) + 
     -      C33phiq3**2*vol4*(57 + 47*z**2))) - 
     -4*beta**3*s*(-8*MW**2*
     -    (-2*C33phiq3*vol2 + z**2 + 
     -      C33phiq3**2*vol4*(4 + z**2)) + 
     -   s*(1 - 3*z**2 + 
     -      C33phiq3*(5*C33phiq3*vol4*(-3 + z**2) + 
     -         vol2*(-2 + 6*z**2)))))*
     -xI3(MT**2,MT**2,s,0d0,MW**2,0d0,musq,ep))/
     -(6.*beta*(-1 + beta**2*z**2)) - 
     -((-7 + 9*beta*z)*(4*MW**2 + s*(1 + beta**2 + 2*beta*z))*
     -(16*MW**4*(-1 + 2*C33phiq3*vol2 + 
     -   beta**2*(-1 + C33phiq3*vol2)**2 - 5*C33phiq3**2*vol4)
     - + 8*MW**2*s*(beta**4*(-1 + C33phiq3*vol2)**2 - 
     -   8*C33phiq3**2*vol4 + 
     -   beta**2*(-1 + 2*C33phiq3*vol2 - 5*C33phiq3**2*vol4) + 
     -   beta**3*(-1 + C33phiq3*vol2)**2*z + 
     -   beta*(-1 + 2*C33phiq3*vol2 - 5*C33phiq3**2*vol4)*z) + 
     -s**2*(-3 + beta**6*(-1 + C33phiq3*vol2)**2 - 
     -   5*C33phiq3*(2*vol2 + 3*C33phiq3*vol4) + 
     -   2*beta**5*(-1 + C33phiq3*vol2)**2*z - 
     -   8*beta**3*C33phiq3**2*vol4*z - 
     -   2*beta*(1 - 2*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*z + 
     -   beta**4*(-3 + 6*C33phiq3*vol2 - 7*C33phiq3**2*vol4 + 
     -      2*(1 - 2*C33phiq3*vol2 + C33phiq3**2*vol4)*z**2) + 
     -   beta**2*(5 + 6*C33phiq3*vol2 - 19*C33phiq3**2*vol4 - 
     -      2*(1 - 2*C33phiq3*vol2 + 5*C33phiq3**2*vol4)*z**2))
     -)*xI4(0d0,0d0,MT**2,MT**2,s,MT**2 - (s*(1+beta*z))/2d0,0d0,0d0,
     -0d0,MW**2,musq,ep))/(-6 + 6*beta**2*z**2)   )

!         print *, "MCFM bx(4) ", bx(4)      
!         print *, "MAKU bx(4) ", MARKUS_bx4
!         pause             

              bx(4) = MARKUS_bx4



!       bx(5) = 
!      .     (-0.25_dp*alpha*genfacbox*gw_sq*mt**2*sigma0*
!      .     (beta*z*(1._dp - z**2) - (4._dp*beta*(1._dp - z**2)*
!      .     (4._dp*beta - z - 2._dp*beta*z**2)*(-xI1(mh**2,musq,
!      .     ep) + xI1(mt**2,musq,ep)))/(s*(1._dp + beta**2 + 2._dp*beta*z)) 
!      .     + ((2._dp*beta*(1._dp - beta**2)*s*(2._dp + 3._dp*beta**2 
!      .     - beta*(3._dp - 4._dp*beta**2)*z - (1._dp 
!      .     + 6._dp*beta**2)*z**2 + beta*(1._dp - 2._dp*beta**2)*z**3 
!      .     + 2._dp*beta**2*z**4) + 2._dp*mh**2*(-6._dp*beta**3 - (5._dp 
!      .     - 8._dp*beta**2 + 4._dp*beta**4)*z - beta*(5._dp 
!      .     - 12._dp*beta**2)*z**2 + (3._dp - 4._dp*beta**2 + 2._dp*beta**4
!      .     )*z**3 + beta*(3._dp - 4._dp*beta**2)*z**4))*xI2(mt**2,mt**2,
!      .     mh**2,musq,ep))/(beta*s*(1._dp + beta*z)) + ((2._dp*mh**2*z*
!      .     (5._dp - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
!      .     - beta*s*(2._dp - beta*(7._dp - 8._dp*beta**2)*z 
!      .     - 2._dp*z**2 + beta*(3._dp - 4._dp*beta**2)*z**3))*
!      .     xI2(s,mt**2,mt**2,musq,ep))/(beta*s) + ((-4._dp*beta*mh**2*
!      .     (beta + z)*(1._dp - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2
!      .     )*z + 2._dp*beta**2*z**2 + beta**3*z**3) - 2._dp*s*(1._dp 
!      .     - 3._dp*beta**6 + beta*(2._dp + 2._dp*beta**2 - 5._dp*beta**4 
!      .     - 4._dp*beta**6)*z + beta**2*(1._dp + 3._dp*beta**2 
!      .     - 6._dp*beta**4)*z**2 + 2._dp*beta**5*(1._dp + beta**2)*z**3 
!      .     + 4._dp*beta**6*z**4 + beta**5*z**5))*xI2(t,mt**2,mh**2,musq,
!      .     ep))/(s*(1._dp + beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) - 
!      .     (1._dp*(8._dp*mh**4 - 4._dp*mh**2*s*(2._dp - 3._dp*beta**2 
!      .     - beta*z) + s**2*(6._dp - 10._dp*beta**2 + 5._dp*beta**4 
!      .     - 2._dp*beta*(1._dp - 2._dp*beta**2)*z + beta**2*z**2))*
!      .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/s + ((s**2*(6._dp 
!      .     - 10._dp*beta**2 + 5._dp*beta**4 + beta*(3._dp - 4._dp*beta**2 
!      .     + 4._dp*beta**4)*z - beta**2*(1._dp - 4._dp*beta**2)*z**2 
!      .     + beta**3*z**3) + (1._dp + beta*z)*(8._dp*mh**4 - 4._dp*mh**2*s*
!      .     (2._dp - 3._dp*beta**2 - beta*z)))*xI3(0._dp,mt**2,t,mt**2
!      .     ,mt**2,mh**2,musq,ep))/s + ((2._dp*mh**4*z*(5._dp 
!      .     - 8._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
!      .     - 2._dp*beta*mh**2*s*(1._dp + beta**2 - 2._dp*beta*(5._dp 
!      .     - 6._dp*beta**2)*z - (1._dp - beta**2)*z**2 
!      .     + 3._dp*beta*(1._dp - beta**2)*z**3) + beta**2*s**2*
!      .     (beta - 2._dp*beta**3 - (7._dp - 13._dp*beta**2 
!      .     + 8._dp*beta**4)*z + beta*(1._dp - 2._dp*beta**2)*z**2 
!      .     - 2._dp*beta**2*(1._dp - beta**2)*z**3))*
!      .     xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2,musq,ep))/(beta*s) 
!      .     + (0.5_dp*(-16._dp*mh**6 + 16._dp*mh**4*s*(1._dp - 2._dp*beta**2 
!      .     - beta*z) - 2._dp*mh**2*s**2*(6._dp - 13._dp*beta**2 
!      .     + 10._dp*beta**4 - 6._dp*beta*(1._dp - 2._dp*beta**2)*z 
!      .     + beta**2*(2._dp + beta**2)*z**2) - s**3*(2._dp + beta**2 
!      .     - 6._dp*beta**4 + 4._dp*beta**6 + beta*(5._dp - 10._dp*beta**2 
!      .     + 8._dp*beta**4)*z + beta**4*(1._dp + 2._dp*beta**2)*z**2 
!      .     + beta**3*z**3))*xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,
!      .     mt**2,mh**2,musq,ep))/s))/(mw**2*pi)
! 


!MARKUS: replaced gw_sq coupl. by vev          
      bx(5) = 
     .     (-0.25_dp*genfacbox/vev**2*mt**2*sigma0*
     .     (beta*z*(1._dp - z**2) - (4._dp*beta*(1._dp - z**2)*
     .     (4._dp*beta - z - 2._dp*beta*z**2)*(-xI1(mh**2,musq,
     .     ep) + xI1(mt**2,musq,ep)))/(s*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + ((2._dp*beta*(1._dp - beta**2)*s*(2._dp + 3._dp*beta**2 
     .     - beta*(3._dp - 4._dp*beta**2)*z - (1._dp 
     .     + 6._dp*beta**2)*z**2 + beta*(1._dp - 2._dp*beta**2)*z**3 
     .     + 2._dp*beta**2*z**4) + 2._dp*mh**2*(-6._dp*beta**3 - (5._dp 
     .     - 8._dp*beta**2 + 4._dp*beta**4)*z - beta*(5._dp 
     .     - 12._dp*beta**2)*z**2 + (3._dp - 4._dp*beta**2 + 2._dp*beta**4
     .     )*z**3 + beta*(3._dp - 4._dp*beta**2)*z**4))*xI2(mt**2,mt**2,
     .     mh**2,musq,ep))/(beta*s*(1._dp + beta*z)) + ((2._dp*mh**2*z*
     .     (5._dp - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     - beta*s*(2._dp - beta*(7._dp - 8._dp*beta**2)*z 
     .     - 2._dp*z**2 + beta*(3._dp - 4._dp*beta**2)*z**3))*
     .     xI2(s,mt**2,mt**2,musq,ep))/(beta*s) + ((-4._dp*beta*mh**2*
     .     (beta + z)*(1._dp - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2
     .     )*z + 2._dp*beta**2*z**2 + beta**3*z**3) - 2._dp*s*(1._dp 
     .     - 3._dp*beta**6 + beta*(2._dp + 2._dp*beta**2 - 5._dp*beta**4 
     .     - 4._dp*beta**6)*z + beta**2*(1._dp + 3._dp*beta**2 
     .     - 6._dp*beta**4)*z**2 + 2._dp*beta**5*(1._dp + beta**2)*z**3 
     .     + 4._dp*beta**6*z**4 + beta**5*z**5))*xI2(t,mt**2,mh**2,musq,
     .     ep))/(s*(1._dp + beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) - 
     .     (1._dp*(8._dp*mh**4 - 4._dp*mh**2*s*(2._dp - 3._dp*beta**2 
     .     - beta*z) + s**2*(6._dp - 10._dp*beta**2 + 5._dp*beta**4 
     .     - 2._dp*beta*(1._dp - 2._dp*beta**2)*z + beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/s + ((s**2*(6._dp 
     .     - 10._dp*beta**2 + 5._dp*beta**4 + beta*(3._dp - 4._dp*beta**2 
     .     + 4._dp*beta**4)*z - beta**2*(1._dp - 4._dp*beta**2)*z**2 
     .     + beta**3*z**3) + (1._dp + beta*z)*(8._dp*mh**4 - 4._dp*mh**2*s*
     .     (2._dp - 3._dp*beta**2 - beta*z)))*xI3(0._dp,mt**2,t,mt**2
     .     ,mt**2,mh**2,musq,ep))/s + ((2._dp*mh**4*z*(5._dp 
     .     - 8._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     - 2._dp*beta*mh**2*s*(1._dp + beta**2 - 2._dp*beta*(5._dp 
     .     - 6._dp*beta**2)*z - (1._dp - beta**2)*z**2 
     .     + 3._dp*beta*(1._dp - beta**2)*z**3) + beta**2*s**2*
     .     (beta - 2._dp*beta**3 - (7._dp - 13._dp*beta**2 
     .     + 8._dp*beta**4)*z + beta*(1._dp - 2._dp*beta**2)*z**2 
     .     - 2._dp*beta**2*(1._dp - beta**2)*z**3))*
     .     xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2,musq,ep))/(beta*s) 
     .     + (0.5_dp*(-16._dp*mh**6 + 16._dp*mh**4*s*(1._dp - 2._dp*beta**2 
     .     - beta*z) - 2._dp*mh**2*s**2*(6._dp - 13._dp*beta**2 
     .     + 10._dp*beta**4 - 6._dp*beta*(1._dp - 2._dp*beta**2)*z 
     .     + beta**2*(2._dp + beta**2)*z**2) - s**3*(2._dp + beta**2 
     .     - 6._dp*beta**4 + 4._dp*beta**6 + beta*(5._dp - 10._dp*beta**2 
     .     + 8._dp*beta**4)*z + beta**4*(1._dp + 2._dp*beta**2)*z**2 
     .     + beta**3*z**3))*xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,
     .     mt**2,mh**2,musq,ep))/s))/(8*pi**2)
     
     
     
c--- use the value of "tevscale" (passed via anomcoup.f)
c--- as an anomalous top Yukawa coupling:
c--- g(top Y) = (tevscale) x g(SM, top Y)
c--- it therefore affects Higgs diagrams as the square
      vrts(5)=vrts(5)*tevscale**2
      slf(5) =slf(5) *tevscale**2
      vrt(5) =vrt(5) *tevscale**2
      bx(5)  =bx(5)  *tevscale**2
    
!     MARKUS: marking contributions that do not scale like gvt,gat,gw
      vrts(3:5) = vrts(3:5) * g_rest 
      slf(3:5)  = slf(3:5)  * g_rest 
      vrt(3:5)  = vrt(3:5)  * g_rest 
      bx(3:5)   = bx(3:5)   * g_rest 
      trih      = trih      * g_rest
    
    
!          print *, "vrts(1)",vrts(1)/(4.7935059840374244d-003 )
!          print *, "vrts(2)",vrts(2)/( 4.5112451732988584d-003)
!          print *, "slf(1)",slf(1)/(9.4776826372363968d-003 )
!          print *, "slf(2)",slf(2)/(-2.7968401948508790d-003)
!          print *, "vrt(1)",vrt(1)/(-3.3877877304341691d-002)
!          print *, "vrt(2)",vrt(2)/(-2.4077932593682856d-003)
!          print *, "bx(1)",bx(1)/(-1.1560399205759183d-002 )
!          print *, "bx(2)",bx(2)/(1.1388151472327629d-002)
!          pause
    
      vew = 
     .     + (vrts(1) + vrts(2) + vrts(3) + vrts(4) + vrts(5))/2._dp
     .     + slf(1) + slf(2) + slf(3) + slf(4) + slf(5)
     .     + vrt(1) + vrt(2) + vrt(3) + vrt(4) + vrt(5)
     .     + bx(1) + bx(2) + bx(3) + bx(4) + bx(5)

      vew = vew + (trih + trizx)/2._dp

      vew = vew/born

!          print *, vrts(1) , vrts(2) , vrts(3) , vrts(4) , vrts(5)
!          print *, slf(1) , slf(2) , slf(3) , slf(4) , slf(5)
!          print *, vrt(1) , vrt(2) , vrt(3) , vrt(4) , vrt(5)
!          print *, bx(1) , bx(2) , bx(3) , bx(4) , bx(5)
!          print *, trih
!          print *, trizx
!          print *, vew
!          pause
      
      end subroutine ggQQb_ew_oneloop

