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
      vrts(3) = 
     .     (-8._dp*alpha*chifacs*genfacschi*sigma0*
     .     (beta**2*s*db0(mt**2,mt**2,mz**2) + 2._dp*
     .     (-xI2(mt**2,mt**2,mz**2,musq,ep) 
     .     + xI2(s,mt**2,mt**2,musq,ep)) 
     .     + (2._dp*mz**2 + beta**2*s)*
     .     xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep)))/pi


!            print *, "MCFM vrts(3)",vrts(3)
!            
!            print *, "MAKU vrts   ",
!      -        1d0/48d0*alpha*chifacs*sigma0/pi*3*       ! empirical prefactor to match MCFM
!      -       (24*(-1 + beta**2)**2*z**2*
!      -    (beta**2*s*db0(mt**2,mt**2,mz**2) - 
!      -      2*xI2(mt**2,mt**2,mz**2,musq,ep) + 2*xI2(s,mt**2,mt**2,musq,ep)+ 
!      -      (2*MZ**2 + beta**2*s)*
!      -    xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep)))/
!      -  (-1 + beta**2*z**2)
!            pause
           

     
!     MARKUS: phi contribution vertex, s-channel gluon 
!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
      vrts(4) = 
     .     (-8._dp*alpha*genfacs*(0.125_dp/sw2)*sigma0*
     .     (-0.03125_dp*(16._dp*beta**2*mb**4*s - 4._dp*beta**2*(1._dp 
     .     - beta**2)*mw**2*s**2 + beta**2*(1._dp - beta**2 )**2*s**3
     .     - 8._dp*mb**2*(2._dp*beta**2*mw**2*s + beta**2*(1._dp 
     .     - beta**2)*s**2))*db0(mt**2,mb**2,mw**2) + 0.0625_dp*
     .     (16._dp*mb**4 - 4._dp*(1._dp - beta**2)*mw**2*s 
     .     - (1._dp - beta**4)*s**2 + mb**2*(-16._dp*mw**2 
     .     + 8._dp*beta**2*s))*(xI2(mt**2,mb**2,mw**2,musq,ep) 
     .     - xI2(s,mb**2,mb**2,musq,ep)) + 0.015625_dp*(64._dp*mb**6 
     .     + 16._dp*(1._dp - beta**2)*mw**4*s + 8._dp*(1._dp 
     .     - beta**4)*mw**2*s**2 + (1._dp - beta**2)**3*s**3 
     .     - 16._dp*mb**4*(8._dp*mw**2 + (1._dp - beta**2)*s) 
     .     + 4._dp*mb**2*(16._dp*mw**4 - (1._dp - beta**2)**2*
     .     s**2))*xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)))
     .     /(mw**2*pi)


!            print *, "MCFM vrts(4)",vrts(4)
!            
!            print *, "MAKU vrts   ",
!      -     1d0/32d0*alpha*(0.125_dp/sw2)*sigma0/(mw**2*pi)*
!      -      (3*(-1 + beta**2)**2*z**2*
!      -    (2*beta**2*s*(4*MW**2 + (-1 + beta**2)*s)*
!      -       db0(mt**2,mb**2,mw**2) - 
!      -      4*(4*MW**2 + s + beta**2*s)*
!      -       (xI2(mt**2,0d0,mw**2,musq,ep) - xI2(s,0d0,0d0,musq,ep))+ 
!      -      (4*MW**2 + (-1 + beta)**2*s)*
!      -       (4*MW**2 + (1 + beta)**2*s)*
!      -       xI3(s,mt**2,mt**2,0d0,0d0,mw**2,musq,ep)))/
!      -  (-1 + beta**2*z**2)
!            pause
     
     
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



      slf(3) = 
     .     (2._dp*alpha*chifac*genfacself*sigma0*
     .     ((0.125_dp*(1._dp - beta**2)*mz**2*s*(1._dp 
     .     - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mt**2,
     .     mz**2))/(1._dp + beta*z) + (0.5_dp*(1._dp - beta**4*z**4 
     .     + (1._dp - beta**2)*(beta**2 - 3._dp*beta**2*z**2))*
     .     (-xI1(mt**2,musq,ep) + xI1(mz**2,musq,ep)))/(1._dp 
     .     + beta**2 + 2._dp*beta*z) - (0.25_dp*mz**2*(2._dp + 2._dp*beta**2 
     .     - 5._dp*beta**4 + 2._dp*beta**6 + beta**3*(3._dp 
     .     - 2._dp*beta**2)*z - 3._dp*beta**2*(2._dp - 3._dp*beta**2 
     .     + beta**4)*z**2 - 3._dp*beta**3*(1._dp - beta**2)*z**3 
     .     - beta**4*(2._dp - beta**2)*z**4 
     .     - beta**5*z**5)*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     (1._dp + beta*z)**2 + (0.125_dp*(1._dp - beta**2)*
     .     (s*(1._dp + beta**2 - beta**4 - 3._dp*beta**2*(1._dp 
     .     - beta**2)*z**2 - beta**4*z**4) + 
     .     (2._dp*beta**2*mz**2*(1._dp - z**2)*(2._dp + beta**2 
     .     - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z)))/(1._dp + beta*z)**2)*
     .     xI2(t,mt**2,mz**2,musq,ep))/
     .     (1._dp + beta**2 + 2._dp*beta*z)))/pi

!          print *, "MCFM slf(3)", slf(3)
!          print *, "MAKU slf(3)",2._dp*alpha*chifac*sigma0/pi/64d0*(
!      - (-8*(-1 + beta)*(1 + beta)*MZ**2*s*(-7 + 9*beta*z)*
!      -     (-1 + beta*(z + beta*(-2 - beta*z + 3*z**2 + 
!      -             beta**2*(2 - 3*z**2 + z**4))))*db0(MT**2,MT**2,MZ**2))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)**2) - 
!      -  (32*(-7 + 9*beta*z)*(-1 + beta**2*(-1 + 3*z**2) + 
!      -       beta**4*(1 - 3*z**2 + z**4))*xI1(mt**2,musq,ep))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
!      -  (32*(-7 + 9*beta*z)*(-1 + beta**2*(-1 + 3*z**2) + 
!      -       beta**4*(1 - 3*z**2 + z**4))*xI1(mz**2,musq,ep))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
!      -  (16*MZ**2*(-7 + 9*beta*z)*(2 + 
!      -       beta**2*(2 - 6*z**2 - 3*beta*z*(-1 + z**2) + 
!      - beta**2*(-5 + 9*z**2 - 2*z**4) + beta**4*(2 - 3*z**2 + z**4) - 
!      - beta**3*z*(2 - 3*z**2 + z**4)))*xI2(mt**2,mt**2,mz**2,musq,ep))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)**3) + 
!      -  (8*(1 - beta**2)*(-7 + 9*beta*z)*
!      - (-s - 2*beta*s*z + beta**3*(-2*(3*MZ**2 + s)*z + 6*(MZ**2 + s)*z**3)+
!      -  2*beta**5*z*(2*MZ**2 + s - 3*(MZ**2 + s)*z**2 + (MZ**2 + s)*z**4) + 
!      -  beta**2*(4*MZ**2*(-1 + z**2) + s*(-1 + 2*z**2)) + 
!      -  beta**4*(s + 2*(-1 + z**2)*(MZ**2 + 2*s*z**2)) + 
!      - beta**6*(s*z**2*(1 - 3*z**2 + z**4) + 2*MZ**2*(2 - 3*z**2 + z**4)))*
!      -     xI2(t,mt**2,mz**2,musq,ep))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z)) )   
!            pause


!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
      slf(4) = 
     .     (alpha*genfacself*(0.125_dp/sw2)*sigma0*
     .     ((-0.03125_dp*(16._dp*mb**2*(mb**2 - mw**2) 
     .     - 4._dp*(1._dp - beta**2)*(2._dp*mb**2 + mw**2)*s 
     .     + (1._dp - beta**2)**2*s**2)*(1._dp - beta**4*z**4 
     .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mb**2,mw**2))/(1._dp + beta*z) 
     .     + (0.5_dp*(4._dp*mb**2 + (1._dp - beta**2)*s)*(1._dp 
     .     - beta**4*z**4 + beta**2*(1._dp - beta**2)*(1._dp 
     .     - 3._dp*z**2))*(-xI1(mb**2,musq,ep) 
     .     + xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp + beta**2 + 2._dp*beta*z)) - (0.0625_dp*(8._dp*(1._dp 
     .     - beta**2)**2*mb**2*s*(1._dp + beta**2 - 2._dp*beta**4 
     .     - beta**2*(2._dp - 3._dp*beta**2)*z**2 
     .     - beta**4*z**4) - 4._dp*(4._dp*mb**2*(-mb**2 + mw**2) 
     .     + (1._dp - beta**2)*mw**2*s)*(-2._dp - 2._dp*beta**2 
     .     + 5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp 
     .     - 2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 
     .     + beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5) 
     .     + beta**2*(1._dp - beta**2)**2*s**2*(1._dp - z**2)*
     .     (2._dp + beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     ((1._dp - beta**2)*s*(1._dp + beta*z)**2) + 
     .     (0.0625_dp*(-8._dp*mb**2*s*(-2._dp - 2._dp*beta**2 + 4._dp*beta**4 
     .     + beta**6 - 2._dp*beta**8 - 2._dp*beta*(2._dp + beta**2 
     .     - 4._dp*beta**4 + 2._dp*beta**6)*z + beta**2*(4._dp 
     .     - 7._dp*beta**2 - beta**4 + 3._dp*beta**6)*z**2 
     .     + 2._dp*beta**3*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**3 
     .     + beta**4*(5._dp - 3._dp*beta**2 - beta**4)*z**4 
     .     + 2._dp*beta**5*(2._dp - beta**2)*z**5 + beta**6*z**6) 
     .     - beta**2*(1._dp - z**2)*(2._dp + beta**2 
     .     - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z))*(16._dp*mb**2*(mb**2 - mw**2) 
     .     - 4._dp*(1._dp - beta**2)*mw**2*s - s**2*(1._dp 
     .     - beta**4 + 2._dp*beta*(1._dp - beta**2)*z)))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta*z)**2*
     .     (1._dp + beta**2 + 2._dp*beta*z))))/(mw**2*pi)


!          print *, "MCFM slf(4)", slf(4)
!          print *, "MAKU slf(4)", alpha*(0.125_dp/sw2)/(mw**2*pi)*sigma0/128d0*
!      -    (  (-4*(-1 + beta)*(1 + beta)*s*(4*MW**2 + (-1 + beta**2)*s)*
!      -     (-7 + 9*beta*z)*(-1 + 
!      -       beta*(z + beta*(-2 - beta*z + 3*z**2 + 
!      -             beta**2*(2 - 3*z**2 + z**4))))*db0(MT**2,mb**2,MW**2))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)**2) + 
!      -  (64*(-7 + 9*beta*z)*(-1 + beta**2*(-1 + 3*z**2) + 
!      -       beta**4*(1 - 3*z**2 + z**4))*xI1(mw**2,musq,ep))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
!      -  (8*(-7 + 9*beta*z)*(4*MW**2*
!      -        (2 + beta**2*(2 - 6*z**2 - 3*beta*z*(-1 + z**2) + 
!      -             beta**2*(-5 + 9*z**2 - 2*z**4) + 
!      -             beta**4*(2 - 3*z**2 + z**4) - beta**3*z*(2 - 3*z**2+z**4)
!      -             )) + beta**2*(-1 + beta**2)*s*(-1 + z**2)*
!      -        (2 + beta*(3*z + beta*(1 + beta*(beta + z)*(-2 + z**2)))))*
!      -   xI2(mt**2,0d0,mw**2,musq,ep))/(3.*(-1 + beta*z)*(1 + beta*z)**3) - 
!      -  (8*beta**2*(-1 + beta**2)*(-7 + 9*beta*z)*(-1 + z**2)*
!      -     (4*MW**2 + s*(1 + beta**2 + 2*beta*z))*
!      -     (2 + beta*(3*z + beta*(1 + beta*(beta + z)*(-2 + z**2))))*
!      -    xI2(t,0d0,mw**2,musq,ep))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)**3*(1 + beta**2 + 2*beta*z)) )
!           pause
         
         
         
         

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




      vrt(3) = 
     .     (2._dp*alpha*chifac*genfacvert*sigma0*
     .     (0.125_dp*(1._dp - beta**2)**2*s - (0.25_dp*(1._dp 
     .     - beta**2)*mz**2*s*(1._dp - beta**4*z**4 
     .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mt**2,mz**2))/
     .     (1._dp - beta*z) + (0.5_dp*(1._dp + 2._dp*beta**2 
     .     - beta**4 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) - 2._dp*(1._dp 
     .     - beta**2)*(2._dp*beta**2*z**2 + 3._dp*beta**3*z**3))*
     .     (xI1(mt**2,musq,ep) - xI1(mz**2,musq,ep)))/
     .     ((1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + (0.125_dp*(-(1._dp - beta**2)**2*s*(1._dp 
     .     - beta*z)*(1._dp + beta**2 + 2._dp*beta*z) 
     .     + 2._dp*mz**2*(1._dp + 8._dp*beta**2 - 11._dp*beta**4 
     .     + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**2*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**4*
     .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mt**2,mz**2,musq,ep))/(1._dp - beta**2*z**2) 
     .     + (0.125_dp*(1._dp - beta**2)*(2._dp*mz**2*(1._dp 
     .     - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     - 8._dp*beta**2 + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp 
     .     + 2._dp*beta**2 - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
     .     + beta*s*(-beta*(1._dp + beta**4) + (1._dp - 7._dp*beta**2 
     .     + beta**4 + beta**6)*z + beta*(3._dp - 12._dp*beta**2 
     .     + 7._dp*beta**4)*z**2 + 6._dp*beta**2*(1._dp - beta**2)*
     .     z**3 + 2._dp*beta**3*(4._dp - 3._dp*beta**2)*z**4 
     .     + 4._dp*beta**4*z**5 + 2._dp*beta**5*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) 
     .     - 0.0625_dp*(1._dp - beta**2)**2*s**2*
     .     (1._dp + beta**2 + 2._dp*beta*z)*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep)))/pi


!          print *, "MCFM vrt(3)",vrt(3)
!          print *, "MAKU vrt(3)",2._dp*alpha*chifac*sigma0/48d0/Pi*
!      -     ((-7 + 9*beta*z)*((2*(-1 + beta**2)*s*(1 + beta*z)*
!      -         ((-1 + beta**2)*(-1 + beta*z) + 
!      -           2*MZ**2*(-1 + beta*
!      -   (z + beta*(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4))))
!      -             *db0(mt**2,mt**2,mz**2)))/(-1 + beta*z) + 
!      -      (8*(1 + beta*z)*(-1 + beta*
!      -            (-z + beta*(-2 + 4*z**2 + beta*(-4*z + 6*z**3) + 
!      -                 beta**3*z*(3 - 6*z**2 + 2*z**4) + 
!      -                 beta**2*(1 - 4*z**2 + 2*z**4))))*xI1(mt**2,musq,ep))/
!      -       ((-1 + beta*z)*(1 + beta**2 + 2*beta*z)) - 
!      -      (8*(1 + beta*z)*(-1 + beta*
!      -            (-z + beta*(-2 + 4*z**2 + beta*(-4*z + 6*z**3) + 
!      -                 beta**3*z*(3 - 6*z**2 + 2*z**4) + 
!      -                 beta**2*(1 - 4*z**2 + 2*z**4))))*xI1(mz**2,musq,ep))/
!      -       ((-1 + beta*z)*(1 + beta**2 + 2*beta*z)) - 
!      -    (2*((-1 + beta**2)**2*s*(-1 + beta*z)*(1 + beta**2 + 2*beta*z) +
!      -           2*MZ**2*(1 + beta*
!      -               (z + beta*(8 - 10*z**2 + beta*(4*z - 6*z**3) + 
!      -                    beta**2*(-11 + 16*z**2 - 4*z**4) + 
!      -                    beta**3*z*(-3 + 6*z**2 - 2*z**4) + 
!      -  2*beta**4*(2 - 3*z**2 + z**4)))))*xI2(mt**2,mt**2,mz**2,musq,ep))/
!      -       (-1 + beta*z) + (2*(-1 + beta**2)*
!      -         (beta*s*(z + beta*(-1 + beta**5*z + 3*z**2 + 
!      -        4*beta**2*z**2*(-3 + 2*z**2) + beta*z*(-7 + 6*z**2) + 
!      -                 beta**3*(z - 6*z**3 + 4*z**5) + 
!      -                 beta**4*(-1 + 7*z**2 - 6*z**4 + 2*z**6))) + 
!      -           2*MZ**2*(1 + beta*
!      -               (z + beta*(-4 + 2*z**2 + beta**2*(-3 + 4*z**2) + 
!      -     beta*(-8*z + 6*z**3) + 2*beta**4*(2 - 3*z**2 + z**4) + 
!      -    beta**3*z*(5 - 6*z**2 + 2*z**4)))))*
!      -      xI2(t,mt**2,mz**2,musq,ep))/
!      -       ((-1 + beta*z)*(1 + beta**2 + 2*beta*z)) - 
!      -      (-1 + beta**2)**2*s**2*(1 + beta*z)*(1 + beta**2 + 2*beta*z)*
!      -    xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep)))/
!      -  (1 + beta*z)**3
!            pause




!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
      vrt(4) = 
     .     (alpha*genfacvert*(0.125_dp/sw2)*sigma0*
     .     (0.125_dp*(1._dp - beta**2)*(4._dp*mb**2 + (1._dp 
     .     - beta**2)*s) + (0.0625_dp*(16._dp*mb**2*(mb**2 
     .     - mw**2) - 4._dp*(1._dp - beta**2)*(2._dp*mb**2 
     .     + mw**2)*s + (1._dp - beta**2)**2*s**2)*(1._dp 
     .     - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mb**2,
     .     mw**2))/(1._dp - beta*z) + (0.5_dp*(4._dp*mb**2 + 
     .     (1._dp - beta**2)*s)*(1._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**4
     .     *z**4*(1._dp + beta*z) - 2._dp*(1._dp - beta**2)*
     .     (2._dp*beta**2*z**2 + 3._dp*beta**3*z**3))*(xI1(mb**2,musq,ep) 
     .     - xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.0625_dp*(-8._dp*(1._dp - beta**2)**2*mb**2*s*(1._dp 
     .     - beta**2 + 4._dp*beta**4 + beta*(1._dp - beta**2)*z 
     .     - 6._dp*beta**4*z**2 + 2._dp*beta**4*z**4) + (-16._dp*mb**4 
     .     + 4._dp*mw**2*(4._dp*mb**2 + (1._dp - beta**2)*s))*
     .     (1._dp + 8._dp*beta**2 - 11._dp*beta**4 + 4._dp*beta**6 + 
     .     beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**2*
     .     (5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 - 6._dp*beta**3*
     .     (1._dp - beta**2)*z**3 - 2._dp*beta**4*
     .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5) 
     .     - (1._dp - beta**2)**2*s**2*(1._dp - 4._dp*beta**2 
     .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
     .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 
     .     + 2._dp*beta**5*z**5))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     ((1._dp - beta**2)*s*(1._dp - beta**2*z**2)) 
     .     + (0.0625_dp*((-16._dp*mb**4 + 4._dp*mw**2*(4._dp*mb**2 
     .     + (1._dp - beta**2)*s))*(1._dp - 4._dp*beta**2 
     .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
     .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
     .     + (1._dp - beta**2)*s**2*(1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 
     .     + beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     + 2._dp*beta**2*(1._dp + 2._dp*beta**2 - 3._dp*beta**4)*z**2 
     .     + 6._dp*beta**3*(1._dp - beta**2)*z**3 + 2._dp*beta**6*z**4 
     .     + 2._dp*beta**5*z**5) + 8._dp*beta*mb**2*s*(-beta*(3._dp 
     .     - 4._dp*beta**2 - beta**4 + 4._dp*beta**6) + (1._dp 
     .     - 11._dp*beta**2 + 13._dp*beta**4 - 7._dp*beta**6)*z 
     .     + beta*(5._dp - 18._dp*beta**2 + 5._dp*beta**4 
     .     + 6._dp*beta**6)*z**2 + 2._dp*beta**2*(5._dp - 11._dp*beta**2 
     .     + 6._dp*beta**4)*z**3 + 2._dp*beta**3*(5._dp - 3._dp*beta**2 
     .     - beta**4)*z**4 + 4._dp*beta**4*(2._dp 
     .     - beta**2)*z**5 + 2._dp*beta**5*z**6))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(1._dp - beta**2*z**2)) 
     .     - 0.25_dp*(1._dp - beta**2)*mb**2*(-4._dp*mb**2 + 3._dp*s 
     .     + beta**2*s + 4._dp*beta*s*z)*
     .     xI3(0._dp,mt**2,t,mb**2,mb**2,mw**2,musq,ep)))/(mw**2*pi)

           
!            print *, "MCFM vrt(4)", vrt(4)
!            
!            print *, "MAKU vrt(4)",alpha*(0.125_dp/sw2)*sigma0/(mw**2*pi)/48d0*
!      -     ((-7 + 9*beta*z)*(((1 + beta*z)*
!      -         (-1 + beta*(z + beta*
!      -        (-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4))))*
!      -  ((-1 + beta**2)*s*(4*MW**2 + (-1+beta**2)*s)*db0(mt**2,mb**2,mw**2)
!      -  -8*xI1(mw**2,musq,ep) + 2*(4*MW**2 + s - beta**2*s)*xI2(mt**2,0d0,mw**2,musq,ep)))/
!      -       (-1 + beta*z) + (-1 + beta**2)*
!      -       (2*(-1 + beta**2)*s*(1 + beta*z) + 
!      - (8*beta**2*(-1 + z)*(1 + z)*(1 + beta*z)*(-1 + beta**2*(-2 + z**2))*
!      -   xI1(mw**2,musq,ep) - (1 + beta**2 + 2*beta*z)*
!      -             (4*MW**2*(1 + beta*
!      -                   (-z + beta*
!      -      (-3 + beta*z + 2*z**2 + 2*beta**2*(2 - 3*z**2 + z**4)))) +
!      -               (-1 + beta**2)*s*
!      -                (-3 + beta*(-z + 
!      -   beta*(-3 + beta*z + 6*z**2 + 2*beta**2*(2 - 3*z**2 + z**4)))))*
!      -    xI2(mt**2,0d0,mw**2,musq,ep) + 
!      -            (4*MW**2 + s*(1 + beta**2 + 2*beta*z))*
!      -             (1 + beta*(z + beta*
!      -  (-4 + 2*z**2 + beta**2*(-3 + 4*z**2) + beta*(-8*z + 6*z**3) + 
!      -  2*beta**4*(2 - 3*z**2 + z**4) + beta**3*z*(5 - 6*z**2 + 2*z**4)
!      -  )))*xI2(t,0d0,mw**2,musq,ep))/
!      -  ((-1 + beta*z)*(1 + beta**2 + 2*beta*z)))))/(1 + beta*z)**3           
!               pause

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


      bx(3) = 
     .     (0.5_dp*alpha*chifac*genfacbox*mt**2*sigma0*
     .     (-beta*z*(1._dp - z**2) + (4._dp*beta*(1._dp 
     .     - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
     .     (xI1(mt**2,musq,ep) - xI1(mz**2,musq,ep)))/(s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + ((-2._dp*beta*(1._dp 
     .     - beta**2)*s*(2._dp - beta**2 + beta*z - z**2 
     .     - beta*z**3) + 2._dp*mz**2*(6._dp*beta**3 + (5._dp 
     .     - 8._dp*beta**2 + 4._dp*beta**4)*z + beta*(5._dp 
     .     - 12._dp*beta**2)*z**2 - (3._dp - 4._dp*beta**2 
     .     + 2._dp*beta**4)*z**3 - beta*(3._dp - 4._dp*beta**2)*z**4)
     .     )*xI2(mt**2,mt**2,mz**2,musq,ep))/(beta*s*(1._dp + beta*z)) 
     .     - (1._dp*(-beta*s*(2._dp + beta*z)*(1._dp - z**2) 
     .     + 2._dp*mz**2*z*(5._dp - 4._dp*beta**2 - (3._dp 
     .     - 2._dp*beta**2)*z**2))*xI2(s,mt**2,mt**2,musq,ep))/(beta*s)
     .     + ((4._dp*beta*mz**2*(beta + z)*(1._dp - 3._dp*beta**2 
     .     + beta*(1._dp - 2._dp*beta**2)*z + 2._dp*beta**2*z**2 
     .     + beta**3*z**3) + 2._dp*s*(1._dp - 4._dp*beta**2 + beta**6 
     .     + beta*(2._dp - 10._dp*beta**2 + 3._dp*beta**4)*z + beta**2*
     .     (3._dp - 5._dp*beta**2)*z**2 + 2._dp*beta**3*(3._dp 
     .     - beta**2)*z**3 + 4._dp*beta**4*z**4 + beta**5*z**5))*
     .     xI2(t,mt**2,mz**2,musq,ep))/(s*(1._dp + beta*z)*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + ((8._dp*mz**4 + 4._dp*beta*mz**2*s*
     .     (beta + z) + s**2*(2._dp - 2._dp*beta**2 + beta**4 
     .     + 2._dp*beta*z + beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/s
     .     - (1._dp*(s**2*(2._dp - 2._dp*beta**2 + beta**4 
     .     + 3._dp*beta*z + 3._dp*beta**2*z**2 + beta**3*z**3) + 4._dp*
     .     (1._dp + beta*z)*(2._dp*mz**4 + beta*mz**2*s*(beta + z)))*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep))/s - 
     .     (1._dp*(-beta**2*s**2*(beta + z)*(1._dp + beta*z) 
     .     + 2._dp*mz**4*z*(5._dp - 8._dp*beta**2 - (3._dp 
     .     - 2._dp*beta**2)*z**2) - 2._dp*beta*mz**2*s*(1._dp + beta**2 
     .     - 2._dp*beta*(1._dp - 2._dp*beta**2)*z - (1._dp 
     .     - beta**2)*z**2 + beta*(1._dp - beta**2)*z**3))*
     .     xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep))/(beta*s) + 
     .     (0.5_dp*(16._dp*mz**6 + 16._dp*beta*mz**4*s*(beta + z) 
     .     + beta*s**3*(beta + z)*(1._dp + beta*z)**2 + 2._dp*mz**2*s**2*
     .     (2._dp - beta**2 + 2._dp*beta**4 + 2._dp*beta*(1._dp 
     .     + 2._dp*beta**2)*z + beta**2*(2._dp + beta**2)*z**2))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,mt**2,mz**2,musq,
     .     ep))/s))/pi


!         print *, "MCFM bx(3) ", bx(3) 
!         
! !       this is the result from   TheRedAmpList[[9]] of DiagGG1
!         print *, "MAKU bx(3) ", 1d0/mt**2/Pi* ! empirical prefactor to match MCFM
!      -   0.5_dp*alpha*chifac*mt**2*sigma0*   ! prefactor from MCFM without the genfacbox
!      -   ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -    ((-2*beta*s*z*(-1 + z**2))/(-1 + beta**2*z**2) - 
!      -   (8*beta*(-1 + z**2)*(z + 2*beta*(-2 + z**2))*xI1(mt**2,musq,ep))/
!      -    ((-1 + beta*z)*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
!      -   (8*beta*(-1 + z**2)*(z + 2*beta*(-2 + z**2))*xI1(mz**2,musq,ep))/
!      -    ((-1 + beta*z)*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
!      -      (4*(beta**5*s + MZ**2*z*(-5 + 3*z**2) + 
!      -           beta**2*z*(8*MZ**2 + s - (4*MZ**2 + s)*z**2) + 
!      -           beta*(2*s - (5*MZ**2 + s)*z**2 + 3*MZ**2*z**4) + 
!      -           beta**4*z*(2*MZ**2*(-2 + z**2) + s*(-1 + z**2)) + 
!      -           beta**3*(s*(-3 + z**2) - 
!      -              2*MZ**2*(3 - 6*z**2 + 2*z**4)))*
!      -         xI2(mt**2,mt**2,mz**2,musq,ep))/
!      -       (beta*(-1 + beta*z)*(1 + beta*z)**2) + 
!      -      (2*(2*MZ**2*z*(5 - 3*z**2) + 2*beta*s*(-1 + z**2) + 
!      -           beta**2*z*(4*MZ**2*(-2 + z**2) + s*(-1 + z**2)))*
!      -         xI2(s,mt**2,mt**2,musq,ep))/(beta*(-1 + beta**2*z**2)) - 
!      -      (4*(s + beta**6*s + 2*beta*(MZ**2 + s)*z + 
!      -           2*beta**3*z*(2*MZ**2*(-1 + z**2) + s*(-5 + 3*z**2)) + 
!      -           beta**2*(2*MZ**2*(1 + z**2) + s*(-4 + 3*z**2)) + 
!      -           beta**5*z*(3*s + (-2 + z**2)*(2*MZ**2 + s*z**2)) + 
!      -           beta**4*(s*z**2*(-5 + 4*z**2) + 2*MZ**2*(-3 + z**4)))*
!      -         xI2(t,mt**2,mz**2,musq,ep))/
!      -       ((-1 + beta*z)*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z))
!      -       - (2*(8*MZ**4 + 4*beta*MZ**2*s*(beta + z) + 
!      -           s**2*(2 + beta**4 + 2*beta*z + beta**2*(-2 + z**2)))*
!      -  xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/(-1 + beta**2*z**2)+
!      -      (2*(8*MZ**4*(1 + beta*z) + 
!      -           4*beta*MZ**2*s*(beta + z)*(1 + beta*z) + 
!      -           s**2*(2 + beta*
!      -               (3*z + beta*(-2 + 3*z**2 + beta*(beta + z**3)))))*
!      - xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep))/(-1 + beta**2*z**2) + 
!      -      (2*(2*MZ**4*z*(5 - 3*z**2) + 2*beta*MZ**2*s*(-1 + z**2) - 
!      -           beta**3*s*(2*MZ**2 + s)*(1 + z**2) - 
!      -           beta**4*s*z*(s - 2*MZ**2*(-4 + z**2)) + 
!      -           beta**2*z*(-s**2 + 4*MZ**4*(-4 + z**2) - 
!      -              2*MZ**2*s*(-2 + z**2)))*
!      -         xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep))/
!      -       (beta*(-1 + beta**2*z**2)) - 
!      -      ((16*MZ**6 + 16*beta*MZ**4*s*(beta + z) + 
!      -           beta*s**3*(beta + z)*(1 + beta*z)**2 + 
!      -           2*MZ**2*s**2*
!      -            (2 + beta*(2*z + 
!      -                 beta*(-1 + 4*beta*z + 2*z**2 + 
!      -                    beta**2*(2 + z**2)))))*
!      -         xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,mt**2,mz**2,musq,
!      .     ep))/(-1 + beta**2*z**2)))/24.
!              pause

             
             

!      MARKUS: replaced gw_sq coupl. by explicit expression  0.125_dp/sw2
      bx(4) = 
     .     (0.25_dp*alpha*genfacbox*(0.125_dp/sw2)*sigma0*
     .     (-0.25_dp*beta*(4._dp*mb**2 + (1._dp - beta**2)*s)*z*(1._dp 
     .     - z**2)*1 + (beta*(4._dp*mb**2 + (1._dp - beta**2)*s)*
     .     (1._dp - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
     .     (xI1(mb**2,musq,ep)*1 - xI1(mw**2,musq,ep)*1))/(s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) - (0.125_dp*((1._dp - beta**2
     .     )**2*s**2*(8._dp*beta - 2._dp*beta**3 - (5._dp 
     .     - 4._dp*beta**2 - 4._dp*beta**4)*z - 9._dp*beta*z**2 + (3._dp 
     .     - 4._dp*beta**2 - 2._dp*beta**4)*z**3 + 3._dp*beta*z**4) 
     .     + 16._dp*beta*(1._dp - beta**2)*mb**2*s*(2._dp 
     .     - 3._dp*beta**2 + beta*(3._dp - 2._dp*beta**2)*z 
     .     - (1._dp - 3._dp*beta**2)*z**2 - beta*(2._dp 
     .     - beta**2)*z**3 - beta**2*z**4) + (16._dp*mb**2*
     .     (mb**2 - mw**2) - 4._dp*(1._dp - beta**2)*mw**2*s)*
     .     (6._dp*beta**3 + (5._dp - 8._dp*beta**2 + 4._dp*beta**4)*z 
     .     + beta*(5._dp - 12._dp*beta**2)*z**2 - (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 - beta*(3._dp 
     .     - 4._dp*beta**2)*z**4))*xI2(mt**2,mb**2,mw**2,musq,ep)*1)/
     .     (beta*s*(1._dp + beta*z)) + (0.125_dp*((16._dp*mb**2*(mb**2 
     .     - mw**2) - 4._dp*(1._dp - beta**2)*mw**2*s)*z*(5._dp 
     .     - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     + 8._dp*beta*mb**2*s*(2._dp + beta*(5._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 - beta*(3._dp - 2._dp*beta**2)*z**3) + (1._dp 
     .     - beta**2)*s**2*(4._dp*beta - (5._dp - 3._dp*beta**2 
     .     - 4._dp*beta**4)*z - 4._dp*beta*z**2 + (3._dp - 3._dp*beta**2 
     .     - 2._dp*beta**4)*z**3))*xI2(s,mb**2,mb**2,musq,ep)*1)/(beta*s) 
     .     - (0.25_dp*(beta*(16._dp*mb**2*(mb**2 - mw**2) 
     .     - 4._dp*(1._dp - beta**2)*mw**2*s)*(beta + z)*(1._dp 
     .     - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2)*z 
     .     + 2._dp*beta**2*z**2 + beta**3*z**3) - (1._dp 
     .     - beta**2)*s**2*(1._dp + beta**2 + 2._dp*beta*z)*
     .     (2._dp - 5._dp*beta**2 + beta**4 + beta*(1._dp - 2._dp*beta**2 
     .     - 2._dp*beta**4)*z + beta**2*(3._dp - 2._dp*beta**2)*z**2 
     .     + beta**3*(2._dp + beta**2)*z**3 + beta**4*z**4) 
     .     - 8._dp*mb**2*s*(1._dp - 6._dp*beta**2 + 3._dp*beta**6 
     .     + beta*(2._dp - 16._dp*beta**2 + 7._dp*beta**4 + 2._dp*beta**6)*z 
     .     + beta**2*(4._dp - 9._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     + beta**3*(9._dp - 4._dp*beta**2 - beta**4)*z**3 
     .     + 2._dp*beta**4*(3._dp - beta**2)*z**4 + beta**5*z**5))*
     .     xI2(t,mb**2,mw**2,musq,ep)*1)/(s*(1._dp + beta*z)*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + (0.125_dp*(64._dp*mb**6 
     .     + 16._dp*mw**4*(4._dp*mb**2 + (1._dp - beta**2)*s) 
     .     + 8._dp*beta*(1._dp - beta**2)*mw**2*s**2*(beta + z) 
     .     + 32._dp*mb**2*mw**2*s*(1._dp + beta*z) + 16._dp*mb**4*
     .     (-8._dp*mw**2 + (1._dp - beta**2)*s - 2._dp*beta*s*z) 
     .     + 4._dp*mb**2*s**2*(3._dp - 2._dp*beta**2 + beta**4 
     .     + 4._dp*beta*(2._dp - beta**2)*z + 2._dp*beta**2*z**2) 
     .     + (1._dp - beta**2)*s**3*(3._dp - 2._dp*beta**2 + beta**4 
     .     + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep)*1)/s 
     .     - (0.125_dp*(64._dp*mb**6*(1._dp + beta*z) 
     .     + 16._dp*(1._dp - beta**2)*mw**4*s*(1._dp + beta*z) 
     .     + 8._dp*beta*(1._dp - beta**2)*mw**2*s**2*(beta + z)*
     .     (1._dp + beta*z) + (1._dp - beta**2)*s**3*(1._dp + beta*z)*
     .     (3._dp - 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp + beta**2)*z 
     .     + 2._dp*beta**2*z**2) - 16._dp*mb**4*(8._dp*mw**2*(1._dp 
     .     + beta*z) - s*(1._dp - beta**2 - beta*(3._dp 
     .     - beta**2)*z - 2._dp*beta**2*z**2)) + mb**2*(4._dp*s**2*
     .     (3._dp - 2._dp*beta**2 + beta**4 + 9._dp*beta*z - 2._dp*beta**3*z 
     .     - beta**5*z + 10._dp*beta**2*z**2 - 4._dp*beta**4*z**2 
     .     + 2._dp*beta**3*z**3) + 32._dp*(1._dp + beta*z)*(2._dp*mw**4 
     .     + mw**2*s*(1._dp + beta*z))))*xI3(0._dp,mt**2,t,mb**2,mb**2
     .     ,mw**2,musq,ep)*1)/s + (0.03125_dp*(-16._dp*(4._dp*mb**6 
     .     - 8._dp*mb**4*mw**2 + mw**4*(4._dp*mb**2 + (1._dp - beta**2
     .     )*s))*z*(5._dp - 8._dp*beta**2 - (3._dp - 2._dp*beta**2
     .     )*z**2) + 64._dp*beta*mb**2*mw**2*s*(1._dp + beta**2 
     .     + 2._dp*beta*z - (1._dp - beta**2)*z**2) 
     .     + 8._dp*(1._dp - beta**2)*mw**2*s**2*(2._dp*beta*(1._dp 
     .     + beta**2) - (5._dp - beta**2 - 8._dp*beta**4)*z 
     .     - 2._dp*beta*(1._dp - beta**2)*z**2 + (3._dp - beta**2 
     .     - 2._dp*beta**4)*z**3) - 16._dp*mb**4*s*(4._dp*beta*(1._dp 
     .     + beta**2) - (5._dp - 17._dp*beta**2 + 8._dp*beta**4)*z 
     .     - 4._dp*beta*(1._dp - beta**2)*z**2 + (3._dp - 9._dp*beta**2 
     .     + 2._dp*beta**4)*z**3) + 4._dp*mb**2*s**2*(8._dp*beta**3*(2._dp 
     .     - beta**2) + (5._dp - 2._dp*beta**2 + 21._dp*beta**4 
     .     - 8._dp*beta**6)*z + 8._dp*beta**3*(2._dp - beta**2)*z**2 
     .     - (3._dp - beta**4 - 2._dp*beta**6)*z**3) + (1._dp 
     .     - beta**2)*s**3*(4._dp*beta*(1._dp + beta**4) - (5._dp 
     .     - 22._dp*beta**2 + 9._dp*beta**4 - 8._dp*beta**6)*z 
     .     - 4._dp*beta*(1._dp - 2._dp*beta**2 - beta**4)*z**2 
     .     + (3._dp - 4._dp*beta**2 + 3._dp*beta**4 - 2._dp*beta**6)*z**3))*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)*1)/(beta*s) 
     .     - (0.03125_dp*(256._dp*mb**8 - 64._dp*(1._dp - beta**2)*
     .     mw**6*s - 16._dp*(1._dp - beta**2)*mw**4*s**2*(1._dp 
     .     + 3._dp*beta**2 + 4._dp*beta*z) - (1._dp - beta**2)*
     .     s**4*(1._dp + beta**2 + 2._dp*beta*z)*(3._dp - 2._dp*beta**2 
     .     + beta**4 + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2) 
     .     - 4._dp*(1._dp - beta**2)*mw**2*s**3*(3._dp + 3._dp*beta**4 
     .     + 4._dp*beta*(1._dp + 2._dp*beta**2)*z + 6._dp*beta**2*z**2) 
     .     - 128._dp*mb**6*(6._dp*mw**2 + beta*s*z*(2._dp + beta*z)) 
     .     + 32._dp*mb**4*(24._dp*mw**4 + 2._dp*mw**2*s*(3._dp + beta**2 
     .     + 8._dp*beta*z + 2._dp*beta**2*z**2) - s**2*(1._dp 
     .     - 3._dp*beta**2 + beta**4 - 6._dp*beta*z + 2._dp*beta**3*z 
     .     - 5._dp*beta**2*z**2 + 2._dp*beta**4*z**2)) - 8._dp*mb**2*
     .     (32._dp*mw**6 + 16._dp*mw**4*s*(1._dp + beta**2 + 2._dp*beta*z) 
     .     + 2._dp*mw**2*s**2*(5._dp + beta**4 + 12._dp*beta*z 
     .     + 2._dp*beta**2*(2._dp + beta**2)*z**2) + s**3*(-2._dp 
     .     + 6._dp*beta**2 - 2._dp*beta**4 + 4._dp*beta*z + 4._dp*beta**3*z 
     .     - 2._dp*beta**5*z + 9._dp*beta**2*z**2 - 4._dp*beta**4*z**2 
     .     + beta**6*z**2 + 2._dp*beta**3*z**3)))*xI4(0._dp,0._dp,mt**2,
     .     mt**2,s,t,mb**2,mb**2,mb**2,mw**2,musq,ep)*1)/s))/(mw**2*pi)

!         print *, "MCFM bx(4) ", bx(4) 
!         
! !       this is the result from   TheRedAmpList[[7]] of DiagGG2
!         print *, "MAKU bx(4) ",0.5d0/mt**2/Pi* ! empirical prefactor to match MCFM
!      -   0.5_dp*alpha*chifac*mt**2*sigma0*(   ! prefactor from MCFM without the genfacbox        
!      -  (beta*(-1 + beta**2)*s*z*(-7 + 9*beta*z)*(-1 + z**2))/
!      -   (6 - 6*beta**2*z**2) + 
!      -  (2*beta*(-1 + beta**2)*(-7 + 9*beta*z)*(-1 + z**2)*
!      -     (z + 2*beta*(-2 + z**2))*xI1(mw**2,musq,ep))/
!      -   (3.*(-1 + beta*z)*(1 + beta*z)*
!      -     (1 + beta**2 + 2*beta*z)) + 
!      -  ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -     (2*beta**5*s + 2*beta**6*s*z*(-2 + z**2) + 
!      -       (4*MW**2 + s)*z*(-5 + 3*z**2) + 
!      -       beta**2*z*(s*(9 - 7*z**2) - 
!      -          16*MW**2*(-2 + z**2)) + 
!      -       2*beta**4*(s*z**3 + 4*MW**2*z*(-2 + z**2)) + 
!      -       beta*(4*MW**2*z**2*(-5 + 3*z**2) + 
!      -          s*(8 - 9*z**2 + 3*z**4)) - 
!      -       beta**3*(8*MW**2*(3 - 6*z**2 + 2*z**4) + 
!      -          s*(10 - 9*z**2 + 3*z**4)))*xI2(mt**2,0d0,mw**2,musq,ep))
!      -    /(12.*beta*(-1 + beta*z)*(1 + beta*z)**2) + 
!      -  ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -     (2*beta**4*s*z*(-2 + z**2) + 4*beta*s*(-1 + z**2) - 
!      -       (4*MW**2 + s)*z*(-5 + 3*z**2) + 
!      -       beta**2*z*(8*MW**2*(-2 + z**2) + 3*s*(-1 + z**2)))
!      -      *xI2(s,0d0,0d0,musq,ep))/(12.*beta*(-1 + beta**2*z**2)) - 
!      -  ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -     (2*s + beta*(4*MW**2 + 5*s)*z + 
!      -       beta**7*s*z*(-2 + z**2) + 
!      -       beta**3*z*(-8*MW**2 - 11*s + 
!      -          8*(MW**2 + s)*z**2) + 
!      -       beta**6*s*(1 - 6*z**2 + 3*z**4) + 
!      -       beta**2*(4*MW**2*(1 + z**2) + s*(-3 + 5*z**2)) + 
!      -       beta**5*z*(4*MW**2*(-2 + z**2) + 
!      -          s*(-2 - z**2 + 2*z**4)) + 
!      -       beta**4*(4*MW**2*(-3 + z**4) + 
!      -          s*(-4 - 3*z**2 + 5*z**4)))*
!      -    xI2(t,0d0,mw**2,musq,ep))/
!      -   (6.*(-1 + beta*z)*(1 + beta*z)**2*
!      -     (1 + beta**2 + 2*beta*z)) - 
!      -  ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -     (16*MW**4 + 8*beta*MW**2*s*(beta + z) + 
!      -       s**2*(3 + beta*
!      -           (2*z + beta*
!      -              (-2 + beta**2 + 2*beta*z + 2*z**2))))*
!      -  xI3(0._dp,0._dp,s,0d0,0d0,0d0,musq,ep))
!      -  /(12.*(-1 + beta**2*z**2))+
!      -  ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -     (16*MW**4 + 8*beta*MW**2*s*(beta + z) + 
!      -       s**2*(3 + beta*
!      -           (2*z + beta*
!      -              (-2 + beta**2 + 2*beta*z + 2*z**2))))*
!      -  xI3(0._dp,mt**2,t,0d0,0d0,mw**2,musq,ep))
!      -    /(12.*(-1 + beta*z)) + 
!      -  ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -     ((4*MW**2 + s)**2*z*(5 - 3*z**2) + 
!      -       2*beta**6*s**2*z*(-4 + z**2) + 
!      -       4*beta*s*(4*MW**2 + s)*(-1 + z**2) - 
!      -       4*beta**5*s**2*(1 + z**2) + 
!      -       beta**4*s*z*(16*MW**2*(-4 + z**2) - 
!      -          3*s*(-3 + z**2)) - 
!      -       8*beta**3*s*(s*z**2 + 2*MW**2*(1 + z**2)) + 
!      -       2*beta**2*z*(16*MW**4*(-4 + z**2) + 
!      -          4*MW**2*s*(-1 + z**2) + s**2*(-11 + 2*z**2)))*
!      -   xI3(s,mt**2,mt**2,0d0,0d0,mw**2,musq,ep))/
!      -   (48.*beta*(-1 + beta**2*z**2)) - 
!      -  ((-1 + beta**2)*(-7 + 9*beta*z)*
!      -     (4*MW**2 + s*(1 + beta**2 + 2*beta*z))*
!      -     (16*MW**4 + 8*beta*MW**2*s*(beta + z) + 
!      -       s**2*(3 + beta*
!      -           (2*z + beta*
!      -              (-2 + beta**2 + 2*beta*z + 2*z**2))))*
!      -   xI4(0._dp,0._dp,mt**2,mt**2,s,t,0d0,0d0,0d0,mw**2,musq,ep)
!      -  )/(48.*(-1 + beta**2*z**2))  )
!      
!      
!      
!      
!         pause


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

