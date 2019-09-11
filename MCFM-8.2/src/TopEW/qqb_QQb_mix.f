      subroutine qqb_QQb_mix(p,msq)
c--- Leading order matrix element for mixed QCD/Z matrix element
c--- for top pair production, q(-p1)+qbar(-p2) -->  Q(p3)+Q~(P4) 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'breit.f'
      include 'first.f'
      include 'anomcoup.f'      
      integer j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),gaq(nf),gvq(nf),
     .     gat,gvt,ss,beta,z,qqb(nf),qbq(nf),cs,sw2,cw2,mz,T3(nf),
     .     gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4

      if(first) then
         first=.false.
         write(6,*) 'Heavy Quark mass:', mass2
      endif

      
!       p= reshape( ! use this point for comparison of vol2>0 vol4=0
!      .(/
!      .    0d0,   0d0,  -190.32693206964728d0,   190.32693206964728d0,   
!      .    0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   
!      .    0d0,   0d0,   0d0,   371.75945810707674d0,  -371.75945810707674d0,   
!      .    0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   
!      .    0d0,  -192.52614214911787d0,   1066.8609386650890d0, 
!      .     -393.87339113072892d0,  -480.46140538524219d0,   0d0,   0d0,  
!      .     0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,  
!      .      0d0,  -192.52614214911787d0,  -1066.8609386650890d0,  
!      .       599.63649322229423d0,   659.75058759191234d0,   0d0,  
!      .        0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0, 
!      .          0d0   /) , (/mxpart,4/))
     
      
      sw2 = xw
      cw2 = 1._dp - sw2
      T3 = (/-0.5_dp,0.5_dp,-0.5_dp,0.5_dp,-0.5_dp/)
      gvq(1:5) = (T3(1:5)-2._dp*sw2*Q(1:5))/2._dp/sqrt(sw2*cw2)
      gaq(1:5) = T3(1:5)/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)


      gvt = gvq(2)
      gat = gaq(2)

      
c**********************************************************************************
c     MARKUS: add dim-6 operator contributions ( variables are in common block of anomcoup.f and set in mdata.f )       

        call ResetEWCouplings(sw2,gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4)!   defined in qqb_QQb_mix.f

c     END MARKUS      
c**********************************************************************************
      
      mz = zmass
      
      call dotem(4,p,s)

      ss = s(1,2)
      beta = sqrt(1._dp - 4._dp*mass2**2/ss)
      z = (s(1,3)-s(2,3))/ss/beta

      qqb = 4._dp*ss/(ss - mz**2)*(
     .     + (2._dp - beta**2*(1._dp - z**2))*gvq*gvt
     .     + 2._dp*beta*z*gaq*gat )

      qbq = 4._dp*ss/(ss - mz**2)*(
     .     + (2._dp - beta**2*(1._dp - z**2))*gvq*gvt
     .     - 2._dp*beta*z*gaq*gat )


      
      msq=zip 
      do j=-nf,nf
         k=-j
         if((j == 0) .and. (k == 0)) then
            msq(j,k)=0._dp
         else if((j > 0) .and. (k < 0)) then
            msq(j,k) = qqb(j)
         else if((j < 0 ) .and. (k > 0)) then
            msq(j,k) = qbq(k)
         end if
      end do
 
! factor of two for the interference, 1/Nc**2/4 is the color and spin average
      msq = msq*2._dp*esq*gsq/Nc**2/4._dp


!       print *, " REAL AMPLITUDE SQUARED:  msqmix="
!       write(13,*) msq
!       pause

      
      
      end subroutine qqb_QQb_mix
