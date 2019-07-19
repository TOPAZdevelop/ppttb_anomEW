      subroutine qqb_QQb_mix_v(p,msq)
c--- Exact electroweak corrections to heavy quark pair production,
c--- including mixed QCD-weak LO process
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'breit.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scheme.f'
!      include 'first.f'
      include 'noglue.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ss,beta,z,
     &     qqbew(5),qbqew(5),ggew1,ggew2,m2(-nf:nf,-nf:nf)
      integer j,k
      
      scheme = 'dred'

!      if(first) then
!         first=.false.
!         write(6,*) 'Heavy Quark mass:', mass2
!      endif

      call dotem(4,p,s)

      ss = s(1,2)
      beta = sqrt(1._dp-4._dp*mt**2/ss)
      z = (s(1,3)-s(2,3))/ss/beta

         print *, "fixing kinematics for comparison"
         ss=   173308.20441330347d0
         beta=  0.55464632643236811d0
         z=  -0.61536288270299611d0      
      
c--- avoid calculating unnecessarily by checking input file flags
      if (ggonly) then
        qqbew(:)=zip
      else
        call qqbQQb_ew_oneloop(qqbew,ss,beta,z)
        call qqbQQb_ew_oneloop(qbqew,ss,beta,-z)
      endif
      
      if (noglue .or. omitgg) then
        ggew1=zip
        ggew2=zip
      else
        call ggQQb_ew_oneloop(ggew1,ss,beta,z)
        call ggQQb_ew_oneloop(ggew2,ss,beta,-z)
      endif
      
      call qqb_QQb(p,m2)

      msq = 0._dp

      do j=-nf,nf
         k=-j
         if((j == 0) .and. (k == 0)) then
            msq(j,k) = ggew1 + ggew2
         else if((j > 0) .and. (k < 0)) then
            msq(j,k) = qqbew(j)
         else if((j < 0 ) .and. (k > 0)) then
            msq(j,k) = qbqew(k)
         end if
      end do

      msq = msq*m2!/gsq**2*16._dp*pi**2!*0.01_dp

      end subroutine qqb_QQb_mix_v

      
      
      
      
      subroutine ResetEWCouplings(sw2,gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4) ! ( variables are in common block of anomcoup.f and set in mdata.f )
      implicit none 
      include 'types.f'      
      include 'anomcoup.f'      
      real(8) :: gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,sw2
      real(8) :: gvtBSM,gatBSM,gwBSM,vol2,vol4
      
      
        vol2 = (vev/Lambda_BSM)**2 * coupl_vol2 
        vol4 = (vev/Lambda_BSM)**4 * coupl_vol4
      
        gvtBSM = (C_phiq_333 - 0.5d0*C_phiu_33)/2d0/dsqrt(sw2*(1d0-sw2))
        gvt_sq = gvt**2 + 2*gvt*vol2*gvtBSM + vol4*gvtBSM**2
        gvt    = gvt + vol2*gvtBSM
        
        gatBSM = (C_phiq_333 + 0.5d0*C_phiu_33)/2d0/dsqrt(sw2*(1d0-sw2))
        gat_sq = gat**2 + 2*gat*vol2*gatBSM + vol4*gatBSM**2
        gat    = gat + vol2*gatBSM
                
        gwBSM = (C_phiq_333)*0.5d0/dsqrt(2d0)/dsqrt(sw2)
        gw_sq = gw**2 + 2*gw*vol2*gwBSM + vol4*gwBSM**2           
        gw    = gw  + vol2*gwBSM
       
        gvt    = gvt    * coupl_gvt
        gat    = gat    * coupl_gat
        gw     = gw     * coupl_gw
        gvt_sq = gvt_sq * coupl_gvt_sq
        gat_sq = gat_sq * coupl_gat_sq
        gw_sq  = gw_sq  * coupl_gw_sq
        g_rest = 1d0    * coupl_rest
      
      return
      end subroutine


