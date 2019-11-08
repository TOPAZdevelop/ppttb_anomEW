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
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scheme.f'
!      include 'first.f'
      include 'noglue.f'
      include 'anomcoup.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ss,beta,z,
     &     qqbew(5),qbqew(5),ggew1,ggew2,m2(-nf:nf,-nf:nf)
      integer j,k
      real(dp) :: Mbbar_ttbar,Mbarb_ttbar,Mqqbar_ttbar(1:2),Mqbarq_ttbar(1:2)
      real(dp) :: aemmz,mw,mz,shat,that,PMG(1:4,1:4),ResMG,till_check
      real(dp) :: sw2,cw2,gvt,gat,gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4
      real(dp) :: gab,gvb,gvu,gau,gw_sq_gat,gw_sq_gvt,gw_qu,gatgvt
      common/em/aemmz
      
      scheme = 'dred'

!      if(first) then
!         first=.false.
!         write(6,*) 'Heavy Quark mass:', mass2
!      endif

      call dotem(4,p,s)

      ss = s(1,2)
      beta = sqrt(1._dp-4._dp*mt**2/ss)
      z = (s(1,3)-s(2,3))/ss/beta
      mw = wmass
      mz = zmass
      
!       print *, "fixing kinematics for comparison in qqb_QQb_mix_v.f"
!       ss=   173308.20441330347d0
!       beta=  0.55464632643236811d0
!       z=  +0.61536288270299611d0      
      
c--- avoid calculating unnecessarily by checking input file flags
      if (ggonly) then
        qqbew(:)=zip
      else
        call qqbQQb_ew_oneloop(qqbew,ss,beta,z)
        call qqbQQb_ew_oneloop(qbqew,ss,beta,-z)
      endif
!       print *, "virtqq/born", qqbew
!       print *, "virtqq/born", qbqew     


      if (noglue .or. omitgg) then
        ggew1=zip
        ggew2=zip
      else
        call ggQQb_ew_oneloop(ggew1,ss,beta,z)
        call ggQQb_ew_oneloop(ggew2,ss,beta,-z)
      endif
      
!        print *, "virtgg/born", ggew1
!        print *, "virtgg/born", ggew2

!     this is for  Cpq3=  13.000000000000000   Cpu=  5.0000000000000000  for vol2 >0 , vol4>0     
!        print *, "ratio1", qqbew/(/-1.4757774020536695d-002, 
!      &                           -4.2609888664684847d-002,  
!      &                           -1.4757774020536695d-002,  
!      &                           -4.2609888664684847d-002,  
!      &                           -1.4757774020536695d-002/)  
!        print *, "ratio2", qbqew/(/-5.8308261197494958d-002, 
!      &                           -1.7475701014336532d-002,  
!      &                           -5.8308261197494958d-002,  
!      &                           -1.7475701014336532d-002,  
!      &                           -5.8308261197494958d-002/)  
!        print *, "ratio3", ggew1/(-5.4430573399833621d-002)  
!        print *, "ratio4", ggew2/(-0.13033209790126796d0  )

!     this is for  Cpq3= -13.000000000000000   Cpu=  5.0000000000000000       for vol2 >0 , vol4>0     
!      print *, "ratio1", qbqew/(/-1.2714991939120214d-003, 
!     &                           2.8438348797801157d-002,  
!     &                           -1.2714991939120214d-003,  
!     &                           2.8438348797801157d-002,  
!     &                           -1.2714991939120214d-003/)  
!      print *, "ratio2", qqbew/(/3.6961212336970986d-002, 
!     &                           6.3731961832221243d-003,  
!     &                           3.6961212336970986d-002,  
!     &                           6.3731961832221243d-003,  
!     &                           3.6961212336970986d-002/)  
!      print *, "ratio3", ggew1/(1.4794182518949119d-002  )  
!      print *, "ratio4", ggew2/( -5.0709055364319311d-004 )  

!     this is for  Cpq3= +13.000000000000000   Cpu= -5.0000000000000000      for vol2 >0 , vol4>0      
!      print *, "ratio1", qbqew/(/3.4858339486468605E-003, 
!     &                           -3.8101006914137314E-002,  
!     &                           3.4858339486468605E-003,  
!     &                           -3.8101006914137314E-002,  
!     &                           3.4858339486468605E-003/)  
!      print *, "ratio2", qqbew/(/-5.5792191441357768E-002, 
!     &                            -3.8900230589952360E-003  ,  
!     &                           -5.5792191441357768E-002,  
!     &                           -3.8900230589952360E-003  ,  
!     &                           -5.5792191441357768E-002/)  
!      print *, "ratio3", ggew1/( -7.6268737890503271E-002 )  
!      print *, "ratio4", ggew2/( -3.2563095544884020E-002 )  

!     this is for  Cpq3= -13.000000000000000   Cpu= -5.0000000000000000     for vol2 >0 , vol4>0      
!      print *, "ratio1", qbqew/(/ 2.1839339216189221E-003, 
!     &                            1.8159055694696054E-002 ,  
!     &                           2.1839339216189221E-003,  
!     &                            1.8159055694696054E-002 ,  
!     &                            2.1839339216189221E-003/)  
!      print *, "ratio2", qqbew/(/2.4689107239455523E-002, 
!     &                            5.1706992849107677E-003  ,  
!     &                          2.4689107239455523E-002,  
!     &                          5.1706992849107677E-003 ,  
!     &                          2.4689107239455523E-002/)  
!      print *, "ratio3", ggew1/( 3.8086211379597447E-004 )  
!      print *, "ratio4", ggew2/(  -5.6888644920950167E-003 )  




! !     this is for  Cpq3=  13.000000000000000   Cpu=  5.0000000000000000    for vol2 >0 , vol4==0        
!        print *, "ratio1", qqbew/(/ -8.5347870772868959d-003, 
!      &                           -3.6386901721435015d-002,  
!      &                            -8.5347870772868959d-003,  
!      &                           -3.6386901721435015d-002,  
!      &                            -8.5347870772868959d-003/)  
!        print *, "ratio2", qbqew/(/-5.2085274254245140E-002, 
!      &                            -1.1252714071086719E-002,  
!      &                           -5.2085274254245140E-002,  
!      &                            -1.1252714071086719E-002,  
!      &                           -5.2085274254245140E-002/)  
!        print *, "ratio3", ggew1/(-3.1390539346056533E-002)  
!        print *, "ratio4", ggew2/( -7.6430442041101998E-002  )
! 
! 
! 
!       pause     


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
      


            
! ! ! ! ! ! ! ! !             swapping momenta because the ME below is calculated with p1+p2-->p3+p4
        s(1,3) = -s(1,3)
        s(1,4) = -s(1,4)
        s(2,3) = -s(2,3)
        s(2,4) = -s(2,4)
! 
!       print*,s(1,2),+2*(p(1,4)*p(2,4)-p(1,1)*p(2,1)-p(1,2)*p(2,2)-p(1,3)*p(2,3))
!       print*,s(1,3),-2*(p(1,4)*p(3,4)-p(1,1)*p(3,1)-p(1,2)*p(3,2)-p(1,3)*p(3,3))
!       print*,s(1,4),-2*(p(1,4)*p(4,4)-p(1,1)*p(4,1)-p(1,2)*p(4,2)-p(1,3)*p(4,3))
!       print*,s(2,3),-2*(p(3,4)*p(2,4)-p(3,1)*p(2,1)-p(3,2)*p(2,2)-p(3,3)*p(2,3))
!       print*,s(2,4),-2*(p(4,4)*p(2,4)-p(4,1)*p(2,1)-p(4,2)*p(2,2)-p(4,3)*p(2,3))
!       print*,s(3,4),+2*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
        
!         print *, "test",p(1,4),p(1,1:3)
!         print *, "test",p(2,4),p(2,1:3)
!         print *, "test",p(3,4),p(3,1:3)
!         print *, "test",p(4,4),p(4,1:3)
        
        
       ss = s(1,2)
       beta = sqrt(1._dp-4._dp*mt**2/ss)
       z = (s(1,3)-s(2,3))/ss/beta        
       sw2 = xw
       cw2 = 1._dp - sw2
       gvt = ((1d0/2d0)-2._dp*sw2*(2d0/3d0))/2._dp/sqrt(sw2*cw2)
       gat = (1d0/2d0)/2._dp/sqrt(sw2*cw2)
       gvu = gvt
       gau = gat
       gvb = ((-1d0/2d0)-2._dp*sw2*(-1d0/3d0))/2._dp/sqrt(sw2*cw2)
       gab = (-1d0/2d0)/2._dp/sqrt(sw2*cw2)       
       gw = 0.5_dp/sqrt(2._dp*sw2)
       call ResetEWCouplings2(sw2,gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4,gatgvt,gw_qu,gw_sq_gvt,gw_sq_gat)


 
!      calculating b bbar --> W(t-channel)  --> t tbar interference with QCD tree
        Mbbar_ttbar=(-8*(esq/4d0/pi)*as*Pi**2*8*gw_sq*
     -    (mt**2*ss*(mt**2 + 2*MW**2 + s(1,3)) + 
     -      (2*MW**2*ss - (mt**2 + 2*MW**2)*s(1,3))*s(2,3)))/
     -  (9*MW**2*ss*(-mt**2 + MW**2 + s(1,3)))*4d0

        Mbarb_ttbar =(-8*(esq/4d0/pi)*as*Pi**2*8*gw_sq*
     -    (mt**2*ss*(mt**2 + 2*MW**2 + s(2,3)) + 
     -      (2*MW**2*ss - (mt**2 + 2*MW**2)*s(2,3))*s(1,3)))/
     -  (9*MW**2*ss*(-mt**2 + MW**2 + s(2,3)))*4d0

     
     
        msq(+5,-5) = msq(+5,-5)  + Mbbar_ttbar
        msq(-5,+5) = msq(-5,+5)  + Mbarb_ttbar

     
!         shat = ss
!         that = -s(1,3) + mt**2
!         till_check = (-2*esq*4*pi*as*(mt**6 + mt**4*(2*MW**2 + shat - 2*that) 
!      & + 2*MW**2*(shat + that)**2 + mt**2*(that**2 - 2*MW**2*(shat + 2*that)))
!      &  *8*gw_sq)/(9*MW**2*shat*(MW**2-that))
!         
! 
! ! ! ! ! ! ! ! ! ! ! ! !       
! 
!       PMG(1:4,1) = -(/p(1,4),p(1,1),p(1,2),p(1,3)/)
!       PMG(1:4,2) = -(/p(2,4),p(2,1),p(2,2),p(2,3)/)
!       PMG(1:4,3) =  (/p(3,4),p(3,1),p(3,2),p(3,3)/)
!       PMG(1:4,4) =  (/p(4,4),p(4,1),p(4,2),p(4,3)/)
!       call coupsm(0)
!       call SBBB_TTB(PMG,ResMG) ! careful the MG object file for the alpha*alpha_s interference has been overwritten 
! 
!       print *, "Mbbar_ttbar ",Mbbar_ttbar
!       print *, "Mbarb_ttbar ",Mbarb_ttbar
!       print *, "ratio1",   -4.8668147577302694d0/Mbbar_ttbar
!       print *, "ratio2",   -2.3832594037514752d0/Mbarb_ttbar
!       
! !       print *, till_check
! !       print *, till_check/Mbbar_ttbar
! !       print *, ResMG
!       print *, "ratio MG SM",ResMG/Mbbar_ttbar      
!       
!       pause

        
        
        
        
        
!      calculating q qbar --> Z/gamma(s-channel)  --> t tbar
       

!      dn quarks 
       Mqqbar_ttbar(1) = (esq/4d0/pi)**2*(4d0*pi)**2*(
     - (-1d0/3d0)**2*(2d0/3d0)**2*(2d0-beta**2*(1-z**2))
     - +ss**2/(ss-MZ**2)**2 *( (gvb**2+gab**2)*
     - (gvt_sq*(2-beta**2*(1d0-z**2))+gat_sq*beta**2*
     - (1d0+z**2))-8*gvb*gab*gatgvt*beta*z )
     - +2*(-1d0/3d0)*(2d0/3d0)*ss/(ss-MZ**2)*(
     - gvb*gvt*(2d0-beta**2*(1d0-z**2))-2d0*gab*gat*beta*z)
     -  )
     
!      up quarks 
       Mqqbar_ttbar(2) = (esq/4d0/pi)**2*(4d0*pi)**2*(
     - (2d0/3d0)**2*(2d0/3d0)**2*(2d0-beta**2*(1d0-z**2))
     - +ss**2/(ss-MZ**2)**2 *( (gvu**2+gau**2)*
     - (gvt_sq*(2d0-beta**2*(1d0-z**2))+gat_sq*beta**2*
     - (1d0+z**2))-8*gvu*gau*gatgvt*beta*z )
     - +2*(2d0/3d0)*(2d0/3d0)*ss/(ss-MZ**2)*(
     - gvu*gvt*(2d0-beta**2*(1d0-z**2))-2d0*gau*gat*beta*z)
     -  )

       
       z = -z 
!      dn quarks 
       Mqbarq_ttbar(1) = (esq/4d0/pi)**2*(4d0*pi)**2*(
     - (-1d0/3d0)**2*(2d0/3d0)**2*(2d0-beta**2*(1-z**2))
     - +ss**2/(ss-MZ**2)**2 *( (gvb**2+gab**2)*
     - (gvt_sq*(2-beta**2*(1d0-z**2))+gat_sq*beta**2*
     - (1d0+z**2))-8*gvb*gab*gatgvt*beta*z )
     - +2*(-1d0/3d0)*(2d0/3d0)*ss/(ss-MZ**2)*(
     - gvb*gvt*(2d0-beta**2*(1d0-z**2))-2d0*gab*gat*beta*z)
     -  )
     
!      up quarks 
       Mqbarq_ttbar(2) = (esq/4d0/pi)**2*(4d0*pi)**2*(
     - (2d0/3d0)**2*(2d0/3d0)**2*(2d0-beta**2*(1d0-z**2))
     - +ss**2/(ss-MZ**2)**2 *( (gvu**2+gau**2)*
     - (gvt_sq*(2d0-beta**2*(1d0-z**2))+gat_sq*beta**2*
     - (1d0+z**2))-8*gvu*gau*gatgvt*beta*z )
     - +2*(2d0/3d0)*(2d0/3d0)*ss/(ss-MZ**2)*(
     - gvu*gvt*(2d0-beta**2*(1d0-z**2))-2d0*gau*gat*beta*z)
     -  )
     
     
        msq(+1,-1) = msq(+1,-1) + Mqqbar_ttbar(1)  
        msq(-1,+1) = msq(-1,+1) + Mqbarq_ttbar(1)
        msq(+3,-3) = msq(+3,-3) + Mqqbar_ttbar(1)
        msq(-3,+3) = msq(-3,+3) + Mqbarq_ttbar(1)  

        msq(+2,-2) = msq(+2,-2) + Mqqbar_ttbar(2)
        msq(-2,+2) = msq(-2,+2) + Mqbarq_ttbar(2)
        msq(+4,-4) = msq(+4,-4) + Mqqbar_ttbar(2)
        msq(-4,+4) = msq(-4,+4) + Mqbarq_ttbar(2)  

        
        
        
        
!       print *, ""
!       print *, "uub-->Z/ga->ttb",Mqqbar_ttbar(1)
!       print *, "ddb-->Z/ga->ttb",Mqqbar_ttbar(2)
! 
!       print *, "ubu-->Z/ga->ttb",Mqbarq_ttbar(1)
!       print *, "dbd-->Z/ga->ttb",Mqbarq_ttbar(2)
        
        
       z = -z ! undo       

! ! ! ! ! ! ! ! ! ! ! ! !       
!       PMG(1:4,1) = -(/p(1,4),p(1,1),p(1,2),p(1,3)/)
!       PMG(1:4,2) = -(/p(2,4),p(2,1),p(2,2),p(2,3)/)
!       PMG(1:4,3) =  (/p(3,4),p(3,1),p(3,2),p(3,3)/)
!       PMG(1:4,4) =  (/p(4,4),p(4,1),p(4,2),p(4,3)/)
!       call coupsm(0)
! !       call SUUB_TTB(PMG,ResMG)
!       call SDBD_TTB(PMG,ResMG)
!       print *, "MZ, sw2",MZ,xw
! !       print *, "Mqqar_ttbar ",Mqqbar_ttbar(2)
! !       print *, "ratio MG SM",ResMG/Mqqbar_ttbar(2)   
!       print *, "Mqqar_ttbar ",Mqbarq_ttbar(1)
!       print *, "ratio MG SM",ResMG/Mqbarq_ttbar(1)
!       pause
! ! ! ! ! ! ! ! ! ! ! ! ! !    
                

                
                
                
                
! !      b quarks 
       Mqqbar_ttbar(1) = (esq/4d0/pi)**2*(4d0*pi)**2*(     ! Z/gam squared  (s-channel)
     - (-1d0/3d0)**2*(2d0/3d0)**2*(2d0-beta**2*(1-z**2))
     - +ss**2/(ss-MZ**2)**2 *( (gvb**2+gab**2)*
     - (gvt_sq*(2-beta**2*(1d0-z**2))+gat_sq*beta**2*
     - (1d0+z**2))-8*gvb*gab*gatgvt*beta*z )
     - +2*(-1d0/3d0)*(2d0/3d0)*ss/(ss-MZ**2)*(
     - gvb*gvt*(2d0-beta**2*(1d0-z**2))-2d0*gab*gat*beta*z)
!      
     - +(4*gw_qu*ss**2*(16*(-1+beta**2)**2*MW**2*ss      ! W squared (t-channel) 
     - +(-1+beta**2)**2*ss**2*(-1-beta*z)**2
     - +64*MW**4*(1-beta*z)**2))
     - /(4*MW**4*(4*MW**2+ss*(1+beta**2+2*beta*z))**2)
! 
     -  +1d0/3d0*(gw_sq*(-1d0/3d0)*(2d0/3d0)*ss*        ! interferences beteen s and t channel (photon)
     -  ((-1 + beta**2)*ss*(2 + 2*beta*z + beta**2*(-1 + z**2))+ 
     -    8*MW**2*(-2 + beta*(beta + 2*z - beta*z**2))))/
     - (MW**2*(4*MW**2 + ss*(1 + beta**2 + 2*beta*z))) 
!      
     - + ((gab + gvb)*ss**2* ! interferences beteen s and t channel (Z)
     - (2*gw_sq_gvt*(8*MW**2+ss)+2*beta**3*(gw_sq_gat-gw_sq_gvt)*ss*z- 
     - 2*beta*(8*(gw_sq_gat+gw_sq_gvt)*MW**2+(gw_sq_gat-gw_sq_gvt)*ss)*z+
     - beta**4*ss*(gw_sq_gat+gw_sq_gvt+(gw_sq_gat-gw_sq_gvt)*z**2)+
     - beta**2*(gw_sq_gat*(8*MW**2 - ss)*(1 + z**2) + 
     - gw_sq_gvt*(ss*(-3 + z**2) + 8*MW**2*(-1 + z**2)))))/
     - (3*MW**2*(MZ**2 - ss)*(4*MW**2 + ss*(1 + beta**2 + 2*beta*z)))   
     -  )        
   
     
       z = -z 
!      b quarks 
       Mqbarq_ttbar(1) =  (esq/4d0/pi)**2*(4d0*pi)**2*(     ! Z/gam squared  (s-channel)
     - (-1d0/3d0)**2*(2d0/3d0)**2*(2d0-beta**2*(1-z**2))
     - +ss**2/(ss-MZ**2)**2 *( (gvb**2+gab**2)*
     - (gvt_sq*(2-beta**2*(1d0-z**2))+gat_sq*beta**2*
     - (1d0+z**2))-8*gvb*gab*gatgvt*beta*z )
     - +2*(-1d0/3d0)*(2d0/3d0)*ss/(ss-MZ**2)*(
     - gvb*gvt*(2d0-beta**2*(1d0-z**2))-2d0*gab*gat*beta*z)
!      
     - +(4*gw_qu*ss**2*(16*(-1+beta**2)**2*MW**2*ss      ! W squared (t-channel)
     - +(-1+beta**2)**2*ss**2*(-1-beta*z)**2
     - +64*MW**4*(1-beta*z)**2))
     - /(4*MW**4*(4*MW**2+ss*(1+beta**2+2*beta*z))**2)
! 
     -  +1d0/3d0*(gw_sq*(-1d0/3d0)*(2d0/3d0)*ss*        ! interferences beteen s and t channel (photon)
     -  ((-1 + beta**2)*ss*(2 + 2*beta*z + beta**2*(-1 + z**2))+ 
     -    8*MW**2*(-2 + beta*(beta + 2*z - beta*z**2))))/
     - (MW**2*(4*MW**2 + ss*(1 + beta**2 + 2*beta*z)))  
!      
     - + ((gab + gvb)*ss**2* ! interferences beteen s and t channel (Z)
     - (2*gw_sq_gvt*(8*MW**2+ss)+2*beta**3*(gw_sq_gat-gw_sq_gvt)*ss*z- 
     - 2*beta*(8*(gw_sq_gat+gw_sq_gvt)*MW**2+(gw_sq_gat-gw_sq_gvt)*ss)*z+
     - beta**4*ss*(gw_sq_gat+gw_sq_gvt+(gw_sq_gat-gw_sq_gvt)*z**2)+
     - beta**2*(gw_sq_gat*(8*MW**2 - ss)*(1 + z**2) + 
     - gw_sq_gvt*(ss*(-3 + z**2) + 8*MW**2*(-1 + z**2)))))/
     - (3*MW**2*(MZ**2 - ss)*(4*MW**2 + ss*(1 + beta**2 + 2*beta*z)))       
     -  )        
       
       z = -z ! undo       
       
       
       
        msq(+5,-5) = msq(+5,-5) + Mqqbar_ttbar(1) 
        msq(-5,+5) = msq(-5,+5) + Mqbarq_ttbar(1)        


!       print *, "b bb-->Z/ga->ttb",Mqqbar_ttbar(1)
!       print *, "bb b-->Z/ga->ttb",Mqbarq_ttbar(1)
!       pause
! 
! ! ! ! ! ! ! ! ! ! ! ! ! !       
!       PMG(1:4,1) = -(/p(1,4),p(1,1),p(1,2),p(1,3)/)
!       PMG(1:4,2) = -(/p(2,4),p(2,1),p(2,2),p(2,3)/)
!       PMG(1:4,3) =  (/p(3,4),p(3,1),p(3,2),p(3,3)/)
!       PMG(1:4,4) =  (/p(4,4),p(4,1),p(4,2),p(4,3)/)
!       call coupsm(0)
!       call SBBB_TTB(PMG,ResMG)
!       print *, "MZ, sw2",MZ,xw
!       print *, "Mqqar_ttbar ",Mqqbar_ttbar(1)
!       print *, "ratio MG SM",ResMG/Mqqbar_ttbar(1)
!       pause
! ! ! ! ! ! ! ! ! ! ! ! ! !    
        
        s(1,3) = -s(1,3)   ! undo the swap from above 
        s(1,4) = -s(1,4)
        s(2,3) = -s(2,3)
        s(2,4) = -s(2,4)
        
        


      end subroutine qqb_QQb_mix_v

      
      
      
      
      subroutine ResetEWCouplings(sw2,gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4) ! ( variables are in common block of anomcoup.f and set in mdata.f )
      implicit none 
      include 'types.f'      
      include 'anomcoup.f'      
      real(8) :: gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,sw2
      real(8) :: gvtBSM,gatBSM,gwBSM,vol2,vol4
      
      
        vol2 = (vev/Lambda_BSM)**2 * coupl_vol2 
        vol4 = (vev/Lambda_BSM)**4 * coupl_vol4
      
        gvtBSM = (C_phiq_333-0.5d0*C_phiu_33)/2d0/dsqrt(sw2*(1d0-sw2))
        gvt_sq = gvt**2 + 2*gvt*vol2*gvtBSM + vol4*gvtBSM**2
        
        gatBSM = (C_phiq_333+0.5d0*C_phiu_33)/2d0/dsqrt(sw2*(1d0-sw2))
        gat_sq = gat**2 + 2*gat*vol2*gatBSM + vol4*gatBSM**2
                
        gwBSM = (C_phiq_333)*0.5d0/dsqrt(2d0)/dsqrt(sw2)
        gw_sq = gw**2 + 2*gw*vol2*gwBSM + vol4*gwBSM**2           

       
!       finally shifting the original values
        gvt = gvt + vol2*gvtBSM
        gat = gat + vol2*gatBSM
        gw  = gw  + vol2*gwBSM

       
       
        gvt    = gvt    * coupl_gvt
        gat    = gat    * coupl_gat
        gw     = gw     * coupl_gw
        gvt_sq = gvt_sq * coupl_gvt_sq
        gat_sq = gat_sq * coupl_gat_sq
        gw_sq  = gw_sq  * coupl_gw_sq
        g_rest = 1d0    * coupl_rest
      
      return
      end subroutine


      
      subroutine ResetEWCouplings2(sw2,gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,vol2,vol4,gatgvt,gw_qu,gw_sq_gvt,gw_sq_gat) ! ( variables are in common block of anomcoup.f and set in mdata.f )
      implicit none 
      include 'types.f'      
      include 'anomcoup.f'      
      real(8) :: gvt,gat,gw,gvt_sq,gat_sq,gw_sq,g_rest,sw2
      real(8) :: gvtBSM,gatBSM,gwBSM,vol2,vol4
      real(8) :: gw_sq_gat,gw_sq_gvt,gw_qu,gatgvt
        
        vol2 = (vev/Lambda_BSM)**2 * coupl_vol2 
        vol4 = (vev/Lambda_BSM)**4 * coupl_vol4
      
        gvtBSM = (C_phiq_333-0.5d0*C_phiu_33)/2d0/dsqrt(sw2*(1d0-sw2))
        gvt_sq = gvt**2 + 2*gvt*vol2*gvtBSM + vol4*gvtBSM**2
        
        gatBSM = (C_phiq_333+0.5d0*C_phiu_33)/2d0/dsqrt(sw2*(1d0-sw2))
        gat_sq = gat**2 + 2*gat*vol2*gatBSM + vol4*gatBSM**2
                
        gwBSM = (C_phiq_333)*0.5d0/dsqrt(2d0)/dsqrt(sw2)
        gw_sq = gw**2 + 2*gw*vol2*gwBSM + vol4*gwBSM**2           

        ! these terms contained dim-10 couplings, i.e. vol2^4 or vol2*vol4, which have been set to zero here
        gatgvt = gat*gvt + gatBSM*gvt*vol2 + gat*gvtBSM*vol2 
     -         + gatBSM*gvtBSM*vol4
        gw_qu = gw**4 + 4*gw**3*vol2*gwBSM + 6*gw**2*vol4*gwBSM**2
        gw_sq_gvt = gvt*gw**2 + gvtBSM*gw**2*vol2 
     -            + 2*gvt*gw*gwBSM*vol2 + 2*gw*gwBSM*gvtBSM*vol4
     -            + gvt*gwBSM**2*vol4
        gw_sq_gat = gat*gw**2 + gatBSM*gw**2*vol2 
     -            + 2*gat*gw*gwBSM*vol2 + 2*gw*gwBSM*gatBSM*vol4
     -            + gat*gwBSM**2*vol4
       
       
!       finally shifting the original values
        gvt = gvt + vol2*gvtBSM
        gat = gat + vol2*gatBSM
        gw  = gw  + vol2*gwBSM

       
       
        gvt    = gvt    * coupl_gvt
        gat    = gat    * coupl_gat
        gw     = gw     * coupl_gw
        gvt_sq = gvt_sq * coupl_gvt_sq
        gat_sq = gat_sq * coupl_gat_sq
        gw_sq  = gw_sq  * coupl_gw_sq
        g_rest = 1d0    * coupl_rest
      
      return
      end subroutine


