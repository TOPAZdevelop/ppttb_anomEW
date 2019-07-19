      real(dp):: delg1_z,delg1_g,lambda_g,lambda_z,
     & h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam,
     & delk_g,delk_z,tevscale
      real(dp):: h1tZ,h2tZ,h3tZ,h4tZ,h1tgam,h2tgam,h3tgam,h4tgam
      logical:: anomtgc
      real(dp) :: hitZ(4), hitgam(4)
      real(dp) :: C_phiq_333,C_phiu_33,vev,Lambda_BSM,
     & coupl_gvt,coupl_gat,coupl_gw,coupl_rest,
     & coupl_gvt_sq,coupl_gat_sq,coupl_gw_sq,
     & coupl_vol2,coupl_vol4
      common/anomcoup/delg1_z,delg1_g,lambda_g,lambda_z,delk_g,delk_z,
     & h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam,
     & tevscale,hitZ,hitgam,
     & C_phiq_333,C_phiu_33,vev,Lambda_BSM,
     & coupl_gvt,coupl_gat,coupl_gw,coupl_rest,
     & coupl_gvt_sq,coupl_gat_sq,coupl_gw_sq,
     & coupl_vol2,coupl_vol4,
     & anomtgc
!$omp threadprivate(/anomcoup/)

