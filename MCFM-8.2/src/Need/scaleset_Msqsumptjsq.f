      subroutine scaleset_Msqsumptjsq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+sum(ptj^2)), where M is the mass of the particle (34)
c--- and the sum is over the pt-squared of each jet in the event
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'breit.f'
      integer:: isub,oldjets
      real(dp):: p(mxpart,4),pjet(mxpart,4),mu0,pt,rcut,sumptsq
      common/rcut/rcut

      if((kcase==kW_1jet) .or.
     &   (kcase==kZ_1jet) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==khjetma) .or.
     &   (kcase==kW_2jet) .or.
     &   (kcase==kZ_2jet) .or.
     &   (kcase==kggfus2)) then

c--- first work out whether this point is real radiation or not
        if (abs(p(npart+2,4)) > 1.e-8_dp) then
          isub=0  ! real term
        else
          isub=1  ! subtraction term
        endif
      
c-- cluster jets but make sure recorded number of jets is not changed
        oldjets=jets
        call genclust2(p,rcut,pjet,isub)
      
c        write(6,*) 'partons:'
c        if (p(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,p)
c        if (p(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,p)  
c      write(6,*) 'jets:'
c        if (pjet(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,pjet)
c        if (pjet(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,pjet)
c        write(6,*)

c--- work out sum of pt-squared
        if     (jets == 1) then
          sumptsq=pt(5,pjet)**2
        elseif (jets == 2) then
          sumptsq=pt(5,pjet)**2+pt(6,pjet)**2
        elseif (jets == 3) then
          sumptsq=pt(5,pjet)**2+pt(6,pjet)**2+pt(7,pjet)**2
        else
          sumptsq=0._dp
        endif
      
c--- restore old value of jets
        jets=oldjets

        mu0=mass3**2+sumptsq
        mu0=sqrt(abs(mu0))
      
      else
        write(6,*) 'dynamicscale sqrt(M^2+sumptj^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
