
      function a6treeQLlc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6treeQLlc

c-----Tree amplitude for leading color contributions
c----- 0 -> q(p1) + qb(p2) + glu(p3) + gam(p4) + lb(p5) + l(p6)
c-----where the photon is coming from the lepton line
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: st
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t

c---- helicity stamps
c     'q+qb-g+ga+lb-l+'=4
c     'q+qb-g+ga-lb-l+'=5

c-----Multiplied by (-i), factor out sqrt(2) 
      if(st==4) then
      a6treeQLlc=
     .za(j2,j5)**2/(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
      elseif(st==5) then
      a6treeQLlc=
     .(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))**2/
     .(t(j4,j5,j6)*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
      else
      write(6,*) 'unimplemented st'
      stop
      endif
c-----done
      return
      end


      function a6treeQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6treeQLslc

c-----Tree amplitude for subleading color contributions
c----- 0 -> q(p1) + qb(p2) + glu(p3) + gam(p4) + lb(p5) + l(p6)
c-----where the photon is coming from  the lepton line
      
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: st
      real(dp):: t

c---- helicity stamps
c     'q+qb-g+ga+lb-l+'=4
c     'q+qb-g+ga-lb-l+'=5

c-----Multiplied by (-i), factor out sqrt(2) 
      if(st==4) then
      a6treeQLslc=
     .-za(j2,j5)**2/(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
      elseif(st==5) then
      a6treeQLslc= 
     .-(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))**2/
     & (t(j4,j5,j6)*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
      else 
      write(6,*) 'unimplemented st'
      stop
      endif
      return
      end



