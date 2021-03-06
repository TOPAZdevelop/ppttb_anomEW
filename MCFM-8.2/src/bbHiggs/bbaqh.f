*********************************************************
*AUTHOR: FABIO MALTONI                                  *
*DATE  : 10/10/2001                                     *
*NOTES : PROGRAM GENERATED BY BHIGGS.M                  *
*********************************************************

      function  BBAQH(I1,I2,I3,I4)                                    
      implicit none
      include 'types.f'
      real(dp)::  BBAQH
* ---------------------------------------------------------------------
* returns the matrix element squared for the process                   
*                                                                      
*          0 -> bbar1 b2 a3 q4 h                                       
*                                                                      
* All momenta outgoing.                                                
* No averaging is performed for initial spins or colors.               
* No STRONG coupling is included.                                      
* No Bottom-Higgs coupling is included.                                
* ---------------------------------------------------------------------

      INCLUDE 'constants.f'
      INCLUDE 'mxpart.f'
      INCLUDE 'sprods_com.f'
      INCLUDE 'zprods_com.f'

*     ARGUMENTS 
      integer:: I1,I2,I3,I4

*     LOCAL 
      real(dp)::  T134,T234
      complex(dp):: AMP(2)
      integer:: L


      t134=s(i1,i3)+s(i1,i4)+s(i3,i4)
      t234=s(i2,i3)+s(i2,i4)+s(i3,i4)

      amp(1)=
     &   (-two*(za(i3,i2)*zb(i2,i1) + za(i3,i4)*zb(i4,i1))*zb(i4,i2))/
     -   (T234*s(i3,i4)) + (two*zb(i1,i4)*
     -     (za(i1,i3)*zb(i2,i1) - za(i4,i3)*zb(i4,i2)))/
     -   (T134*s(i3,i4))

      amp(2)=
     &   (-two*(za(i4,i2)*zb(i2,i1) + za(i4,i3)*zb(i3,i1))*zb(i3,i2))/
     -   (T234*s(i3,i4)) + (two*zb(i3,i1)*
     -     (za(i4,i1)*zb(i2,i1) - za(i4,i3)*zb(i3,i2)))/
     -   (T134*s(i3,i4))


       

      BBAQH=0._dp
      do l=1,2                                         
         BBAQH=BBAQH+abs(amp(l))**2                  
      enddo                                            

      BBAQH=two*BBAQH*V/four                            

      END
