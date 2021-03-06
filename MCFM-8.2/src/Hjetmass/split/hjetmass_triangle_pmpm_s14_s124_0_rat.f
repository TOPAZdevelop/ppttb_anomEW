
      complex*32 function hjetmass_triangle_pmpm_s14_s124_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = za(i2, i3)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = t2 * t3
      t6 = t1 * t4
      t7 = t5 + t6
      if ( real(t7) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t7 = cg * sqrt(t7 ** 2) + t5 + t6
      t8 = za(i1, i2)
      t9 = zb(i2, i1)
      t10 = t7 / 2
      t11 = t8 * t9
      t12 = t11 * (-t10 + t5)
      t13 = za(i2, i4)
      t14 = za(i3, i4)
      t15 = zb(i4, i2)
      t16 = zb(i4, i3)
      t17 = t2 * t9
      t10 = t8 * (t10 * t15 - t17 * t16)
      t5 = -t7 * t8 * t3 * (t14 * t15 + t17) / 2 + t5 * t11 * (t14 * t16
     # + t5 + t6)
      t2 = 0.1q1 / t2
      t5 = 0.1q1 / t5
      t6 = 0.1q1 / t16
      t8 = 0.1q1 / t8
      t11 = 0.1q1 / t14
      t9 = 0.1q1 / t9
      t10 = 0.1q1 / t10
      t15 = t3 * t5
      t16 = t14 * t15 + t10
      t17 = t3 * t6
      t18 = t9 ** 2 * t2
      ret = 16 * t18 * (t17 * (t16 * t9 * t8 * t2 * t12 ** 2 - t3 * t12 
     #* t16) + t8 * t12 * (t11 * (t12 * t10 ** 2 + t17 * t10) + t3 ** 2 
     #* t5 * (-t12 * t14 * t5 + t6)) * (t3 * za(i1, i4) + t4 * t13) * t1
     #) - 8 * t18 * t17 * t7 * t13 * t12 * t8 * (t10 * t11 + t15)

      hjetmass_triangle_pmpm_s14_s124_0_rat = ret/32q0/(0,1q0)
      return

      end function
