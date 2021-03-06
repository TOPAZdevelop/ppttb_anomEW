
      complex*32 function hjetmass_box_pppm_0_0_s34_mhsq_s12_s134
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i4)
      t2 = zb(i3, i1)
      t3 = zb(i3, i2)
      t4 = za(i1, i2)
      t5 = za(i1, i3)
      t6 = zb(i4, i3)
      t7 = zb(i2, i1)
      t8 = zb(i4, i1)
      t9 = zb(i4, i2)
      t10 = za(i3, i4)
      t11 = za(i2, i3)
      t12 = za(i2, i4)
      t13 = t11 * t2
      t14 = t12 * t8
      t15 = t14 + t13
      t16 = t4 * t7
      t17 = t16 * t15
      t18 = t1 * t9
      t19 = t3 * t5 + t18
      t20 = mt ** 2
      t21 = t5 * t2
      t22 = t1 * t8
      t23 = t10 * t6
      t24 = t16 * (t23 + t21 + t22)
      t25 = 16 * t17 * t20 * t19 + 4 * t24 ** 2
      t25 = sqrt(t25)
      t24 = 2 * t24
      t26 = -t24 - t25
      t17 = 0.1q1 / t17
      t13 = t14 + t13
      t14 = t16 * t17
      t13 = t14 * (t21 * t13 + t22 * t13 + t23 * t13)
      t15 = t25 * t17 * t15 / 2
      t21 = 2 * t23 + 2 * t21 + 2 * t22
      t22 = -t15 + t21 - t13
      t23 = -t24 + t25
      t13 = t15 + t21 - t13
      t15 = 0.1q1 / t19
      t19 = t12 * t26 * t17 / 4
      t21 = t18 / 2
      t24 = t21 * t22 * t15
      t25 = t8 * (t1 + t19) - t24
      t24 = -t19 * t8 + t24
      t21 = t21 * t13 * t15
      t27 = t12 * t23
      t28 = t27 * t17 / 4
      t29 = -t28 * t8 + t21
      t21 = t8 * (t28 + t1) - t21
      t19 = t22 * t15 * t1 * t3 / 2 - t19 * t2
      t30 = t3 * t13
      t28 = t30 * t15 * t1 / 2 - t28 * t2
      t31 = 0.1q1 / t4
      t8 = 0.1q1 / t8
      t10 = 0.1q1 / t10
      t21 = 0.1q1 / t21
      t5 = 0.1q1 / t5
      t25 = 0.1q1 / t25
      t32 = 0.1q1 / t6
      t33 = t1 * t2
      t34 = t33 - t28
      t35 = t20 * t31
      t36 = -t35 + t7
      t37 = t8 * t7
      t6 = t37 * t6
      t38 = t3 * t7
      t39 = t2 * t8
      t40 = t1 * t31
      t41 = t40 * t5
      t42 = t5 * t8
      t43 = t42 * t7
      t44 = t39 * t32
      t45 = t15 * t1
      t46 = t22 + t13
      t47 = t1 ** 2
      t48 = t7 ** 2
      t49 = t3 * t15
      t50 = t32 * t10
      t51 = t48 * t10
      t52 = t51 * t47
      t53 = t52 * t15
      t6 = t37 * (t7 * (t50 * (t19 + t28) * t2 - t49 * (t13 * t28 * t21 
     #+ t22 * t19 * t25) + t45 * (t22 * (t20 * t2 * t26 * t17 * t25 - t3
     #) - t30) * t5) + t53 * t5 * t46 + t2 * t3 * t29 * t5 * t32 + t35 *
     # t33 * t5 * (-t18 * t13 * t15 * t21 + t50 * t29)) + t7 * (t45 * (t
     #20 * (-t42 * t38 + t39 * (t7 * t23 * t5 * t17 - t3 * t31) + t41 * 
     #(-t6 + t3)) + t39 * t38 + t43 * t34 * t9) * t21 * t13 + t44 * (t5 
     #* (t24 * t3 + t9 * (-t19 - t28) + (-t24 - t29) * t10 * t7 * t1) + 
     #t35 * (t1 * t24 * t5 - t19 - t28) * t10) + t45 * (t3 * (t39 * t36 
     #+ t20 * (-t37 + t40) * t5) + t42 * (-t19 * t7 + t33 * t36) * t9 - 
     #t6 * t41 * t20) * t25 * t22)
      t30 = 0.1q1 / t1
      t36 = t2 ** 2
      t41 = t47 * t10
      t45 = t41 * t7
      t54 = t40 * t20
      t55 = t32 * t8
      t56 = -t19 * t26 - t23 * t28
      t57 = t1 * t11
      t58 = t57 * t5 - t12
      t59 = t26 + t23
      t60 = t15 ** 2
      t61 = t26 * t24
      t62 = t23 * t29
      t63 = t23 * t13
      t64 = t26 * t22
      t65 = t3 * t12
      t66 = t65 * t59
      t67 = t9 * t12
      t68 = t26 * t22 ** 2
      t69 = t17 * t7
      t70 = t13 ** 2
      t71 = t5 * t2
      t72 = t1 * t5
      t4 = t69 * (t23 * (t52 * t42 * t70 * t60 + t71 * t32 * (t57 * t39 
     #* t7 * t10 - t28) + t7 * (t1 * (-t50 + t21) * t8 * t36 + t39 * (t2
     #1 * (t65 - t28) + (t67 * t21 - 1) * t5 * t1) - t1 * t21 * t5 * (t6
     #5 + t28)) * t15 * t13) + t26 * (-t19 * t2 * t5 * t32 + (t1 * (-t50
     # + t25) * t8 * t36 + t39 * (t25 * (t65 - t19) + (t67 * t25 - 1) * 
     #t5 * t1) - t72 * t25 * (t65 + t19)) * t15 * t22 * t7 + t42 * t4 * 
     #t25 * (t33 - t19) * t15 * t22 * t48)) + t69 * (t68 * t52 * t42 * t
     #60 + t55 * (t5 * (t67 * t59 + t61 + t62) + (t26 * t58 - t27) * t10
     # * t7) * t36 + t5 * (-t66 + (t8 * (t11 * t56 + t12 * (t62 + t61)) 
     #+ (t64 + t63) * t15 * t47) * t10 * t7) * t32 * t2 + t63 * t42 * t4
     #8 * t15 * t34 * t21 * t4 + t50 * t35 * t12 * t2 * (t59 * t5 * t1 +
     # t39 * t59))
      t27 = t56 * t30
      t33 = t59 * t2
      t34 = t23 ** 2
      t52 = t26 ** 2
      t56 = t34 * t13
      t57 = t52 * t22
      t59 = t57 + t56
      t34 = t52 + t34
      t52 = t40 * t58
      t58 = t23 * t70
      t14 = t69 * (-t42 * t14 * t12 * t30 * t32 * t34 * t36 + t39 * t10 
     #* t15 * t7 * (t64 * (t32 * (t72 * (-t67 - t24) + t19 + t65) + t52)
     # + t63 * (t32 * (t72 * (-t67 - t29) + t28 + t65) + t52)) - t38 * t
     #60 * t47 * t5 * (t58 * t21 + t68 * t25) - t55 * t12 * t17 * t30 * 
     #t34 * t2 * t36)
      t24 = t5 * (t16 * t8 - t1) + t39
      t29 = t58 + t68
      t18 = t17 * t48 * t2 * (-t50 * t5 * t34 * t17 * t12 ** 2 + (t44 * 
     #t11 * t10 * t5 * t34 + t56 * t15 * t24 * t21 + t57 * t24 * t25 * t
     #15) * t17 * t12 + t50 * t8 * t60 * t1 * (t18 * t29 * t5 - t29 * t3
     #))
      ret = -t51 * t32 * t17 ** 2 * t15 * t12 * t2 * (t39 * t59 - t72 * 
     #t59) / 16 + t14 / 4 - t18 / 8 - 3 * t43 * t20 * t17 * t36 * t32 * 
     #(t26 + t23) + 8 * t32 * t2 * (t45 * t35 * t5 + t39 * (t35 * (-t1 *
     # t7 * t10 + t3) - t38)) - 2 * t6 - 4 * t2 * (t5 * (t32 * (t7 * (-t
     #3 * (t20 * t8 + t1) + t45) + t54 * t3) - t54 * t37) - t38 * t8 * t
     #32 * (t19 + t28) * t30) - 4 * t55 * t1 * ((-t35 + t7) * t5 * t9 - 
     #t51) * t36 + 4 * t53 * t42 * t35 * (t22 + t13) + t4 / 2 + t7 * t2 
     #* (t69 * t15 * t47 * t5 * (t63 * t21 + t64 * t25) + t55 * (t17 * (
     #t2 * (t66 * t30 + t27 + t33) + t16 * (t33 + t27) * t5) + t20 * (-t
     #41 * t9 * t15 * t5 * t46 * t31 + t40 * (t49 * t46 + t71 * (-t26 - 
     #t23) * t17 * t11) * t10)))

      hjetmass_box_pppm_0_0_s34_mhsq_s12_s134 = ret/32q0/(0,1q0)
      return

      end function
