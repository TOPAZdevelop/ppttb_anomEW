
      complex*32 function hjetmass_box_pppm_0_s34_0_mhsq_s134_s234
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i1, i3)
      t3 = za(i2, i3)
      t4 = zb(i2, i1)
      t5 = zb(i3, i1)
      t6 = zb(i3, i2)
      t7 = za(i1, i4)
      t8 = zb(i4, i1)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t3 * t5 + t8 * t9
      t14 = mt ** 2
      t15 = t2 * t6
      t16 = t7 * t10
      t17 = t15 + t16
      t18 = t4 ** 2
      t19 = t1 ** 2 * t18
      t20 = t19 * t17
      t21 = t3 * t6
      t22 = t9 * t10
      t23 = t22 + t21
      t24 = t1 * t11 * t4 * t12
      t25 = t7 * t8
      t26 = t2 * t5
      t27 = t1 * t4
      t23 = t27 * (t25 * t23 + t26 * t23 - t24)
      t28 = 16 * t14 * t1 * t4 * t13 * t20 + 4 * t23 ** 2
      t28 = sqrt(t28)
      t23 = 2 * t23
      t29 = t23 - t28
      t20 = 0.1q1 / t20
      t30 = t26 + t25
      t30 = t15 * t30 + t16 * t30
      t31 = t22 + t21
      t25 = t25 * t31 + t26 * t31 - t24
      t26 = (0.1q1 / 0.4q1)
      t17 = t26 * t27 * t28 * t20 * t17
      t19 = t19 * t20 * (t21 * t30 + t24 * (-t16 - t15) + t22 * t30) / 2
      t21 = -t17 + t19 - t25
      t13 = t27 * t13
      t24 = 0.1q1 / t1
      t13 = 0.1q1 / t13
      t27 = 0.1q1 / t4
      t30 = t31 * t24
      t32 = t30 * t27 * t7
      t33 = t21 * t13
      t34 = t33 * t9
      t35 = t32 + t34
      t16 = t16 * t26
      t36 = -t16 * t29 * t20 + t8 * t35
      t23 = t23 + t28
      t28 = t5 * t31
      t27 = t27 * t28
      t37 = t26 * t23 * t20
      t38 = -t37 * t1 * t6 + t27
      t17 = t17 + t19 - t25
      t19 = t17 * t13
      t25 = t19 * t9
      t32 = t32 + t25
      t8 = -t16 * t23 * t20 + t8 * t32
      t16 = t26 * t29 * t20
      t1 = -t16 * t1 * t6 + t27
      t16 = -t16 * t7 * t6 + t5 * t35
      t27 = -t37 * t7 * t6 + t5 * t32
      t32 = t30 * t7
      t35 = t34 * t4 + t32
      t32 = t25 * t4 + t32
      t12 = 0.1q1 / t12
      t8 = 0.1q1 / t8
      t37 = 0.1q1 / t11
      t3 = 0.1q1 / t3
      t36 = 0.1q1 / t36
      t39 = 0.1q1 / t9
      t40 = t27 * t8
      t41 = t16 * t36
      t42 = t41 + t40
      t43 = t8 + t36
      t44 = t5 ** 2
      t45 = t1 * t4
      t46 = t4 * t5 * t6
      t47 = t4 * t38
      t48 = t9 * t37
      t49 = t38 * t12
      t50 = t1 * t12
      t51 = t50 * t36
      t52 = t48 * t12
      t53 = t6 * t9
      t54 = t3 * t12
      t55 = t3 * t24
      t56 = t14 - t31
      t57 = t14 * t6
      t58 = t14 * t9 * t24
      t59 = t11 * t44
      t60 = t59 * t12
      t61 = t6 * t24
      t62 = t61 * t12
      t63 = t11 * t36
      t64 = t12 * t37
      t65 = t3 * t9
      t66 = t10 * t27
      t67 = t42 * t31
      t68 = t24 * (-t41 * t54 * t31 * t56 + t40 * t12 * t31 * (-t3 * t56
     # - t6) + t65 * (t64 * (t31 * t32 - t35 * t56) + t57 * (t14 * t12 *
     # t8 + t63))) + t54 * (-t8 * (t57 + t66) * t9 + t67) * t4
      t69 = t51 * t30 * t4 * t16 * t3
      t70 = t47 * t40
      t71 = t70 * t54 * t24 * (-t31 * (t7 * t23 * t39 * t20 + 1) - t14)
      t42 = t5 * t68 + t3 * (t4 * (t12 * (t24 * (t53 * (-t38 - t1) - t41
     # * t14 * t1) + t9 * (t36 * (-t10 * t16 - t57) + t37 * (t35 + t32))
     # * t5) + t58 * t36 * (t16 + t53) + t58 * (t27 + t53) * t8) + t60 *
     # (-t10 * t42 - t57 * t43)) + t12 * (t40 * (-t46 + t22 * (t47 - t28
     #) * t3 * t24) + t41 * ((-t45 * t7 * t29 * t31 * t39 * t20 + t22 * 
     #(t45 - t28)) * t3 * t24 - t46) + t32 * t3 * (-t47 * t48 * t24 + t4
     #4)) + t53 * t14 ** 2 * t5 * t24 * t3 * t36 * t12 + t54 * (-t45 * t
     #48 * t24 + t44) * t35 + t55 * (t5 * (t6 * (t11 * (t8 * (-t49 + t9)
     # - t51) - t9 * t12 * t43 * t31) - t52 * t32 + t22 * t12 * t42) - t
     #52 * (t16 + t27) * t4) * t14 + t62 * (t41 * (t5 * t56 + t45) + t40
     # * (t14 * t5 + t47)) - t60 * t6 * t42 * t39 - t69 + t71 + t67 * t5
     #9 * t54 * t39
      t43 = t31 * t16
      t44 = t31 * t27
      t46 = t40 * t23
      t51 = t41 * t29
      t56 = t4 * t7
      t58 = t16 * t29 + t23 * t27
      t60 = t14 * t7
      t67 = t27 * t39
      t68 = t3 * t23
      t69 = t1 * t36
      t71 = t3 * t7 * t6
      t72 = t20 * t4
      t73 = t30 + t4
      t74 = -t49 + t9
      t75 = t55 * t6
      t76 = t53 * t3
      t77 = t4 * t27
      t15 = t55 * t15
      t78 = t19 * t12
      t30 = t72 * (t23 * (t8 * (t7 * (t12 * (t5 * (t3 * (-t57 - t66) - t
     #67 * t6) + t75 * (t14 - t31) * t38) + t61 * (t11 * t38 * t3 - t27)
     # + t77 * t3 + t76 * t73) + t15 * t27 * t74) + t78 * (-t70 * t3 + (
     #-t65 + t40) * t6 * t5) - t54 * t7 * (t53 * t30 + t77) * t37) + t29
     # * (t7 * (t36 * (t3 * (-t5 * t10 * t12 + t4) - t61) - t54 * t4 * t
     #37) * t16 + t76 * (-t33 * t5 * t12 + t7 * t36 * t73)))
      t66 = t6 * t11
      t70 = t6 * t16
      t53 = t53 * t14
      t76 = t53 * t24
      t15 = t72 * (t29 * (t15 * t41 * t9 + (-t45 * t41 * t3 * t39 + t75 
     #* (-t63 * t1 ** 2 * t39 - t48 * t31 - t69 * t31) + t36 * (t39 * (t
     #3 * (-t66 * t1 + t43) - t70) - t57 * t3) * t5) * t12 * t7 + t33 * 
     #t12 * (t36 * (t3 * (t16 * (t22 - t31) + t53) + t70) * t5 - t16 * t
     #3 * (t69 + t48) * t4)) + t68 * (t7 * t8 * (t76 + (t27 * t28 + t66 
     #* (-t24 * t38 ** 2 - t38 * t5)) * t39 * t12) + t78 * (-t77 * t48 +
     # (t53 - t44) * t8 * t5)))
      t14 = t54 * t6 * (t9 * (-t48 * t14 * t4 * t24 + t5 * t73) + t59)
      t28 = t29 ** 2
      t33 = t23 ** 2
      t9 = t71 * t20 ** 2 * t18 * (t13 * ((t36 * (-t50 + t9) - t52) * t2
     #1 * t28 + t17 * t33 * (t74 * t8 - t52)) + t7 * (t28 * (t36 * (t50 
     #* t39 - 1) + t64) + t33 * (t8 * (t49 * t39 - 1) + t64)))
      ret = -t26 * t9 - 3 * t56 * t62 * t20 * ((-t41 * t39 + t3) * t1 * 
     #t29 + t23 * t38 * (-t40 * t39 + t3)) + 8 * t14 + 4 * t42 + 2 * t56
     # * t55 * t20 * (t31 * (t51 + t46) + t12 * (t51 * t1 + t46 * t38) *
     # t10 + t64 * (t23 * (-t32 * t38 - t44) + t29 * (-t1 * t35 - t43)))
     # + t15 + t30 + t72 * (t12 * (t37 * (t6 * (t7 * t24 * t58 + t65 * (
     #t24 * (-t2 * t58 + t60 * (-t29 - t23)) - t56 * (t29 + t23))) + t3 
     #* (t23 * t32 * (-t25 + t7) + (-t34 + t7) * t35 * t29) * t5) + t68 
     #* (t19 * t22 * t5 * t27 + t47 * (-t67 - t6) * t7) * t8 + t69 * t3 
     #* t29 * t6 * (t24 * (-t16 * t2 + t60) - t56)) + t65 * (t60 * t61 *
     # t29 * t36 + (t46 * t17 + t51 * t21) * t13 * t4) + t71 * (t5 * t23
     # * t8 + (t1 * t24 + t5) * t36 * t29) * t11) - 16 * t76 * t54 * t5

      hjetmass_box_pppm_0_s34_0_mhsq_s134_s234 = ret/32q0/(0,1q0)
      return

      end function
