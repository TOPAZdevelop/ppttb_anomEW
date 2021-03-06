
      complex*32 function hjetmass_box_pppm_0_0_s12_mhsq_s34_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i2, i1)
      t3 = zb(i3, i2)
      t4 = za(i1, i2)
      t5 = za(i1, i3)
      t6 = zb(i4, i1)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = za(i1, i4)
      t10 = zb(i4, i2)
      t11 = zb(i3, i1)
      t12 = t9 * t11
      t13 = t1 * t3
      t14 = t12 + t13
      t15 = mt ** 2
      t16 = za(i2, i3)
      t17 = t5 * t6
      t18 = t16 * t10
      t19 = t18 + t17
      t20 = t7 * t8
      t21 = t20 * t19
      t22 = t4 * t2
      t23 = t9 * t6
      t24 = t1 * t10
      t25 = t20 * (t24 + t22 + t23)
      t26 = 16 * t15 * t14 * t21 + 4 * t25 ** 2
      t26 = sqrt(t26)
      t25 = 2 * t25
      t27 = -t25 - t26
      t21 = 0.1q1 / t21
      t28 = t24 + t23
      t29 = t20 * t21
      t28 = t29 * (t17 * t28 + t22 * (t18 + t17) + t18 * t28)
      t30 = (0.1q1 / 0.2q1)
      t19 = t30 * t26 * t21 * t19
      t22 = 2 * t24 + 2 * t22 + 2 * t23
      t31 = t22 - t19 - t28
      t14 = 0.1q1 / t14
      t32 = (0.1q1 / 0.4q1)
      t33 = t11 * t31 * t14
      t34 = t33 * t30
      t17 = t17 * t32
      t35 = t17 * t27 * t21
      t36 = t9 * (t6 - t34) + t35
      t25 = -t25 + t26
      t19 = t22 + t19 - t28
      t22 = t11 * t19 * t14 * t30
      t17 = t17 * t25 * t21
      t26 = t9 * (t6 - t22) + t17
      t28 = t32 * t27 * t21
      t37 = t31 * t30 * t14 * t9 * t3 - t28 * t5 * t10
      t35 = t34 * t9 - t35
      t38 = t32 * t25 * t21
      t39 = t19 * t30 * t14 * t9 * t3 - t38 * t5 * t10
      t17 = t22 * t9 - t17
      t28 = -t28 * t16 * t6 + t34 * t1
      t22 = -t38 * t16 * t6 + t22 * t1
      t34 = t9 * t16
      t38 = -t1 * t5 + t34
      t40 = t27 * t31
      t41 = t25 * t19
      t4 = 0.1q1 / t4
      t42 = 0.1q1 / t6
      t43 = 0.1q1 / t5
      t44 = t9 ** 2
      t45 = t11 * t4
      t46 = t15 * t4
      t47 = t4 ** 2
      t48 = t2 * t8
      t49 = t48 * t1
      t50 = t49 * t9
      t51 = t13 * t9 * t43
      t52 = t47 * t15
      t36 = 0.1q1 / t36
      t26 = 0.1q1 / t26
      t53 = 0.1q1 / t9
      t54 = t25 * t39
      t55 = t27 * t37 + t54
      t56 = t25 * t26
      t57 = t36 * t27
      t58 = t57 + t56
      t59 = t27 * t35
      t60 = t17 * t25
      t61 = t60 + t59
      t62 = t27 + t25
      t63 = t43 * t4
      t64 = t24 * t53
      t65 = t3 * t6
      t66 = t2 * t43
      t67 = t21 * t7
      t68 = t3 * t35
      t69 = t37 * t43 + t3
      t70 = t15 + t35
      t71 = t4 * t8
      t72 = t39 * t43
      t73 = t3 * t17
      t74 = t3 * t11
      t75 = (t11 * t37 - t68) * t43
      t76 = t71 * t1
      t34 = t21 * (t7 * (t2 * (t57 * (t63 * t70 * t9 * t8 + t74 + t75) +
     # t56 * (t11 * (t72 + t3) + t43 * (t71 * t9 - t3) * t15 - t73 * t43
     #)) + t4 * (t63 * t55 * t8 * t1 - t74 * t56 * t15) + t45 * t15 * t9
     # * t43 * t58 * t10 + t76 * (t25 * (-t66 * t39 - t73 * t4) - t2 * t
     #69 * t27) * t53 * t42) + t34 * t48 * t52 * t42 * t43 * t62) + t67 
     #* (-t57 * t15 * t3 * (t66 + t45) + t4 * (t65 * t15 * t43 * t58 + t
     #48 * (-t57 * t11 + t56 * (t17 * t43 - t11))) * t9 + t4 * (t8 * (-t
     #13 * (t2 * t25 + t59 * t4) * t53 + t63 * (t25 * (-t17 ** 2 + t22 *
     # t39) + t27 * (t28 * t37 - t35 ** 2) + t24 * t61)) + t64 * (t3 * t
     #62 + t43 * t55) * t11 - t13 * t10 * t43 * t53 * t61) * t42)
      t55 = t41 + t40
      t58 = t25 ** 2
      t77 = t27 ** 2
      t78 = t77 * t35
      t79 = t58 * t17
      t80 = t79 + t78
      t81 = t77 * t36
      t82 = t58 * t26
      t83 = t82 + t81
      t84 = t24 * t11
      t79 = t79 * t26
      t85 = t3 * t7
      t86 = t71 * t21
      t73 = t43 * (t11 * t39 - t73) + t74
      t87 = t77 + t58
      t88 = t10 * t11
      t18 = t67 * (t21 * (t58 * (-t8 * t17 * t47 + (t53 * (t10 * t73 + t
     #48 * (-t72 - t3)) + t71 * (t48 - t88 + t65)) * t26 * t7) + t77 * (
     #-t8 * t35 * t47 + (t8 * (t4 * (t48 + t65) - t2 * t3 * t53) + t88 *
     # (t53 * t69 - t71)) * t36 * t7) + t64 * t8 * t11 * t47 * t87 * t42
     # * t5) + t51 * t8 * t47 * t55 * t14) + t67 * (t43 * (t8 * (t14 * (
     #t41 * (t9 * (t22 * t3 + t84) - t13 * t17) + t40 * (t9 * (t28 * t3 
     #+ t84) - t13 * t35)) - t50 * t14 * t55 + t18 * t80 * t21) * t42 * 
     #t47 - t81 * t67 * t53 * (t68 * t10 + t48 * t37)) + t48 * t21 * t4 
     #* (t11 * t83 + t76 * (-t77 - t58) * t53 * t42) * t5 + t86 * (t23 *
     # t83 * t2 - t85 * t53 * (t78 * t36 + t79)))
      t23 = t11 * t62
      t47 = t9 * t62
      t55 = t21 * t4
      t64 = t27 * t28
      t6 = t55 * (t7 * (t25 * t73 + t27 * (t75 + t74) + t11 ** 2 * t62 *
     # t42 * t10) + t48 * (t64 * t9 * t4 * t42 + t79 * t67)) + t86 * (t7
     # * (t42 * (t4 * (t11 * (-t60 - t59) + t3 * (t22 * t25 + t64)) - t2
     #3 * t2) + t2 * (t82 * t15 + t81 * t70) * t21 - t66 * t57 * t33 * t
     #44 + t49 * t4 * t42 * t62 - t47 * t66) - t63 * t2 * t44 * t62 * t1
     #6 - t72 * t82 * t21 * t7 ** 2 * t6 + t4 * t2 * t42 * (t25 * (-t1 *
     # t17 + t22 * t9) - t59 * t1)) + t55 * (t76 * (t47 * t2 + t85 * t62
     # + t7 * (t80 * t21 * t53 - t23) * t42 * t10) + t7 * (t29 * (t79 * 
     #(t39 * t53 - t10) + t81 * (t35 * (t37 * t53 - t10) - t37 * t6)) + 
     #t14 * t11 * t44 * t8 * (-t40 * t4 + t41 * (-t2 * t26 - t4)) + (t48
     # * t62 - t88 * t62) * t42 * t15) * t43)
      t23 = t48 - t88 - t65
      t26 = t25 * t58 * t26
      t29 = t27 * t77 * t36
      t20 = t20 * t21 ** 2
      t33 = t29 + t26
      t8 = t25 * t41 * t21 * t14 * t8 * t38 + t27 * t40 * t21 * t14 * t8
     # * t38
      t25 = t46 * t42
      ret = -t30 * t6 + t32 * t18 + 8 * t46 * t43 * (-t13 * t2 * t42 + t
     #45 * t44) + t20 * t4 * (t45 * t5 * t87 + t67 * (t26 * t17 + t29 * 
     #t35) * t53 * t10 + t7 * (t82 * t23 * t19 + t81 * t31 * t23) * t43 
     #* t14 * t9 + t4 * (t48 * t87 - t88 * t87) * t42 * t16) / 8 + t67 *
     # t4 * (t63 * t42 * (t48 * t8 - t88 * t8) + t20 * (t48 * t33 - t88 
     #* t33) * t53 * t5) / 16 - 4 * t52 * (t42 * (t11 * (t43 * (t1 * (t3
     #7 + t39) + t9 * (-t24 + t35 + t17)) + t12 + t13) - t50 * t43) - t5
     #1) - t34 + 2 * t4 * (t43 * (t67 * (t15 * t11 * (t57 * t37 + t56 * 
     #t39) + t71 * (-t61 * t9 + (t54 * t17 + t59 * t37) * t53 * t42 * t1
     #)) + t25 * t3 * (t1 * (-t35 - t17) + t9 * (t28 + t22))) + t25 * t4
     #9 * t21 * t62)

      hjetmass_box_pppm_0_0_s12_mhsq_s34_s124 = ret/32q0/(0,1q0)
      return

      end function
