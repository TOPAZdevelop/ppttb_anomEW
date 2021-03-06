
      double complex function hjetmass_bubble_pmpm_s34_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double precision mt
          double complex ret

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = zb(i3, i2)
      t4 = za(i1, i2)
      t5 = za(i1, i3)
      t6 = zb(i2, i1)
      t7 = za(i1, i4)
      t8 = zb(i4, i2)
      t9 = 0.1D1 / t8
      t10 = 0.1D1 / t5
      t11 = t3 * t9
      t12 = t7 * t10
      t13 = t11 + t12
      t14 = za(i2, i3)
      t15 = 0.1D1 / t14
      t16 = t1 * t15 + t11
      t17 = zb(i4, i1)
      t18 = 0.1D1 / t17
      t19 = t18 * t2 + t12
      t20 = t7 * t2
      t21 = t1 * t3
      t22 = t20 + t21
      t23 = t5 * t17
      t24 = t14 * t8 + t23
      t25 = t5 * t2
      t26 = t14 * t3
      t27 = t7 * t17
      t28 = t1 * t8
      t29 = -t27 - t28 + t25 + t26
      t29 = t29 ** 2
      t30 = 0.4D1 * t22 * t24 + t29
      t30 = cdsqrt(t30)
      t31 = -t30 - t27 - t28 + t25 + t26
      t32 = 0.1D1 / t24
      t33 = 0.1D1 / 0.2D1
      t34 = t33 * t17
      t35 = -t34 * t31 * t32 + t2
      t36 = 0.1D1 / t24
      t37 = t1 * t22
      t38 = t22 * t36
      t39 = t4 * t6
      t40 = t28 + t27
      t41 = t36 ** 2
      t42 = t41 ** 2
      t43 = t36 * t41
      t44 = t32 * t1
      t45 = t44 * (t39 + t27 + t28)
      t46 = t40 * t41
      t47 = -t39 + t25 + t26
      t48 = t31 ** 2
      t49 = t31 * t48
      t50 = t14 * t48 * t41
      t51 = t14 * t49 * t43 * t24
      t52 = 0.1D1 / 0.4D1
      t53 = -0.1D1 / 0.8D1
      t54 = t38 * t31 * t14 - t33 * t31 * (t46 * t31 * t14 + t45) + t52 
     #* t50 * t47 + t53 * t51 + t37
      t29 = 0.4D1 * t22 * t24 + t29
      t29 = cdsqrt(t29)
      t55 = -t27 - t28 + t25 + t26 - t29
      t56 = t33 * t31 * t32
      t57 = -t56 + t11
      t55 = 0.1D1 / t55
      t58 = 0.2D1 * t22
      t59 = -t58 * t55 - t56
      t30 = t30 - t27 - t28 + t25 + t26
      t60 = t32 * (t31 - t30)
      t61 = -t33 * t14 * t30 * t32 - t1
      t29 = -t27 - t28 + t25 + t26 + t29
      t62 = t33 * t30 * t32
      t63 = -t62 - t12
      t29 = 0.1D1 / t29
      t58 = -t58 * t29 - t62
      t34 = -t34 * t30 * t32 + t2
      t64 = -t56 * t14 - t1
      t56 = -t56 - t12
      t65 = t1 ** 2
      t66 = t30 ** 2
      t67 = t1 * t2
      t68 = t14 * t17
      t69 = t65 * t3
      t70 = t2 * t8
      t71 = t3 * t17
      t72 = t71 + t70
      t73 = t14 ** 2
      t74 = t14 * t2
      t75 = t1 * t17
      t40 = t75 * t40
      t76 = t66 * t41
      t77 = t1 * t72
      t78 = t77 * t41
      t79 = t25 * t17 * t43
      t44 = t44 * t2 * (t39 + t25)
      t80 = t1 * t5
      t81 = t7 * t14
      t82 = t39 * t14
      t83 = t17 * (-t81 + t80) + t82
      t84 = t43 * t30 * t66
      t85 = t26 + t25
      t86 = 0.3D1 / 0.4D1
      t87 = 0.3D1 / 0.2D1
      t32 = t20 * t87 * t32 + 0.2D1 * t21 * t36
      t72 = -t33 * t30 * (t14 * (t78 * t30 + t79 * t66) - t44) - t52 * t
     #76 * (t72 * t36 * t30 * t73 + t40 + t39 * (t75 - t74)) + t53 * t84
     # * t17 * t83 + t86 * t74 * t76 * t85 + t68 * t66 ** 2 * t42 * t24 
     #/ 0.16D2 + t2 * (t7 * (-t68 * t66 * t41 + t67) + t69) + t74 * t30 
     #* t32
      t88 = -t71 - t70
      t89 = t48 * t41
      t32 = t33 * t31 * (t14 * (-t78 * t31 - t79 * t48) + t44) + t52 * t
     #89 * (t88 * t36 * t31 * t73 - t40 + t39 * (-t75 + t74)) + t53 * t4
     #9 * t43 * t17 * t83 + t86 * t50 * t2 * t85 + t2 * (t7 * (-t50 * t1
     #7 + t67) + t69) + t68 * t48 ** 2 * t42 * t24 / 0.16D2 + t74 * t31 
     #* t32
      t40 = t81 * t10
      t42 = t1 - t40
      t44 = t76 * t14
      t49 = -t62 + t11
      t62 = t26 * t9 + t1
      t78 = -t71 * t9 + t2
      t79 = t25 - t39
      t85 = (-t27 - t28 - t39) * t1
      t90 = t79 * t36
      t91 = t73 * t3
      t28 = (-t28 - t27) * t36
      t65 = t65 * t8 * t17
      t92 = t67 * (t39 + t25)
      t93 = (t1 * (-t7 * t17 ** 2 - t39 * t17) - t65 + t82 * t2) * t36
      t70 = t71 + t70
      t71 = t14 * t1
      t94 = t76 * t17
      t95 = t10 ** 2
      t96 = t10 * t95
      t97 = t7 ** 2
      t12 = t12 * t74
      t98 = t27 * t10 + t2
      t99 = t3 ** 2
      t100 = t9 ** 2
      t101 = -0.1D1 / t56
      t102 = -0.1D1 / t63
      t6 = 0.1D1 / t6
      t19 = 0.1D1 / t19
      t4 = 0.1D1 / t4
      t103 = t1 * t98
      t104 = t42 * t2
      t60 = 0.1D1 / t60
      t105 = 0.1D1 / t49
      t56 = 0.1D1 / t56
      t59 = 0.1D1 / t59
      t58 = 0.1D1 / t58
      t106 = 0.1D1 / t57
      t63 = 0.1D1 / t63
      t107 = t60 ** 2
      t108 = t61 * t72
      t109 = t108 * t10 * t63
      t37 = (t38 * t30 * t14 - t33 * t30 * (t46 * t30 * t14 + t45) + t52
     # * t44 * t47 + t53 * t84 * t14 * t24 + t37) * t105
      t38 = t64 * t32 * t10 * t56
      t45 = t55 * t59
      t46 = t4 * t6
      t47 = t46 * t107
      t52 = 0.1D1 / t13
      t16 = 0.1D1 / t16
      t49 = -0.1D1 / t49
      t13 = 0.1D1 / t13
      t53 = -0.1D1 / t57
      t57 = t74 * t18 + t1
      t110 = t78 ** 2
      t111 = t62 * t16
      t52 = t52 * t9
      t112 = t52 * t2
      t16 = t21 * t100 * t16
      t113 = t75 + t74
      t114 = t46 * t10
      t48 = t114 * t36 * (t104 * t97 * t95 * t19 * t102 * t101 * t57 + (
     #t113 * t63 * t29 * t58 * t61 * t66 - t45 * t48 * t64 * t56 * t113)
     # * t60 * t41 * t22)
      t3 = t46 * (t16 * (t110 * (t20 * t1 - t21 * (t27 + t39) * t9 + t79
     # * t14 * t100 * t99 - t68 * t5 * t3 * t99 * t9 * t100 + 0.2D1 * t1
     #1 * t81 * t78) * t41 * t49 ** 2 * t53 ** 2 + (t68 * t99 * t100 - 0
     #.2D1 * t11 * t74 - t67) * t10 * t13) * t15 + t20 * t95 * t19 * t18
     # * (t42 * (-t20 * t39 * t1 * t10 + t69 * t2 + (t39 * (-t75 + t74) 
     #- t65) * t95 * t97 + t39 * t68 * t7 * t97 * t96 + t97 ** 2 * t73 *
     # t17 * t8 * t95 ** 2 - 0.2D1 * t95 * t14 * t97 * (t40 * t88 + t77)
     # + t12 * t3 * (-0.4D1 * t1 + 0.3D1 * t40)) * t41 * t102 ** 2 * t10
     #1 ** 2 + t52 * (t68 * t97 * t95 + 0.2D1 * t12 - t67)))
      t8 = t34 ** 2
      t12 = t35 ** 2
      t39 = t58 ** 2
      t40 = t8 * t9 * t105
      t8 = (t37 * t8 * t9 + t109) * t29
      t65 = t10 * t56
      t66 = t45 * t31
      t23 = t66 * (t65 * (t32 * (t64 * (t56 + t59) - t14) - t64 * (-t93 
     #* t31 - t33 * t51 * t17 - t86 * t89 * t17 * (t17 * (t81 - t80) - t
     #82) + t87 * t89 * t73 * t70 - t92 - 0.4D1 * t74 * (-t27 * t31 * t3
     #6 + t21) - 0.3D1 * t74 * (-t23 * t89 + t20 + (t26 + t25) * t36 * t
     #31) + 0.2D1 * t71 * t31 * t36 * t70)) + (t54 * t106 ** 2 + t106 * 
     #(t54 * t59 + t90 * t31 * t14 - t86 * t50 * t24 + t91 * t31 * t36 +
     # t85 + 0.2D1 * t14 * (t28 * t31 + t20 + t21))) * t9 * t12)
      t5 = t47 * t43 * t22 * (t30 * (t8 * t39 + (t10 * (t63 * (-t14 * t7
     #2 - t61 * (-t93 * t30 - t33 * t84 * t68 * t24 + t86 * t94 * t83 + 
     #t87 * t76 * t73 * t70 + 0.3D1 * t74 * (t94 * t5 - t20 + (-t26 - t2
     #5) * t36 * t30) - t92 + 0.4D1 * t74 * (t27 * t30 * t36 - t21) + 0.
     #2D1 * t71 * t30 * t36 * t70)) + t108 * t63 ** 2) + t40 * (t37 + t9
     #0 * t30 * t14 - t86 * t44 * t24 + t91 * t30 * t36 + t85 + 0.2D1 * 
     #t14 * (t28 * t30 + t20 + t21))) * t29 * t58) + t23)
      t14 = t61 ** 2
      t21 = t114 * t60 * t41 * t22 * (t30 * (t1 * t34 * t58 * t29 * t63 
     #* t61 + t2 * t58 * t29 * t63 * t14) - t66 * t64 * t56 * (t1 * t35 
     #+ t2 * t64))
      t23 = t29 ** 2
      t14 = t46 * t60 * t41 * t22 ** 2 * (t30 * (t34 * t10 * t39 * t23 *
     # t63 * t14 + t40 * t39 * t23 * t61) - t31 * t64 * t35 * t55 ** 2 *
     # t59 ** 2 * (t35 * t9 * t106 + t65 * t64))
      ret = -0.12D2 * t20 * t42 * t4 * t95 * t6 * t19 * t18 * t36 * t102
     # * t101 * (-t104 + t103) - 0.128D3 * t47 * t41 * t22 * ((t37 * t34
     # * t9 * (t17 * t30 * t36 - t34) - t109) * t29 * t58 + t45 * (t35 *
     # t54 * t9 * t106 * (t17 * t31 * t36 - t35) - t38)) - 0.8D1 * t46 *
     # (t7 * (t2 ** 2 * t42 * t98 * t18 ** 2 * t95 * (t42 * t36 * t101 *
     # t102 + t52) * t19 ** 2 + t112 * t18 * t95 * (-t104 + t103) * t19)
     # + t112 * t96 * t57 * t19 * t97 + t16 * (t111 * t1 * t110 * t15 **
     # 2 * t36 * t49 * t53 + (t15 * (-t2 * t62 + (t111 * t15 + 0.1D1) * 
     #t78 * t1) + t11 * (-t75 * t15 - t2)) * t13 * t10)) - 0.20D2 * t48 
     #+ 0.16D2 * t3 + 0.64D2 * t5 - 0.256D3 * t46 * t60 * t107 * t43 * t
     #22 * (-t8 * t30 * t58 + t66 * (t12 * t54 * t9 * t106 + t38)) + 0.2
     #4D2 * t21 + 0.32D2 * t14
      
      hjetmass_bubble_pmpm_s34_dp = ret/16d0*(0,1d0)

      end function
