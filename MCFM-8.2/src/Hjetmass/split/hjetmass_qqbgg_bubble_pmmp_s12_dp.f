
      double complex function hjetmass_qqbgg_bubble_pmmp_s12_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i2, i3)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = 0.1D1 / t4
      t6 = t2 * t3
      t7 = t6 * t5 + t1
      t8 = zb(i4, i1)
      t9 = za(i1, i4)
      t10 = za(i3, i4)
      t11 = 0.1D1 / t2
      t12 = t3 * t5
      t13 = t1 * t11
      t14 = t13 + t12
      t15 = za(i2, i4)
      t16 = 0.1D1 / t9
      t17 = t15 * t16
      t18 = t17 + t12
      t19 = zb(i4, i3)
      t20 = zb(i4, i2)
      t21 = 0.1D1 / t20
      t22 = t21 * t8 + t17
      t23 = t15 * t2
      t24 = t23 * t16
      t25 = t24 - t1
      t26 = t15 * t20
      t27 = t9 * t3
      t28 = t15 * t4
      t29 = t27 + t28
      t30 = t29 * t1
      t31 = t23 * (t28 * t16 + t3)
      t32 = t30 - t31
      t33 = t3 * t20
      t34 = t33 * t5 - t8
      t30 = t30 + t6 * (t27 * t5 + t15)
      t35 = t9 * t20
      t36 = t2 * t4
      t37 = t35 + t36
      t38 = t1 * t3
      t39 = t15 * t8
      t40 = t39 + t38
      t41 = t1 * t4
      t42 = t9 * t8
      t43 = -t41 + t42 - t26 + t6
      t43 = t43 ** 2
      t44 = 4 * t40 * t37 + t43
      t44 = cdsqrt(t44)
      t45 = -t41 + t42 - t26 + t6 + t44
      t46 = 0.1D1 / t37
      t47 = (0.1D1 / 0.2D1)
      t48 = t47 * t45 * t46
      t49 = t48 - t12
      t44 = t41 - t42 + t26 - t6 + t44
      t50 = t47 * t44 * t46
      t51 = -t50 - t12
      t52 = t1 ** 2
      t53 = t3 ** 2
      t54 = t41 * t26
      t55 = t52 * t4 ** 2
      t56 = t2 * t53
      t57 = t9 * t1
      t58 = -t23 + t57
      t59 = -t41 - t6
      t60 = t10 * t19
      t61 = t4 * t8
      t62 = 2
      t63 = -t48 * t20 + t8
      t64 = 0.1D1 / t37
      t65 = t35 * t1
      t66 = t59 * t2
      t67 = t6 * (t42 + t6)
      t68 = t2 * t45 * t64
      t69 = t9 * (t68 - t1) + t23
      t68 = t68 + t1
      t37 = 4 * t40 * t37 + t43
      t37 = cdsqrt(t37)
      t43 = -t41 + t42 - t26 + t6 + t37
      t43 = 0.1D1 / t43
      t70 = t62 * t40
      t71 = -t70 * t43 - t48
      t72 = t46 * (t45 + t44)
      t73 = -t48 * t2 - t1
      t74 = -t41 + t6
      t75 = t64 ** 2
      t76 = t46 * (t60 * (-t23 - t57) + t57 * (t41 - t6) + t23 * t74)
      t77 = t42 * t58
      t78 = t26 * t58
      t79 = t1 * (t15 * (t60 - t6) + t27 * t1)
      t80 = t9 * (-t41 + t60) + t23 * t4
      t81 = t75 * t45 ** 2 * t2
      t82 = (0.1D1 / 0.4D1)
      t83 = t62 * t39 * t58
      t84 = -t47 * t45 * (t35 * (-t23 + t57) * t75 * t45 + t76) + t82 * 
     #t81 * t80 + t79 + (t77 - t78) * t64 * t45 + t83
      t48 = -t48 - t17
      t85 = t2 ** 2
      t86 = t1 * (t15 * (t61 + t33) + t60 * t3)
      t87 = t52 * t4
      t88 = t61 + t33
      t46 = t46 * (-t23 * t88 + t57 * t88 + t60 * t74)
      t74 = t60 * t4 + t88 * t9
      t88 = t62 * t3 * (t2 * (t39 + t38) + t87)
      t53 = -t47 * t45 * (t4 * (t2 * (-t41 - t6) - t65) * t75 * t45 + t4
     #6) + t82 * t81 * t74 - t86 - (t53 * t85 - t41 * (t41 + t26) + t42 
     #* t6) * t64 * t45 - t88
      t81 = t16 ** 2
      t89 = t16 * t81
      t90 = t15 ** 2
      t3 = -t62 * t31 * t1 + t29 * t52 + t85 * (t4 * t15 * t90 * t81 + t
     #3 * t90 * t16)
      t29 = -t50 + t17
      t31 = t5 ** 2
      t91 = t5 * t31
      t92 = t30 * t8 + t33 * (-t1 * t15 + t12 * (-t23 - t57) - t56 * t9 
     #* t31)
      t93 = t50 * t20 + t8
      t94 = t2 * t44 * t64
      t95 = t9 * (t94 + t1) - t23
      t94 = t94 - t1
      t37 = t41 - t42 + t26 - t6 + t37
      t37 = 0.1D1 / t37
      t70 = t70 * t37 + t50
      t50 = t50 * t2 - t1
      t96 = t75 * t44 ** 2 * t2
      t76 = t47 * t44 * (t35 * (t23 - t57) * t75 * t44 + t76) + t82 * t9
     #6 * t80 + t83 + t79 + (-t77 + t78) * t64 * t44
      t46 = t47 * t44 * (t4 * (t2 * (t41 + t6) + t65) * t75 * t44 + t46)
     # + t82 * t96 * t74 - t88 - t86 - (t54 + t55 - t67) * t64 * t44
      t27 = t2 * (t27 - t28) + t41 * t62 * t9
      t28 = t60 + t6
      t42 = t42 + t26
      t47 = -t41 + t6
      t74 = t57 * t47
      t47 = t23 * t47
      t79 = t36 * t58 * t64
      t80 = t35 * t58 * t64
      t82 = 0.1D1 / t51
      t83 = 0.1D1 / t29
      t22 = 0.1D1 / t22
      t10 = 0.1D1 / t10
      t86 = 0.1D1 / t49
      t19 = 0.1D1 / t19
      t88 = 0.1D1 / t18
      t96 = -0.1D1 / t48
      t18 = 0.1D1 / t18
      t14 = 0.1D1 / t14
      t97 = t34 ** 2
      t98 = (t26 * t16 + t8) * t21
      t99 = t98 * t22
      t100 = t7 * t11
      t101 = t34 * t14
      t10 = t10 * t19
      t19 = 0.1D1 / t72
      t70 = 0.1D1 / t70
      t49 = -0.1D1 / t49
      t71 = 0.1D1 / t71
      t51 = -0.1D1 / t51
      t29 = -0.1D1 / t29
      t48 = 0.1D1 / t48
      t72 = t49 ** 2
      t102 = t71 ** 2
      t103 = t51 ** 2
      t104 = t29 ** 2
      t105 = t43 ** 2
      t106 = t63 ** 2
      t107 = t48 ** 2
      t108 = t63 * t73
      t109 = t53 * t72 * t31
      t110 = t73 * t84 * t81 * t107
      t111 = t50 * t93
      t112 = t93 ** 2 * t31 * t103
      t113 = t70 * t37
      t114 = t113 * t44
      t47 = (t114 * (t112 * (t46 * t70 + t62 * (t54 + t55 + t4 * (t66 - 
     #t65) * t64 * t44 - t67) - t61 * t95 - t60 * (t4 * t94 + t6) - t33 
     #* t95) + (t76 * (t2 * t93 + t50 * (-t70 * t93 + t20)) + t111 * (t6
     #2 * (-t80 * t44 - t77 + t78) + t47 - t79 * t44 - t60 * (-t9 * t94 
     #+ t23) - t74)) * t104 * t81) + (t71 * (-t106 * (-t62 * (t54 + t55 
     #+ t4 * (-t66 + t65) * t64 * t45 - t67) - t61 * t69 - t60 * (t4 * t
     #68 - t6) - t33 * t69) * t72 * t31 + (t84 * (t2 * t63 + t20 * t73) 
     #+ t108 * (t62 * (t80 * t45 - t77 + t78) - t60 * (t68 * t9 + t23) +
     # t47 + t79 * t45 - t74)) * t107 * t81) + t63 * (t109 * t63 - t110)
     # * t102) * t43 * t45) * t75 * t19
      t65 = t10 * t64
      t66 = t65 * t19 * t40
      t67 = t19 ** 2
      t68 = t64 * t45
      t44 = t44 * t64
      t43 = t71 * t43
      t69 = t43 * t63
      t71 = t10 * t75
      t44 = t71 * t67 * t40 * (t69 * (t53 * (t68 * t63 * t31 * t49 * t72
     # + (-t68 * t20 + t63) * t31 * t72) - t110 * (t68 * t48 + 1)) + t11
     #3 * t112 * t46 * (t44 * t51 - 1) + t113 * (-t44 * t20 * t46 * t31 
     #* t103 + (-t44 * t29 * t104 + t104) * t81 * t76 * t50) * t93)
      t53 = t18 ** 2
      t68 = t96 ** 2
      t72 = t83 ** 2
      t74 = t32 * t31
      t77 = t25 * t3
      t78 = t77 * t75 * t68 * t72
      t79 = t74 * t53 + t78
      t80 = t82 ** 2
      t84 = t86 ** 2
      t94 = t88 ** 2
      t95 = t92 * t14
      t107 = t97 * t80 * t84 * t75
      t112 = t30 * t81
      t3 = t10 * (t8 * (t98 * t15 * t89 * t79 * t22 ** 2 + t81 * (-t17 *
     # t79 + t98 * ((t17 * (-t2 * t3 - t25 * (t23 * t28 - t57 * t28 + t6
     #2 * (t85 * t90 * t4 * t16 + t23 * t42 - t57 * t42) + t87 * t9 - 3 
     #* t41 * t23)) - t77) * t68 * t75 * t72 - t31 * t53 * (t17 * t27 + 
     #t32))) * t22) + t13 * t14 * t31 * (-t34 * (t34 * t92 * t75 * t84 *
     # t80 + t112 * t94) + t12 * (t107 * (t95 + t62 * (t56 * (t35 * t5 +
     # t2) + t54 + t55) + t26 * t6 - t61 * t58 - t60 * t59 - t38 * (-4 *
     # t36 - 3 * t35)) + (-t27 * t34 + t30 * (t101 - t20)) * t94 * t81))
     #)
      t4 = t66 * t5 * (t43 * t1 * t106 * t49 + t113 * t93 * t51 * (t1 * 
     #t93 + t50 * t8) + t43 * t108 * t8 * t49)
      t6 = t114 * t93
      t9 = t10 * t64 * t75 * t19 * t67 * t40 * (t45 * (t43 * t109 * t106
     # - t69 * t110) + t6 * (-t93 * t46 * t31 * t103 + t76 * t50 * t81 *
     # t104))
      t2 = t1 * t20 + t2 * t8
      t15 = t10 * (t38 * (t101 * t91 * (-t20 * t92 * t80 * t84 * t75 + t
     #112 * t88 * t94) * t11 - t95 * t107 * t91 * (t82 + t86) * t11) + t
     #99 * t39 * t89 * (t74 * t18 * t53 + t78 * (t83 + t96)))
      t20 = t65 * t101 * t38 * t31 * t86 * t82 * (t13 * t20 + t8)
      ret = -8 * t10 * (-t52 * t7 * t97 * t11 ** 2 * t14 ** 2 * t5 * t64
     # * t86 * t82 + t101 * t52 * t11 * t16 * t88 * t5 * (t100 * t14 + 1
     #) + t5 * t16 * t1 * (t22 * t18 * (t98 + t17) + (t100 + t12) * t14 
     #* t88) * t8 + t16 * t21 * t22 * (t98 * t25 ** 2 * t83 * t22 * t64 
     #* t96 + (t25 * (-t99 + 1) + t24) * t18 * t5) * t8 ** 2 + t33 * t52
     # * t11 * t16 * t14 * t88 * t31) + 64 * t66 * (t40 * (t73 ** 2 * t6
     #3 * t16 * t105 * t48 * t102 - t111 * t37 ** 2 * t70 ** 2 * (t50 * 
     #t16 * t29 + t93 * t5 * t51) + t73 * t106 * t5 * t105 * t49 * t102)
     # + t47) + 128 * t44 - 16 * t3 + 12 * t13 * t65 * t101 * t5 * t86 *
     # t82 * (t1 * t34 + t7 * t8) + 48 * t4 - 256 * t9 - 40 * t71 * t5 *
     # t19 * t40 * (t69 * t2 * t49 * t45 - t6 * t51 * t2) - 32 * t15 + 2
     #0 * t20

      hjetmass_qqbgg_bubble_pmmp_s12_dp = -ret/16d0*(0,1d0)
      return

      end function
