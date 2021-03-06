
      complex*32 function hjetmass_triangle_ppmm_s23_mhsq_s14
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

          parameter (cg = 1q0)

      t1 = za(i1, i4)
      t2 = zb(i4, i1)
      t3 = za(i1, i2)
      t4 = zb(i2, i1)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t3 * t4
      t14 = t5 * t6
      t15 = t9 * t10
      t16 = t11 * t12
      t17 = t14 + t15 + t16 + t13
      t17 = -4 * t1 * t2 * t7 * t8 + t17 ** 2
      t17 = cg * sqrt(t17) + t13 + t14 + t15 + t16
      t18 = 2 * t7 * t8
      t19 = t18 + t17
      t20 = (0.1q1 / 0.2q1)
      t18 = -t1 * t18 * t2 + t17 ** 2 * t20
      t21 = (0.1q1 / 0.4q1)
      t22 = t17 ** 2
      t23 = t1 * t7
      t24 = t23 * t8 * t2
      t25 = t22 * t21 - t24
      t25 = 0.1q1 / t25
      t26 = t11 * t6 + t4 * t9
      t22 = t22 * t25
      t27 = t22 * t26
      t28 = t17 * t21
      t23 = t23 * t20
      t29 = t11 * t28 + t23 * t4
      t30 = t17 * t8 * t25
      t31 = -t30 * t29
      t32 = t17 * t2 * t25
      t29 = t32 * t29
      t33 = t14 + t13
      t24 = t24 * t20 * t17 * t25
      t34 = t21 * t22 * t33 - t24
      t35 = -t10 * t23 + t28 * t5
      t36 = t30 * t35
      t37 = t32 * t1
      t26 = t37 * t26
      t38 = t10 * t3 + t12 * t5
      t39 = t22 * t38
      t40 = -t30 * (t12 * t23 + t28 * t3)
      t41 = t11 * t2 * t20 * t8 + t28 * t4
      t42 = t17 * t7 * t25
      t43 = -t42 * t41
      t44 = t10 * t11 + t4 * t5
      t45 = t30 * t7
      t46 = t45 * t44
      t47 = t15 + t16
      t48 = t21 * t22 * t47 - t24
      t14 = t14 + t16
      t16 = t22 * t21
      t49 = t16 * t7 * t8
      t50 = -t14 * t20 * t45 + t49
      t16 = t16 * t1 * t2
      t33 = -t20 * t33 * t37 + t16
      t5 = -t2 * t20 * t5 * t8 + t10 * t28
      t10 = t17 * t1 * t25
      t25 = t10 * t5
      t3 = t10 * (t2 * t20 * t3 * t8 + t12 * t28)
      t12 = t37 * t38
      t13 = t15 + t13
      t15 = -t13 * t20 * t45 + t49
      t38 = -t2 * t20 * t8 * t9 + t28 * t6
      t45 = t42 * t38
      t6 = -t23 * t6 + t28 * t9
      t9 = t32 * t6
      t23 = 2 * t1 * t2 + t17
      t14 = t14 * t21 * t22 - t24
      t28 = t22 * t44
      t41 = t10 * t41
      t16 = -t20 * t37 * t47 + t16
      t13 = t13 * t21 * t22 - t24
      t6 = t30 * t6
      t22 = t32 * t35
      t10 = t10 * t38
      t5 = t42 * t5
      t24 = 0.1q1 / t2
      t30 = 0.1q1 / t43
      t32 = 0.1q1 / t8
      t35 = 0.1q1 / t33
      t37 = 0.1q1 / t36
      t17 = 0.1q1 / t17
      t3 = 0.1q1 / t3
      t38 = 0.1q1 / t7
      t42 = 0.1q1 / t1
      t12 = 0.1q1 / t12
      t40 = 0.1q1 / t40
      t44 = 0.1q1 / t9
      t47 = 0.1q1 / t34
      t49 = 0.1q1 / t45
      t39 = 0.1q1 / t39
      t51 = 0.1q1 / t29
      t52 = mt ** 2
      t18 = 0.1q1 / t18 ** 2
      t53 = t23 ** 2
      t54 = t52 * t27
      t55 = t54 * t32 * t17 * t38
      t7 = t7 * t8
      t8 = t7 * t26
      t56 = t8 * t53 * t18
      t57 = t11 ** 2 * t38 * t42
      t58 = t4 ** 2 * t24 * t32
      t59 = t57 + t58
      t60 = t50 ** 2
      t61 = t25 ** 2
      t62 = t43 ** 2
      t63 = t44 ** 2
      t64 = t44 * t63
      t65 = t33 ** 2
      t66 = t49 ** 2
      t67 = t49 * t66
      t68 = t35 ** 2
      t69 = t34 ** 2
      t70 = t40 ** 2
      t71 = t40 * t70
      t72 = t29 ** 2
      t73 = t3 ** 2
      t74 = t3 * t73
      t53 = t7 * t53
      t75 = t53 * t41
      t1 = t1 * t2
      t2 = t1 * t19 ** 2
      t76 = t2 * t31
      t77 = t23 * t41
      t78 = t77 * t19
      t79 = t41 * (-t55 - t56)
      t80 = t26 * t18
      t81 = t80 * t72
      t55 = t55 * t31
      t82 = t75 * t33 * t28 * t24 * t18 * t42
      t83 = t80 * t19
      t84 = t46 * t33
      t85 = t43 * t25
      t1 = t1 * t19
      t86 = t29 * t23
      t87 = t34 * t46
      t88 = t33 * t19
      t89 = t27 * t14
      t90 = t1 * t31 * t32 * t38
      t91 = t33 * t25
      t92 = t18 * t19
      t93 = t92 * t27
      t94 = t2 * t18
      t95 = t72 * t35
      t96 = t47 * t62
      t97 = t44 * t95 + t49 * t96
      t98 = t26 ** 2
      t99 = t14 ** 2
      t100 = t15 * t41
      t101 = t43 * t28
      t102 = t29 * t12
      t103 = t102 * t3
      t104 = t36 * t26
      t4 = t11 * t4
      t11 = t52 * t29
      t105 = t92 * (t27 * (-t41 * t97 + t61 * (t102 * t30 * t33 * t50 * 
     #t73 - t103)) + t34 * (t26 * (-t95 * t37 * t63 * (t13 * t31 + t100)
     # - t101 * t36 * t40 * t39 * t51) + t95 * t15 * t13 * t37 * t64 * t
     #98) + t104 * t69 * t50 * t28 * t70 * t39 * t51) * t23
      t106 = t53 * t18
      t107 = t106 * t27 * t65 * t61 * t99 * t42 * t24 * t30 * t74 * t12
      t55 = t35 * (t81 * (-t78 * t45 + (t13 * t24 * t42 * t75 + t15 * t3
     #2 * t38 * t76) * t37 * t34) * t63 + t29 * (-t18 * t26 * t32 * t38 
     #* t43 * t76 + t24 * t29 * t42 * t79) * t44) - t91 * t23 * t18 * (t
     #24 * t25 * t30 * t42 * t7 * t86 * t89 + t28 * t6 * t88) * t12 * t7
     #3 + (t33 * (-t23 * t25 * t46 * t83 + t25 * t59 - t82) + t2 * t25 *
     # t32 * t18 * t38 * (-t84 + t85) * t27) * t12 * t3 + (-t56 * t69 * 
     #t36 * t14 * t28 * t70 * t39 * t51 + (t55 * t66 * t9 + t49 * t79) *
     # t47 * t62 - t55 * t29 * t43 * t47 * t49) * t42 * t24 + t81 * t19 
     #* t34 * (-t77 + t90) * t44 * t68 + t93 * (t1 * t32 * t36 * t38 * t
     #39 * t40 * t87 - t31 * t43 * t47 * t49 * t86) + t94 * t27 * t65 * 
     #t61 * t60 * t38 * t32 * t30 * t74 * t12 + t107 + t11 * t26 * t42 *
     # t38 * t32 * t24 * t35 * t44 * (-t104 * t17 * t29 * t35 + t4) + t1
     #05
      t79 = 0.1q1 / t25
      t81 = t27 ** 2
      t105 = t47 ** 2
      t107 = t31 * t34 * t28
      t108 = t27 * t26
      t109 = t81 * t33
      t110 = t108 * t13
      t111 = t84 * t27
      t112 = t27 * t50
      t113 = t27 * t36
      t114 = t113 * t29
      t115 = t43 * t46
      t116 = t87 * t29
      t117 = t116 * t44
      t118 = t34 * t41
      t119 = t108 * t25 * t49
      t120 = t106 * t42 * t24
      t121 = t94 * t38 * t32
      t103 = t52 * t81 * t62 * t25 * t17 * t42 * t38 * t32 * t24 * t105 
     #* t49 + t120 * (t28 * (t27 * (t12 * t14 * t25 * t30 * t65 * t73 - 
     #t103 * t30 * t91) + t118 * t29 * t26 * t35 * t37 * t44) + t62 * (t
     #110 * t47 * t66 + t105 * t119) - t119 * t29 * t43 * t35 * t47) + t
     #121 * (t26 * (t35 * (-t114 * t43 * t44 * t47 + t15 * t27 * t63 * t
     #72) + t113 * t72 * t44 * t68) + t98 * (t35 * (t116 * t15 * t37 * t
     #63 - t115 * t44) + t117 * t68) + t111 * t31 * t43 * t79 * t47 * t4
     #9) - t59 * t44 * t35 * t26 * t29 + t92 * (t29 * (-t13 * t35 * t37 
     #* t63 * t87 * t98 - t107 * t26 * t35 * t37 * t44) - t111 * t43 * t
     #41 * t79 * t47 * t49 + t96 * (-t109 * t13 * t15 * t67 * t79 - t108
     # * t15 * t66) - t112 * t65 * t25 * t28 * t30 * t73 * t12 - t110 * 
     #t95 * t63) * t23
      t108 = t86 * t7 * t42 * t24
      t110 = t19 * t43
      t111 = -t110 + t108
      t116 = t34 * t13
      t119 = t116 * t37 * t44 - 1
      t122 = t15 ** 2
      t123 = t26 * t46
      t124 = t23 * t13
      t125 = t19 * t15
      t126 = t43 * t18
      t127 = t27 * t43
      t128 = t29 * t46
      t57 = t105 * (t85 * t23 * t18 * t81 * t111 * t49 + t109 * t23 * t1
     #8 * t62 * (-t124 * t24 * t42 * t7 + t125) * t66) + t47 * (t127 * (
     #t82 * t79 - t57 - t58 + t92 * (-t28 * t31 * t33 * t79 + t123) * t2
     #3) * t49 - t109 * t62 * t79 * t18 * (t13 ** 2 * t24 * t42 * t53 + 
     #t122 * t2 * t32 * t38) * t67 + t126 * t81 * (t124 * t111 + t125 * 
     #(t1 * t32 * t38 * t43 - t86)) * t66) - t120 * t91 * t81 * t62 * t4
     #7 * t105 * t49 + t80 * t44 * t35 * (t26 * (t53 * t29 * t24 * t42 *
     # t119 * t28 + t19 * (-t29 * t34 * (t15 * t37 * t44 + t35) + t43) *
     # t28 * t23) - t128 * t2 * t27 * t38 * t32)
      t58 = t87 * t26
      t67 = t51 * t58 + t113
      t82 = t118 + t104
      t109 = t12 ** 2
      t113 = t12 * t109
      t129 = t94 * t46
      t130 = t52 * t28 * t24 * t17 * t42
      t131 = t31 * t32 * t38
      t132 = t26 * t50
      t133 = t23 * t14
      t134 = t92 * t36
      t135 = t134 * t34
      t88 = t88 * t18
      t136 = t14 * t3 * t33 - t29
      t137 = t65 * t28
      t138 = t26 * t43
      t139 = t136 * t25
      t140 = t42 * t38
      t141 = t140 * t17 * t32 * t24
      t142 = t141 * t52
      t6 = t26 * (t131 * t17 * t24 * t42 * t52 * t97 + t35 * t41 * t43 *
     # t44 * t86 * t92) + t65 * (t3 * (-t92 * t48 * t25 * t23 * t28 * t1
     #09 + t131 * (t130 + t129) * t12) - t92 * t89 * t23 * t61 * t50 * t
     #30 * t74 * t12) + t135 * (t1 * t32 * t38 * t50 * t67 - t133 * t67)
     # * t39 * t70 + t88 * t61 * (t23 * (-t132 + t89) - t112 * t1 * t32 
     #* t38) * t12 * t73 + (t36 * (-t121 * t43 * t67 + t36 * t86 * t93 -
     # t34 * t59) + t120 * t34 * t82 * t28) * t39 * t40 + t142 * (t29 * 
     #t26 * t41 * t35 * t44 * (-t29 * t44 * t45 + t43) - t137 * t6 * t25
     # * t73 * t12) + t120 * t12 * t3 * t25 * (t12 * t137 * t16 + t139 *
     # t26) + t92 * (t31 * t26 * t97 + (t137 * t31 + t138 * t61) * t12 *
     # t3) * t23
      t45 = 0.1q1 / t26
      t59 = 0.1q1 / t27
      t67 = t1 * t48 * t32 * t38
      t74 = t23 * t16
      t89 = t74 - t67
      t97 = t74 * t7 * t24 * t42
      t112 = t19 * t48
      t121 = -t112 + t97
      t137 = t39 ** 2
      t143 = t36 ** 2
      t131 = t131 * t52 * t17
      t144 = t19 * t89
      t145 = t1 * t50 * t32 * t38
      t146 = t133 - t145
      t147 = t34 * t28
      t82 = t82 * t40
      t148 = t52 * t32
      t149 = t148 * t17 * t38
      t150 = t133 * t7 * t24 * t42
      t151 = t19 * t50
      t152 = t39 * t40
      t153 = t152 * t36
      t154 = t153 * t34
      t155 = t31 * t14
      t156 = t91 * t3
      t157 = t156 * t12
      t14 = t33 * (t120 * t29 * t61 * t16 * t3 * t109 + t85 * (-t76 * t3
     #2 * t18 * t38 + (t56 * t59 + t149) * t42 * t24 * t41) * t12 * t3) 
     #+ t65 * (t18 * (t112 * t146 + t74 * (-t150 + t151)) * t109 * t73 *
     # t61 - t78 * t18 * t25 * t5 * t3 * t109) + t154 * ((t106 * (t147 *
     # t16 * t39 - t14 * t82) - t149 * t87 * t10 * t40) * t42 * t24 + t9
     #2 * (t23 * (-t147 * t39 * t48 + t50 * t82) + t90 * (-t34 * t40 * t
     #50 + t43))) + t152 * t148 * t140 * t34 * t24 * (t17 * t41 * t87 - 
     #t36 * t4) + t92 * (t69 * (-t132 * t14 * t143 * t39 * t51 * t71 + t
     #155 * t36 * t39 * t70 + t152 * t41 * t46) - t153 * t118 * t43 - t1
     #57 * (t102 * t25 * t48 + t138 * t31 * t59)) * t23
      t56 = t34 * t98
      t78 = t98 * t43
      t82 = t157 * t28
      t90 = t33 * t49 * t105
      t102 = t33 * t13 * t79
      t106 = t23 * t18
      t9 = t92 * (t27 * (t157 * t46 * t30 * t136 + ((-t102 - t9) * t66 *
     # t47 - t90) * t62 * t31) + t98 * (t35 * (t72 * t15 * t63 + t128 * 
     #t44) + t72 * t36 * t44 * t68) + t117 * t26 * t41 * t35 * t37) * t2
     #3
      t9 = t106 * (t19 * (-t27 * (t100 * t96 * t33 * t66 * t79 + t153 * 
     #t147) + t78 * t29 * t13 * t35 * t63 - t116 * t72 * t98 * t68 * t63
     #) + t7 * (t26 * (t95 * t26 * t13 * t63 * t119 + t82) + t27 * (-t29
     # * t41 * t47 * t49 * t43 + (t102 * t66 * t47 + t90) * t41 * t62)) 
     #* t42 * t24 * t23) + (t94 * (t27 * (-t46 * t50 * t65 * t25 * t30 *
     # t73 * t12 + t96 * (t15 * t33 * t66 * t79 - t49) * t31) + t29 * (t
     #35 * (-t31 * t37 * t44 * t58 - t15 * t63 * t78) - t78 * t36 * t44 
     #* t68) + t72 * (t35 * (t122 * t37 * t56 * t64 - t27 * t31 * t44) +
     # t56 * t15 * t63 * t68 + t56 * t36 * t44 * t35 * t68)) + t54 * t42
     # * t24 * (t4 * t43 * t47 * t49 + t82 * t17)) * t38 * t32 + t141 * 
     #t11 * t98 * t46 * t35 * t44 + t9
      t11 = t88 * t115
      t1 = t27 * (t47 * (t49 * (-t127 * (t130 + t129) * t38 * t32 + t106
     # * (-t101 * t8 * t23 * t24 * t42 + (t128 - t101) * t27 * t19)) + t
     #11 * t27 * t79 * (t1 * t15 * t32 * t38 - t124) * t66) - t11 * t23 
     #* t27 * t105 * t49 + t83 * t86 * t28 * t35 * t44)
      t8 = t53 * t16 ** 2 * t42 * t24 + t2 * t48 ** 2 * t38 * t32
      t11 = t19 * t31
      t15 = t18 * (t69 * (t11 * t23 * t22 * t40 * t137 * t36 + (t112 * (
     #-t133 + t145) + t74 * (-t151 + t150)) * t137 * t70 * t143) + t110 
     #* (t74 - t67) * t137 * t40 * t143 * t34 + t85 * t109 * t3 * t59 * 
     #t65 * (t77 * (-t112 + t97) + t11 * (-t74 + t67)))
      t2 = t135 * (-t104 * t146 * t51 * t43 + t87 * t23 * t10) * t39 * t
     #70 + t152 * (t106 * (-t118 * t114 * t19 * t45 + t11 * t69 * t28 - 
     #t143 * t111 * t26) + t34 * (t76 * t18 * (t114 * t45 + t87) + (-t10
     #4 * t46 + t107) * t42 * t17 * t24 * t52) * t38 * t32) - t80 * t69 
     #* t143 * t51 * (t2 * t60 * t32 * t38 + t53 * t99 * t24 * t42) * t3
     #9 * t71 + t12 * t3 * t33 * ((t94 * t91 * (t46 * t48 * t12 - t31 * 
     #t50 * t3) + t52 * t24 * t42 * (t84 * t41 * t17 - t4 * t25)) * t38 
     #* t32 + t106 * (t19 * (-t31 * t29 * t25 + t84 * (-t25 * t16 * t12 
     #+ t41) + t156 * (t41 * t50 + t155)) - t139 * t77 * t7 * t42 * t24)
     #)
      t4 = -128
      t10 = 512
      t11 = t92 * t23 * (-1024 * t65 * t43 * t61 * t48 * t16 * t59 * t3 
     #* t113 + t157 * (-6 * t27 * t28 + 96 * t43 * t41) + t154 * (t29 * 
     #(-4096 * t34 * t36 * t48 * t16 * t45 * t137 - 192 * t31) + 48 * t1
     #23))
      ret = -t20 * t1 + t21 * t106 * t49 * t47 * t28 * t81 * ((-t125 * t
     #79 * t49 + t7 * (t13 * t79 * t49 + t47) * t42 * t24 * t23) * t33 *
     # t43 - t108) + 8 * t55 + 16 * t6 + t4 * (t34 * (t24 * t29 * t42 * 
     #(-t18 * t75 + t131) * t39 * t40 * t36 - t93 * t29 * t45 * t89 * t1
     #37 * t40 * t143) + t91 * t109 * t3 * (t142 * t33 * t41 * t5 + t85 
     #* (t121 * t23 * t26 * t59 + t144) * t18) - t134 * t46 * t89 * t137
     # * t40 * t69) + 64 * t14 + t10 * ((-t131 * t22 * t24 * t42 + (t77 
     #* t121 - t144 * t31) * t18 * t45 * t29) * t137 * t40 * t36 * t69 +
     # t126 * t65 * t61 * t59 * t3 * t113 * t8 + t86 * t18 * t121 * t137
     # * t40 * t143 * t34) - 256 * t15 - 32 * t2 + 2048 * t18 * t29 * t6
     #9 * t143 * t45 * t39 * t137 * t40 * t8 - 2 * t103 + t57 - 4 * t9 +
     # t11

      hjetmass_triangle_ppmm_s23_mhsq_s14 = ret/32q0/(0,1q0)
      return

      end function
