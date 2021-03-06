
      complex*32 function hjetmass_triangle_pmpm_s23_mhsq_s14
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

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = za(i3, i4)
      t10 = zb(i4, i3)
      t11 = za(i2, i3)
      t12 = zb(i3, i2)
      t13 = t1 * t2
      t14 = t3 * t4
      t15 = t7 * t8
      t16 = t9 * t10
      t17 = t13 + t14 + t15 + t16
      t17 = -4 * t5 * t11 * t12 * t6 + t17 ** 2
      t17 = cg * sqrt(t17) + t13 + t14 + t15 + t16
      t18 = 2 * t5 * t6
      t19 = t18 + t17
      t20 = (0.1q1 / 0.2q1)
      t18 = -t18 * t11 * t12 + t20 * t17 ** 2
      t21 = 2 * t11 * t12 + t17
      t22 = (0.1q1 / 0.4q1)
      t23 = t17 ** 2
      t24 = t5 * t11
      t25 = t24 * t12 * t6
      t26 = t23 * t22 - t25
      t26 = 0.1q1 / t26
      t27 = t2 * t7 + t4 * t9
      t23 = t23 * t26
      t28 = t23 * t27
      t29 = t13 + t14
      t30 = t17 * t5
      t31 = t30 * t6 * t26
      t32 = t23 * t22
      t33 = t32 * t5 * t6
      t34 = -t20 * t31 * t29 + t33
      t35 = t17 * t22
      t24 = t24 * t20
      t36 = -t24 * t4 + t35 * t7
      t37 = t17 * t6 * t26
      t38 = t37 * t36
      t39 = -t20 * t7 * t12 * t6 + t35 * t4
      t40 = t17 * t11 * t26
      t41 = t40 * t39
      t14 = t14 + t16
      t42 = t40 * t12
      t32 = t32 * t11 * t12
      t43 = -t20 * t42 * t14 + t32
      t44 = t20 * t9 * t12 * t6 + t35 * t2
      t45 = -t40 * t44
      t46 = t1 * t8 + t10 * t3
      t47 = t31 * t46
      t30 = t30 * t26
      t39 = t30 * t39
      t48 = t20 * t1 * t12 * t6 + t35 * t10
      t49 = t30 * t48
      t50 = t17 * t12 * t26
      t36 = t50 * t36
      t25 = t25 * t20 * t17 * t26
      t26 = t22 * t23 * t29 - t25
      t16 = t15 + t16
      t29 = -t20 * t31 * t16 + t33
      t2 = t24 * t2 + t35 * t9
      t9 = t37 * t2
      t27 = t31 * t27
      t31 = t23 * t46
      t33 = t35 * t1 + t24 * t10
      t46 = -t50 * t33
      t1 = t1 * t4 + t10 * t7
      t10 = t23 * t1
      t14 = t22 * t23 * t14 - t25
      t51 = t30 * (-t20 * t3 * t12 * t6 + t35 * t8)
      t13 = t13 + t15
      t15 = t22 * t23 * t13 - t25
      t13 = -t20 * t42 * t13 + t32
      t16 = t22 * t23 * t16 - t25
      t2 = -t50 * t2
      t3 = t50 * (-t24 * t8 + t35 * t3)
      t1 = t42 * t1
      t8 = -t40 * t48
      t23 = t37 * t33
      t24 = 0.1q1 / t26
      t25 = 0.1q1 / t34
      t32 = 0.1q1 / t11
      t33 = 0.1q1 / t5
      t35 = 0.1q1 / t49
      t37 = 0.1q1 / t6
      t3 = 0.1q1 / t3
      t40 = 0.1q1 / t45
      t42 = 0.1q1 / t12
      t47 = 0.1q1 / t47
      t48 = 0.1q1 / t51
      t17 = 0.1q1 / t17
      t50 = 0.1q1 / t38
      t51 = mt ** 2
      t18 = 0.1q1 / t18 ** 2
      t11 = t11 * t12
      t12 = t11 * t19 ** 2
      t52 = t12 * t27
      t53 = t52 * t18
      t54 = t28 * t51 * t32 * t42 * t17 + t53
      t55 = t19 * t21
      t56 = t55 * t18
      t57 = t51 * t32
      t58 = t57 * t33 * t37 * t42
      t59 = t58 * t17
      t60 = t56 + t59
      t61 = t19 * t15
      t62 = t21 * t13
      t63 = t61 * t11 * t33 * t37 - t62
      t64 = t21 * t36
      t65 = t64 * t5 * t6 * t32 * t42
      t66 = t19 * t39
      t67 = t66 - t65
      t68 = t34 * t41
      t69 = -t68 * t14 * t35 * t40 + t38
      t70 = t25 ** 2
      t71 = t27 ** 2
      t72 = t21 ** 2
      t73 = t46 ** 2
      t74 = t24 ** 2
      t75 = t3 ** 2
      t76 = t3 * t75
      t77 = t48 ** 2
      t78 = t48 * t77
      t79 = t15 ** 2
      t80 = t12 * t18
      t81 = t80 * t38
      t82 = t81 * t39
      t83 = t37 * t33
      t84 = t5 * t6
      t85 = t84 * t72
      t86 = t85 * t41
      t87 = t56 * t38
      t88 = t18 * t34
      t89 = t88 * t28 * t21
      t90 = t21 * t41
      t91 = t19 * t38
      t92 = t56 * t71 * t26 * t15
      t93 = t85 * t18
      t94 = t93 * t38
      t95 = t83 * t41
      t96 = t95 * t10
      t97 = t56 * t28
      t98 = t71 * t1
      t99 = t47 * t40
      t100 = t99 * t34
      t101 = t13 ** 2
      t102 = t28 * t38
      t103 = t102 * t46
      t104 = t26 * t1
      t105 = t104 * t27
      t106 = t84 * t26 * t71 * t73
      t30 = t104 * t30 * t44
      t44 = t41 * t46
      t107 = t19 * t27
      t108 = t48 * t15
      t109 = t34 * t24
      t72 = t72 * (t25 * (t84 * t41 * t71 * t73 * t13 * t32 * t42 * t50 
     #* t75 - t106 * t101 * t32 * t42 * t50 * t76 + t84 * t36 * t32 * t4
     #2 * (t105 + t103) * t3) - t106 * t13 * t32 * t75 * t42 * t70) - t1
     #09 * t83 * t12 * t28 * t39 * t49 * t48 * (t24 * t41 + t108) + (-t9
     #1 * t28 * t39 * t46 * t3 + t107 * t46 * (-t44 * t27 * t15 * t50 + 
     #t30) * t75) * t25 * t21
      t5 = t18 * t72 + t49 * (t24 * (t89 * (t13 * t67 + t61 * t36) * t77
     # + t28 * (-t87 * t36 + t86 * t36 * t32 * t42 * t18 + t83 * (-t57 *
     # t7 * t4 * t42 + t82)) * t48) + t68 * t56 * t28 * t36 * t74 * t48)
     # + t73 * (t25 * (-t80 * t71 * t26 * t79 * t33 * t37 * t50 * t76 + 
     #t19 * t71 * t18 * t63 * t75) + t70 * (t21 * t71 * t18 * (t90 * t5 
     #* t6 * t32 * t42 - t91) * t3 + t92 * t75) - t94 * t71 * t26 * t32 
     #* t42 * t25 * t70 * t3) - t98 * t60 * t3 * t25 * t46 + t100 * (t93
     # * t68 * t28 * t1 * t43 * t32 * t42 * t40 * t35 + t97 * t69 * t1 -
     # t96 * t54)
      t6 = 0.1q1 / t28
      t72 = 0.1q1 / t46
      t106 = 0.1q1 / t9
      t31 = 0.1q1 / t31
      t110 = t19 * t29
      t111 = t110 * t83 * t11
      t112 = t21 * t16
      t113 = t112 - t111
      t114 = t41 ** 2
      t115 = t34 ** 2
      t116 = t40 ** 2
      t117 = t40 * t116
      t118 = t47 ** 2
      t119 = t47 * t118
      t120 = t106 ** 2
      t121 = t93 * t36
      t122 = t121 * t32 * t42
      t123 = t59 * t39
      t124 = t56 * t36
      t125 = t59 * t36
      t126 = t38 * t36
      t127 = t43 * t26 * t72
      t128 = t43 * t27
      t129 = t10 * t26
      t130 = t21 * t43
      t131 = t19 * t18
      t132 = t56 * t39
      t133 = t31 * t106
      t134 = t133 * t26
      t135 = t134 * t38
      t136 = t31 ** 2
      t137 = t26 ** 2
      t138 = t38 ** 2
      t139 = t41 * t39
      t140 = t38 * t49
      t141 = t99 * t68
      t142 = t138 * t26
      t143 = t115 * t43
      t144 = t19 * t14
      t12 = t83 * t12
      t58 = t18 * (t21 * (t144 * t36 * t138 * t137 * t120 * t31 * t72 + 
     #t84 * (-t143 * t16 * t32 * t116 * t42 * t118 * t114 + t134 * t126 
     #* t32 * t42 * t41) * t21) + t12 * (t137 * (-t138 * t39 * t14 * t72
     # * t31 * t120 + t38 * t29 * t10 * t136 * t106) + t141 * (t139 * t2
     #7 * t6 + t140 * t29 * t47) - t142 * t27 * t14 * t120 * t31)) + t11
     #5 * (t131 * (t113 * t14 + t130 * t29) * t118 * t116 * t114 - t132 
     #* t41 * t8 * t40 * t118) + t34 * (t114 * (t99 * (-t124 * t27 * t6 
     #- t122 + t123) - t125 * t9 * t47 * t116) + (-t87 * t49 * t16 * t11
     #8 + t126 * t60 * t47) * t40 * t41) + t135 * (-t127 * t94 * t36 * t
     #32 * t42 * t106 + t56 * (t38 * (t106 * (t39 * (t127 + t45) + t128)
     # - t127 * t27 * t14 * t120) - t129 * t16 * t31) - t58 * t7 * t4)
      t94 = t112 * t84 * t32 * t42
      t127 = t55 * t36
      t145 = t66 * t83 * t11
      t146 = t40 * t118
      t147 = 0.1q1 / t41
      t148 = -t91 * t83 * t11 + t90
      t149 = t62 * t84 * t32 * t42
      t33 = t7 ** 2 * t32 * t33 + t4 ** 2 * t37 * t42
      t37 = t28 ** 2
      t150 = t49 ** 2
      t151 = t56 * t49 * t27
      t152 = t21 * (t61 - t149)
      t153 = t18 * t49
      t154 = t88 * t37 * t150
      t155 = t13 * t3
      t156 = t26 * t25
      t157 = t131 * t37
      t158 = t27 * t46
      t159 = t158 * t18
      t160 = t83 * t80
      t52 = t24 * (t28 * (-t12 * t88 * t39 * t10 - t151 * t1 - t33 * t49
     #) * t48 + t154 * t147 * (t85 * t101 * t32 * t42 + t12 * t79) * t78
     # + t153 * t28 * ((-t91 * t63 * t147 + t152) * t49 * t28 - t55 * t3
     #4 * t10 * t2) * t77) + t74 * (t157 * t150 * t148 * t48 + t154 * t1
     #9 * t63 * t77) + t159 * t3 * t25 * (t55 * (t50 * (t155 * t26 - t41
     #) + t156) * t10 * t27 + t83 * t52 * (-t26 * t15 * t3 * t50 + 1) * 
     #t10 + t85 * t28 * t1 * t32 * t42) + t160 * t68 * t37 * t150 * t24 
     #* t74 * t48
      t79 = 0.1q1 / t27
      t101 = t85 * t43 ** 2 * t32 * t42 + t12 * t14 ** 2
      t154 = t145 - t64
      t161 = -t130 * t84 * t32 * t42 + t144
      t162 = t110 - t94
      t163 = t104 * t72
      t164 = -t102 * t79 - t163
      t51 = t83 * t51
      t165 = t39 * t26
      t166 = t47 * t18
      t167 = t138 * t137
      t168 = t133 * t38
      t9 = t115 * (t146 * t90 * t1 * t18 * t162 + t166 * t116 * t35 * (-
     #t130 * t67 + t144 * t154) * t114) + t34 * (-t124 * t9 * t114 * t11
     #6 * t47 + t95 * (-t82 * t99 + t57 * (t39 * t36 * t24 * t48 * t17 +
     # t99 * t4 * t7) * t42)) + t156 * t123 * t126 * t3 + t90 * t38 * t2
     #6 * t27 * t18 * t161 * t31 * t120 + t167 * t27 * t72 * t18 * t101 
     #* t31 * t106 * t120 + t168 * (t26 * (t51 * t27 * t1 * t17 + t121 *
     # t164) * t42 * t32 + t131 * (-t165 * t21 * t164 - t158 * t148))
      t12 = t85 * t16 ** 2 * t32 * t42 + t12 * t29 ** 2
      t124 = t18 * t79
      t148 = t36 * t10
      t156 = t21 * t28
      t164 = t99 * t114
      t51 = t51 * t17
      t169 = t38 * t25
      t30 = t122 * t3 * t25 * t46 * (t26 * (t169 + t155) - t41) - t56 * 
     #t26 * (t39 * t46 * t13 * t25 * t75 + t133 * t41 * t10) + t30 * t59
     # * t46 * t25 * t75
      t155 = t99 * t41
      t4 = t115 * (t166 * t28 * t114 * t35 * t101 * t117 - t99 * t96 * t
     #80 * t39 * t35) + t27 * t30 + t34 * (-t81 * t95 * t28 * t14 * t116
     # * t47 - t155 * t33) + t25 * (t27 * (t132 * (-t104 + t44) + t83 * 
     #(-t82 * t46 + t57 * (-t104 * t39 * t17 + t7 * t4 * t46) * t42)) * 
     #t3 + t159 * t61 * t26 * t154 * t75 + t92 * t13 * t73 * t50 * t76) 
     #- t158 * t38 * (t132 * t26 + t158 * t59) * t3 * t70 + t131 * (t41 
     #* (t49 * (-t154 * t27 - t156 * t39) * t48 * t24 + t99 * t21 * (t14
     #8 * t115 * t35 + t102 * (t34 * t43 * t40 - t49))) + t38 * t27 * t1
     #0 * t72 * (-t144 * t83 * t11 + t130) * t31 * t120 * t137 + t133 * 
     #t129 * t102 * t21) + t28 * (t51 * (-t139 * t49 * t24 * t48 + t168 
     #* t129) + t93 * (-t141 * t1 + t168 * t104 + t164 * t49)) * t42 * t
     #32
      t7 = t83 * t11
      t2 = t18 * (t27 * (t25 * (t55 * t28 * t73 * t15 * t75 - t129 * t19
     # * t154 * t3) - t102 * t85 * t73 * t32 * t42 * t70 * t3) + t71 * (
     #-t104 * t85 * t46 * t32 * t42 * t70 * t3 + t152 * t104 * t46 * t50
     # * t75 * t25) + t100 * t10 * t28 * t19 * (t130 * t68 * t40 * t35 +
     # t7 * t69 * t19)) + t24 * (t56 * t37 * t34 * t15 * t13 * t150 * t1
     #47 * t78 - t107 * t28 * t150 * t18 * t63 * t77 + t34 * t28 * (-t12
     #5 * t10 + t56 * (t1 * t39 - t148)) * t48) - t95 * t28 * t150 * t54
     # * t48 * t74 + t109 * t28 * t32 * t42 * t48 * (t51 * t2 * t49 * t1
     #0 * t48 - t121 * t1) + (t73 * (t86 * t28 * t24 * t32 * t42 * t18 *
     # t3 - t93 * t28 * t13 * t32 * t42 * t75) - t33 * t3 * t46 + t81 * 
     #t83 * t28 * t150 * t24 * t48) * t25 * t27 + t98 * t86 * t18 * t46 
     #* t32 * t42 * t25 * t50 * t3
      t11 = t21 * t162
      t8 = t38 * (-t11 * t1 * t18 * t136 * t106 * t137 + t139 * t134 * t
     #60) + t142 * (t106 * (-t156 * t124 * t46 * t162 * t136 + t83 * (t5
     #7 * t36 * t42 * t17 - t80 * t39) * t31) - t123 * t45 * t31 * t120)
     # + t146 * t68 * (t123 * t34 * t8 + t153 * (-t107 * t113 * t6 + t11
     #) * t41)
      t11 = t93 * t1
      t13 = t28 * (t49 * (t24 * (t48 * (t28 * (t51 * t10 + t11) * t42 * 
     #t32 + t131 * (t7 * t107 * t10 + t156 * (-t38 * t1 * t147 + t10))) 
     #+ t89 * t1 * t147 * (t61 - t149) * t77) + t97 * t34 * t1 * t48 * t
     #74) - t158 * t56 * t10 * t25 * t3)
      t15 = (t122 * t115 * t1 * t35 * t47 + t160 * t115 * t29 * t10 * t1
     #18) * t40 * t41
      t3 = t114 * (t83 * t53 * t34 * t14 * t116 * t47 + t151 * t99) + t3
     #8 * (t158 * t125 * t25 * t3 + t133 * (t87 * t28 * t46 + t26 * t33)
     # + t128 * t11 * t137 * t32 * t42 * t120 * t31 * t72) + t18 * (t133
     # * (t129 * t91 * (t7 * (t165 * t72 + t27) * t19 - t64 * t26 * t72)
     # - t86 * (t105 + t103) * t42 * t32) + t127 * (t68 * t39 * t24 * t4
     #8 + t169 * (t158 + t165) * t3) + t26 * t38 * t21 * (-t163 * t107 *
     # t14 - t102 * t161) * t31 * t120) + t40 * (-t56 * t115 * t41 * t10
     # * t16 * t118 + t41 * (-t140 * t83 * t53 - t56 * t34 * t1 * (t34 *
     # t39 * t35 + t27) - t59 * t34 * t27 * t1) * t47) - t143 * t97 * t1
     #14 * t14 * t117 * t35 * t47 + t88 * t114 * t21 * (t43 * (-t156 * t
     #84 * t32 * t42 - t107) + t144 * t28) * t47 * t116 + t15
      t11 = -512
      t1 = t56 * (t34 * (-6 * t155 * t28 * t10 + 96 * t164 * t39) - 1024
     # * t115 * t114 * t49 * t16 * t29 * t6 * t40 * t119 + t135 * (t38 *
     # (-4096 * t26 * t29 * t46 * t16 * t79 * t136 - 192 * t36) + 48 * t
     #27 * t1))
      ret = t20 * t13 + t22 * t157 * t48 * t24 * t10 * t49 * (-t62 * t34
     # * t147 * t48 + t7 * (t147 * (t108 * t34 - t38) + t109) * t19) + 6
     #4 * t58 + 256 * t18 * (t38 * (-t127 * t23 * t106 * t136 * t137 + t
     #90 * t46 * (-t110 + t94) * t136 * t106 * t26) + t146 * t6 * t114 *
     # t115 * (t110 * (t64 - t145) + t112 * (t66 - t65)) + (t144 * (t112
     # - t111) + t130 * (t110 - t94)) * t136 * t120 * t137 * t138) + 32 
     #* t9 + t11 * (t138 * (t131 * t46 * t113 * t136 * t106 * t26 + t124
     # * (-t110 * t154 + t112 * t67) * t136 * t106 * t137) - t153 * t115
     # * t114 * t6 * t40 * t119 * t12 + t126 * t59 * t137 * t23 * t106 *
     # t136) + 8 * t4 - 128 * t8 + 16 * t3 + 2048 * t124 * t167 * t46 * 
     #t106 * t31 * t136 * t12 + 4 * t5 - t52 + 2 * t2 + t1

      hjetmass_triangle_pmpm_s23_mhsq_s14 = ret/32q0/(0,1q0)
      return

      end function
