
      complex*32 function hjetmass_bubble_pmpm_s23
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
      t2 = za(i2, i4)
      t3 = zb(i2, i1)
      t4 = za(i1, i3)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = zb(i4, i2)
      t8 = za(i3, i4)
      t9 = zb(i3, i1)
      t10 = zb(i4, i3)
      t11 = 0.1q1 / t1
      t12 = 0.1q1 / t9
      t13 = t4 * t11
      t14 = t3 * t12
      t15 = t14 + t13
      t16 = 0.1q1 / t10
      t17 = t7 * t16
      t18 = t13 + t17
      t19 = 0.1q1 / t2
      t20 = t19 * t8 + t17
      t21 = t2 * t10
      t22 = t1 * t9
      t23 = t21 + t22
      t24 = t4 * t3
      t25 = t8 * t7
      t26 = t25 + t24
      t27 = t1 * t3
      t28 = t4 * t9
      t29 = t2 * t7
      t30 = t8 * t10
      t31 = -t30 + t27 - t28 + t29
      t32 = 0.4q1
      t31 = t31 ** 2
      t33 = t32 * t26 * t23 + t31
      t33 = sqrt(t33)
      t34 = -t30 + t27 - t28 + t29 + t33
      t33 = -t30 + t27 - t28 + t29 - t33
      t35 = t9 * t7
      t36 = t35 * t16
      t37 = t36 - t3
      t38 = 0.1q1 / t23
      t39 = 0.1q1 / 0.2q1
      t40 = t39 * t34 * t38
      t41 = -t17 + t40
      t42 = t39 * t33 * t38
      t43 = -t17 + t42
      t44 = t4 * t8
      t45 = t44 * t26
      t46 = t25 * t26
      t47 = t2 * t4
      t48 = t47 * t11
      t49 = t48 - t8
      t50 = t28 * t11 + t3
      t51 = t29 * t16
      t52 = t51 + t8
      t53 = t1 * t8
      t54 = t53 + t47
      t55 = t7 ** 2
      t56 = t16 ** 2
      t57 = t1 * t2
      t58 = t7 * (t57 * t55 * t56 + t17 * t54 + t44)
      t59 = 0.2q1
      t60 = t8 * t32
      t61 = t4 * t10
      t62 = t1 * t7
      t63 = t11 ** 2
      t64 = t4 ** 2
      t65 = t21 * t64 * t63
      t66 = t4 * (t65 - t25 + t13 * (-t30 + t29))
      t67 = t61 * t11
      t68 = t64 * t2
      t69 = t8 ** 2
      t70 = t2 ** 2
      t71 = t69 * t7
      t67 = -t59 * t68 * t8 * t11 * (t67 + t7) + t4 * (t4 * t64 * t70 * 
     #t10 * t11 * t63 + t64 * t70 * t7 * t63 + t67 * t69 + t71)
      t72 = t13 + t40
      t73 = t13 + t42
      t74 = t47 * t3
      t75 = t44 * t3
      t17 = t7 * (t57 * t9 * t7 * t55 * t16 * t56 + t17 * (t8 * (t28 - t
     #27) - t74) - t75 + (-t57 * t3 + t54 * t9) * t56 * t55)
      t54 = t44 * t9
      t76 = t75 * t10
      t71 = t71 * t10
      t77 = t5 * t6
      t78 = -t40 * t9 + t3
      t79 = 0.1q1 / t23
      t80 = t3 * t10
      t81 = t35 + t80
      t82 = t79 ** 2
      t83 = t82 ** 2
      t84 = t79 * t82
      t85 = t10 ** 2
      t86 = t34 ** 2
      t87 = t86 ** 2
      t88 = t34 * t86
      t89 = t38 * t7 * (t8 * (-t77 + t27) - t74)
      t90 = t53 * t81 * t82
      t91 = t2 * t8 * t85
      t92 = t25 + t24
      t93 = t7 * t92
      t94 = t4 * t81
      t95 = t94 * t79
      t96 = t29 * t10
      t97 = t77 * t30
      t98 = t70 * t55
      t99 = t29 * (t77 + t27) - t97 + t98
      t100 = t86 * t82
      t101 = t30 + t28
      t102 = t77 + t28
      t103 = t21 * t88 * t84
      t104 = 0.3q1 / 0.4q1
      t105 = 0.1q1 / 0.4q1
      t106 = 0.1q1 / 0.8q1
      t107 = 0.1q1 / 0.16q2
      t108 = 0.3q1 / 0.8q1
      t109 = t25 * t59
      t110 = t30 * t34 * t79
      t86 = t39 * t34 * (t91 * t86 * t84 - t90 * t34 + t89) - t105 * t10
     #0 * t99 + t104 * t100 * t30 * t101 + t106 * t103 * t102 + t107 * t
     #21 * t87 * t83 * t23 + t8 * (-t96 * t86 * t82 - t95 * t34 + t93) +
     # t110 * (t100 * t22 * t108 - t109)
      t31 = t32 * t26 * t23 + t31
      t31 = sqrt(t31)
      t111 = t31 - t30 + t27 - t28 + t29
      t111 = 0.1q1 / t111
      t112 = t59 * t26
      t113 = -t112 * t111 - t40
      t114 = t38 * (t34 - t33)
      t115 = t2 * t34 * t79
      t97 = t29 * (t77 + t27) - t97 + t98
      t98 = t7 * (t8 * (t77 - t27) + t74)
      t116 = t97 * t79
      t81 = t81 * t79
      t117 = t25 * t10
      t118 = -0.9q1 / 0.4q1
      t119 = t100 * t53
      t40 = -t40 * t2 - t8
      t120 = -t42 * t2 - t8
      t121 = -t42 * t9 + t3
      t31 = -t31 - t30 + t27 - t28 + t29
      t31 = 0.1q1 / t31
      t42 = -t112 * t31 - t42
      t112 = t2 * t33 * t79
      t122 = t30 * t33 * t79
      t123 = t33 ** 2
      t124 = t123 ** 2
      t125 = t33 * t123
      t126 = t21 * t125 * t84
      t127 = t123 * t82
      t89 = t104 * t127 * t30 * t101 - t105 * t127 * t99 + t106 * t126 *
     # t102 + t107 * t21 * t124 * t83 * t23 + t39 * t33 * (t91 * t123 * 
     #t84 - t90 * t33 + t89) + t8 * (-t127 * t96 - t95 * t33 + t93) + t1
     #22 * (t127 * t22 * t108 - t109)
      t5 = 0.1q1 / t5
      t90 = t2 * t3
      t91 = t8 * t9
      t93 = -t91 + t90
      t95 = 0.1q1 / t7
      t6 = 0.1q1 / t6
      t96 = t95 ** 2
      t99 = t96 ** 2
      t128 = t95 * t96
      t129 = t2 * t9
      t130 = t8 * t3
      t131 = t6 * t95
      t132 = -t59 * t8 * (-t61 + t62) + t47 * t7
      t133 = t30 * t96
      t134 = t32 * t10
      t135 = t59 * t9
      t136 = t5 * t6
      t137 = t90 * t59
      t138 = 0.1q1 / t33
      t139 = 0.1q1 / t34
      t140 = t138 + t139
      t141 = t138 ** 2
      t142 = t138 * t141
      t143 = t23 ** 2
      t144 = t143 ** 2
      t145 = t23 * t144
      t146 = t23 * t143
      t147 = t139 * t138
      t148 = t147 * t146 * (t139 * t140 + t141)
      t149 = t138 + t139
      t150 = t147 * t143
      t151 = t150 * t149
      t152 = 0.16q2
      t153 = t147 * t32
      t24 = t4 * (t8 * (t77 - t28 + t29 - t30) - t74) + t53 * (t24 * t32
     # + t109)
      t109 = t139 ** 2
      t154 = t109 ** 2
      t155 = t139 * t109
      t156 = t109 * t6 * t141
      t157 = t156 * t143
      t35 = t35 + t80
      t158 = t59 * t44 * t35 + t7 * (t8 * (t77 - t27) + t74) + t71 * t32
      t159 = t77 + t27 + t29
      t160 = t57 * t125 * t84
      t161 = t91 + t90
      t162 = t161 * t38 * t64
      t38 = t44 * (-t77 - t29 + t30) * t38
      t163 = t69 * t10
      t164 = t163 * t1
      t165 = t25 * t79
      t92 = t4 * t92
      t166 = t27 + t29
      t167 = t127 * t53
      t64 = t9 * t64
      t168 = t77 * t53
      t75 = t75 * t59 * t1
      t124 = t107 * t57 * t124 * t83 * t23 + t108 * t53 * t125 * t84 * t
     #23 - t106 * t160 * t159 - t104 * t167 * t166 - t105 * t127 * (t2 *
     # (t4 * (t77 - t30) + t64) + t168) + t39 * t33 * (t164 * t33 * t82 
     #+ t162 + t38) - t8 * (t1 * (-t127 * t28 + t165 * t33) + t92) - t75
     # * t33 * t79
      t125 = t44 * t1
      t169 = t2 * (t4 * (t77 - t30) + t64) + t168
      t170 = t4 * (t8 * (t77 - t28 + t29 - t30) - t74)
      t171 = t169 * t79
      t172 = t1 * t69
      t173 = t57 * t88 * t84
      t38 = -t104 * t119 * t166 + t105 * t100 * (t2 * (t4 * (-t77 + t30)
     # - t64) - t168) - t106 * t173 * t159 + t107 * t57 * t87 * t83 * t2
     #3 + t108 * t53 * t88 * t84 * t23 + t39 * t34 * (t164 * t34 * t82 +
     # t162 + t38) - t8 * (t1 * (-t100 * t28 + t165 * t34) + t92) - t75 
     #* t34 * t79
      t64 = t131 * t5
      t75 = t80 * t95
      t80 = t147 * t136
      t83 = t3 * t23
      t84 = t3 ** 2
      t87 = t9 ** 2
      t88 = t6 * t82
      t92 = t152 * t144
      t105 = t23 * t139
      t106 = 0.96q2 * t105
      t107 = t109 * t143
      t108 = t105 * t152
      t139 = t23 * t138
      t144 = t141 * t143
      t162 = 0.1q1 / t4
      t164 = 0.1q1 / t18
      t165 = 0.1q1 / t72
      t168 = 0.1q1 / t43
      t18 = 0.1q1 / t18
      t174 = 0.1q1 / t41
      t15 = 0.1q1 / t15
      t175 = 0.1q1 / t73
      t20 = 0.1q1 / t20
      t176 = t162 ** 2
      t177 = t162 * t176
      t178 = t18 ** 2
      t179 = t37 ** 2
      t180 = t165 ** 2
      t181 = t174 ** 2
      t182 = t175 ** 2
      t183 = t168 ** 2
      t184 = t1 ** 2
      t185 = t164 ** 2
      t186 = t49 * t67 * t82 * t180 * t182
      t14 = t14 * t15
      t156 = t5 * (t3 * t9 * (-0.64q2 * t156 * t46 * t146 * t96 * t149 +
     # t157 * t96 * (-0.32q2 * t46 * t10 * t95 + t152 * t158)) + t30 * (
     #(t9 * t17 * t183 * t181 * t82 + t58 * t164 * t185 * t63) * t20 * t
     #19 * t6 * t128 * t37 + t88 * t17 * t183 * t181 * (t174 + t168) * t
     #20 * t19 * t128 * t179)) + t83 * t80 * t184 * t69 * t177 + (t5 * (
     #-t152 * t157 * (t45 * (-t91 + t90) + t130 * t24) - 0.64q2 * t156 *
     # t130 * t45 * t146 * t149 + t14 * t6 * t50 * (t66 * t18 * t178 * t
     #56 - t186 * (t175 + t165))) + t136 * t76 * t59 * t96 - t64 * (t3 *
     # t132 * t95 + t54)) * t177 * t1
      t157 = t20 ** 2
      t187 = t66 * t178 * t56
      t188 = t15 * t3
      t189 = t58 * t185 * t63
      t22 = t30 * t128 * (-t179 * (-t32 * t7 * (t2 * (t22 * t56 * t55 + 
     #t25) + t54) + t59 * (t2 * t55 * t16 * (t27 - t29) + t76 - t71) - 0
     #.3q1 * t7 * (-t53 * t3 + t47 * t36) - t7 * (t77 * (t51 + t8) - t74
     #) - 0.5q1 * t53 * t9 * t55 * t16) * t82 * t181 * t183 + (-t37 * (t
     #4 * (t30 * t59 + t29) + t62 * (t60 + 0.3q1 * t51)) - t58 * t9) * t
     #63 * t185) * t20 + t30 * t37 * t128 * (-t37 * t17 * t82 * t181 * t
     #183 + t189) * t157
      t4 = t136 * (t162 * (t1 * (t95 * t161 * t162 + t188 * (((t49 * (t5
     #9 * (t62 * t69 + t68 * t50 - t44 * (t27 + t28)) + 0.3q1 * t4 * (t1
     #3 * t70 * t7 + t163) + t77 * t4 * t49 + t47 * (-0.7q1 * t44 * t10 
     #* t11 + t65 * t32 - 0.5q1 * t25)) + t67 * (-t15 * t49 + t2)) * t18
     #2 * t82 * t180 + t178 * t56 * (t66 * t15 - t7 * (-t53 * t59 + t47)
     # - t61 * (0.3q1 * t48 - t60))) * t12 * t50 + t186 - t187) * t176) 
     #- t10 * t96 * t161) + t19 * t22)
      t13 = t18 * t16
      t18 = t9 * t19
      t22 = -t2 * (0.8q1 * t83 * t147 * t131 * (-t83 * t140 + t9) - t147
     # * t134 * t84 * t23 * t6 * t96) + t84 * t6 * t96 * (t82 * t144 * (
     #t32 * t107 * (t32 * (-t32 * t53 * t35 + 0.6q1 * t30 * t101 - 0.8q1
     # * t21 * t25 - t59 * t97) + t105 * (t106 * t46 - 0.32q2 * t158)) +
     # t139 * (-0.32q2 * t107 * (-t108 * t46 + t158 * t32) + 0.384q3 * t
     #46 * t146 * t109 * t138)) - t134 * (-0.64q2 * t46 * t145 * t109 * 
     #t141 * t149 + t92 * t158 * t109 * t141) * t95 * t82 + 0.96q2 * t10
     #7 * t46 * t85 * t96 * t141) + (-t6 * (0.32q2 * t107 * t129 * t45 *
     # t141 + t130 * t144 * (t32 * t107 * (t32 * (-0.6q1 * t53 * t166 - 
     #t169 * t59 + t53 * (t30 * t32 + 0.8q1 * t28)) - t105 * (t106 * t45
     # + 0.32q2 * t24)) - t139 * (0.32q2 * t107 * (t108 * t45 + t24 * t3
     #2) + 0.384q3 * t45 * t146 * t109 * t138)) * t82) + t59 * t88 * t93
     # * (0.64q2 * t45 * t145 * t109 * t141 * t149 + t92 * t24 * t109 * 
     #t141)) * t176
      t24 = t164 * t11
      t25 = t37 * t168 * t174 * t79
      t12 = -t162 * (t136 * (t3 * (-0.6q1 * t133 * (-t61 * t95 + t1) - t
     #134 * t132 * t128) + t135 * t96 * (-t44 * t59 * t10 + t132)) * t16
     #2 + t131 * (-t5 * t93 * t10 * t95 - t130 * t85 * t5 * t96 + t129 *
     # t5)) + t5 * t22 + t6 * (t13 * t188 * t5 * t176 * t1 * (t49 - t2 *
     # t12 * (t50 + t3) - t8) + t184 * t84 * t49 * t50 * t177 * t5 * t12
     # ** 2 * (-t49 * t165 * t79 * t175 + t13) * t15 ** 2) + t176 * (-t6
     #4 * (t8 * (t75 - t9) + t90) - 0.8q1 * t80 * t130 * t23 * (t8 * t23
     # * t140 + t2) + t153 * t136 * t69 * t9 * t23) * t1 + t24 * t136 * 
     #t133 * (t18 * (t52 + t8) + t3 - t37) * t20 + t136 * t69 * t37 * t5
     #2 * t19 ** 2 * (t24 - t25) * t157 * t128 * t85
      t13 = -0.1q1 / t72
      t15 = -0.1q1 / t43
      t22 = -0.1q1 / t41
      t24 = 0.1q1 / t42
      t27 = 0.1q1 / t113
      t29 = 0.1q1 / t114
      t30 = -0.1q1 / t73
      t35 = t121 ** 2
      t36 = t27 ** 2
      t41 = t24 ** 2
      t42 = t31 ** 2
      t43 = t78 * t22
      t44 = t40 * t78
      t47 = t136 * t29
      t48 = t13 ** 2
      t49 = t22 ** 2
      t51 = t40 * t38 * t63
      t54 = t29 ** 2
      t55 = t78 ** 2
      t58 = t30 ** 2
      t61 = t15 ** 2
      t62 = t9 * t120
      t65 = t120 * t121
      t66 = t35 * t56
      t67 = t66 * t61
      t68 = t78 * (t78 * t86 * t56 * t49 + t51 * t48) * t155 * t111
      t70 = t24 * t142 * t31
      t33 = t70 * (t67 * (-t89 * t24 - t104 * t21 * t123 * t82 * t102 + 
     #t32 * t117 * (t112 + t8) - t39 * t126 * t23 + t59 * t8 * (t81 * t3
     #3 * t1 + t94) - 0.3q1 * t122 * (t10 * (t112 + t8) + t28) + t116 * 
     #t33 + t98 + t53 * t118 * t9 * t10 * t123 * t82) + (t124 * (t121 * 
     #(-t120 * t24 + t2) + t62) + t65 * (t104 * t127 * t57 * t159 + t118
     # * t167 * t23 + t32 * t125 * (-t9 * t33 * t79 + t3) - t39 * t160 *
     # t23 - t59 * t172 * (t10 * t33 * t79 - t7) + t171 * t33 + t170 + 0
     #.3q1 * t53 * t33 * t79 * t166)) * t63 * t58)
      t71 = t136 * t54 * t26
      t72 = t71 * t23
      t1 = t72 * (-t68 * t36 + (t55 * (-t104 * t100 * t21 * t102 + t32 *
     # t117 * (t115 + t8) - t39 * t103 * t23 + t59 * t8 * (t81 * t34 * t
     #1 + t94) - 0.3q1 * t110 * (t10 * (t115 + t8) + t28) + t116 * t34 +
     # t98 + t119 * t118 * t9 * t10) * t49 * t56 + (t38 * (t2 * t78 + t4
     #0 * t9) + t44 * (t104 * t100 * t57 * t159 + t118 * t119 * t23 + t3
     #2 * t125 * (-t9 * t34 * t79 + t3) - t39 * t173 * t23 - t59 * t172 
     #* (t10 * t34 * t79 - t7) + t171 * t34 + t170 + 0.3q1 * t53 * t34 *
     # t79 * t166)) * t63 * t48) * t155 * t111 * t27 + t33)
      t7 = t136 * t29 * t54 * t26 * t23 * (t70 * t65 * t124 * t63 * t58 
     #+ t70 * t67 * t89 - t68 * t27)
      t10 = t86 * t49 * t56
      t21 = t120 * t124
      t28 = t72 * (t70 * t66 * t89 * t15 * t61 + t78 * t111 * t155 * t27
     # * (t51 * t13 * t48 + t10 * (t43 - t9)) + t70 * (t21 * t30 * t58 *
     # t63 - t56 * t61 * t89 * t9) * t121)
      t29 = t121 * t31
      t10 = t71 * t143 * (t44 * t38 * t63 * t111 * t154 * t48 * t27 + t1
     #0 * t55 * t111 * t154 * t27 + t29 * t141 ** 2 * t24 * (t121 * t56 
     #* t61 * t89 + t21 * t58 * t63))
      t21 = t29 * t141 * t15 * t24
      t24 = t47 * t16 * t26 * t23
      t22 = t24 * (t2 * t55 * t111 * t109 * t22 * t27 + t21 * (-t121 * t
     #2 + t62) - t44 * t9 * t111 * t109 * t22 * t27)
      t14 = t136 * ((t17 * t179 * t181 * t183 * t19 * t20 * t82 * t99 - 
     #t189 * t19 * t20 * t37 * t99) * t85 * t8 + t14 * t184 * t50 * t176
     # ** 2 * (-t187 + t186))
      t17 = t25 * t136 * t133 * t20
      t19 = t17 * (-t18 * t52 + t37)
      t20 = t24 * (t109 * t111 * t161 * t27 * t43 - t161 * t21)
      t21 = t64 * t150 * (-t147 * t46 * t87 * t95 + t130 * t140 * t9 + t
     #140 * t2 * t84)
      ret = t152 * t4 + t32 * ((-t152 * t136 * t8 * (-t137 + t91) * t151
     # + 0.32q2 * t136 * t69 * t3 * t148 + t153 * t136 * t2 * (-t60 * t9
     # + t137) * t23) * t162 + t8 * (0.32q2 * t131 * t84 * t148 + t153 *
     # t23 * (t59 * t131 * (t84 * t85 * t96 + t87) - t134 * t3 * t9 * t6
     # * t96) - t152 * t131 * t3 * (t135 - t75) * t151) * t5) + 0.512q3 
     #* t47 * t26 ** 2 * t143 * (t120 ** 2 * t121 * t11 * t42 * t142 * t
     #30 * t41 - t44 * t111 ** 2 * t155 * t36 * (t11 * t13 * t40 + t16 *
     # t43) + t120 * t35 * t16 * t42 * t142 * t15 * t41) + 0.1024q4 * t1
     # - 0.4096q4 * t7 - 0.2048q4 * t28 + 0.6144q4 * t10 + 0.192q3 * t22
     # - 0.48q2 * t14 - 0.128q3 * t64 * t147 * t23 * (t161 * t75 - t161 
     #* t9) - 0.20q2 * t17 * (t18 * t8 + t3) + 0.12q2 * t19 + 0.320q3 * 
     #t20 - 0.256q3 * t21 + 0.32q2 * t156 + 0.8q1 * t12 + t136 * t130 * 
     #t184 * t177 * (0.768q3 * t107 * t45 * t162 * t141 - 0.40q2 * t95)

      hjetmass_bubble_pmpm_s23 = ret/16q0*(0,1q0)
      return

      end function
