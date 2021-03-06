
      complex*32 function hjetmass_bubble_pmpm_s14
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          double precision mt
          complex*32 ret

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i2)
      t5 = zb(i4, i2)
      t6 = zb(i2, i1)
      t7 = zb(i4, i3)
      t8 = za(i3, i4)
      t9 = 0.1q1 / t3
      t10 = 0.1q1 / t5
      t11 = 0.1q1 / t4
      t12 = t8 * t9
      t13 = t6 * t10
      t14 = t11 * t10
      t15 = za(i2, i3)
      t16 = za(i1, i2)
      t17 = t6 * t7
      t18 = t17 * t10
      t19 = -t18 + t2
      t20 = t10 ** 2
      t21 = t10 * t20
      t22 = t6 ** 2
      t23 = t6 * t22
      t24 = t22 * t11
      t25 = t24 * t20 * t19
      t26 = 0.1q1 / t6
      t27 = t9 ** 2
      t28 = t9 * t27
      t29 = t8 ** 2
      t30 = t8 * t29
      t31 = t1 * t6
      t32 = t8 * t2
      t33 = t32 + t31
      t34 = t16 * t5
      t35 = t3 * t7
      t36 = t35 + t34
      t37 = t16 * t6
      t38 = t3 * t2
      t39 = t1 * t5
      t40 = t8 * t7
      t41 = -t39 - t40 + t37 + t38
      t42 = 0.4q1
      t41 = t41 ** 2
      t43 = t42 * t33
      t44 = t43 * t36 + t41
      t44 = sqrt(t44)
      t45 = t44 - t39 - t40 + t37 + t38
      t46 = 0.1q1 / t33
      t47 = 0.1q1 / 0.2q1
      t48 = t47 * t1
      t49 = -t48 * t45 * t46 + t16
      t50 = t47 * t2
      t51 = -t50 * t45 * t46 - t7
      t52 = 0.1q1 / t33
      t53 = t16 ** 2
      t54 = t6 * t45
      t55 = t54 * t52
      t56 = t53 * t8
      t4 = t15 * t4
      t57 = t1 * t3
      t58 = t57 * t7
      t59 = t3 ** 2
      t60 = t2 * t59
      t61 = (-t37 + t4) * t3
      t62 = t16 * t8
      t63 = t4 * t62
      t64 = t3 * (t16 * (-t39 + t37 + t38 - t4) + t58)
      t65 = t39 + t40
      t66 = t62 * t45
      t67 = t66 * t52
      t68 = t52 ** 2
      t69 = t68 ** 2
      t70 = t52 * t68
      t71 = t45 ** 2
      t72 = t71 ** 2
      t73 = t45 * t71
      t74 = t62 * t71 * t68
      t75 = t40 + t4 + t39
      t76 = t71 * t68
      t77 = t76 * t8
      t78 = t1 * t8
      t79 = t78 * t73 * t70
      t80 = t62 * t3
      t81 = 0.2q1
      t82 = 0.3q1 / 0.4q1
      t83 = 0.9q1 / 0.4q1
      t41 = t42 * t33 * t36 + t41
      t41 = sqrt(t41)
      t84 = t41 - t39 - t40 + t37 + t38
      t85 = 0.1q1 / t8
      t86 = t3 * t85
      t87 = t47 * t45 * t46
      t88 = -t87 + t86
      t44 = t44 + t39 + t40 - t37 - t38
      t89 = t46 * (t45 + t44)
      t84 = 0.1q1 / t84
      t90 = t81 * t36
      t91 = -t90 * t84 - t87
      t92 = t16 * t2
      t93 = t1 * t7
      t94 = t93 + t92
      t95 = t16 * t3
      t96 = t35 + t34
      t97 = t34 * t52
      t98 = t3 * t96
      t99 = 0.3q1 / 0.8q1
      t100 = 0.1q1 / 0.4q1
      t101 = 0.1q1 / 0.8q1
      t102 = 0.1q1 / 0.16q2
      t54 = -t47 * t45 * (t94 * t46 * t59 - t56 * t54 * t68 + t95 * (-t3
     #9 - t4 + t37) * t46) - t82 * t74 * t65 - t100 * t76 * (t1 * ((-t37
     # + t4) * t3 + t60) + t63) - t99 * t62 * t73 * t70 * t33 + t101 * t
     #79 * t75 + t102 * t78 * t72 * t69 * t33 + t16 * (t8 * (t38 * t76 +
     # t97 * t45) - t98) + t35 * t67 * t81
      t103 = t2 * t5
      t104 = t103 + t17
      t105 = t46 * t5 * (t16 * (t40 - t4) - t58)
      t106 = t16 * t1
      t107 = t106 * t22
      t96 = t5 * t96
      t108 = t3 * t104
      t109 = t108 * t52
      t110 = t31 * t76
      t111 = t37 + t38
      t112 = t5 ** 2
      t113 = t1 ** 2
      t114 = t113 * t112
      t115 = t37 * t4
      t116 = t38 + t4
      t73 = t31 * t73 * t70
      t117 = t34 * t81
      t118 = t55 * t16
      t66 = -t100 * t76 * (t39 * (t40 + t4) + t114 - t115) - t101 * t73 
     #* t116 + t102 * t31 * t72 * t69 * t33 - t47 * t45 * (t66 * t104 * 
     #t68 + t107 * t71 * t70 + t105) + t82 * t37 * t76 * t111 + t16 * (t
     #109 * t45 - t110 * t5 + t96) + t118 * (-t77 * t99 * t2 + t117)
      t71 = t5 * t26
      t72 = -t87 - t71
      t76 = t57 * t85
      t87 = t16 - t76
      t119 = t3 * t6
      t120 = t119 * t85
      t121 = t85 ** 2
      t122 = t53 * t5
      t120 = -t81 * t106 * t59 * t85 * (t5 + t120) + t3 * (t3 * t59 * t1
     #13 * t6 * t85 * t121 + t59 * t113 * t5 * t121 + t120 * t53 + t122)
      t123 = t38 * t85 + t7
      t124 = 0.1q1 / t2
      t125 = t124 * t7 + t86
      t126 = t47 * t44 * t46
      t127 = -t126 - t86
      t128 = t53 * t6
      t129 = t103 * t26
      t130 = -t129 + t7
      t131 = t26 ** 2
      t132 = t95 * t7
      t133 = t5 * (t78 * t2 * t5 * t112 * t26 * t131 - t132 + t71 * (t16
     # * (-t40 + t38) - t58) + (t8 * (-t93 + t92) + t57 * t2) * t131 * t
     #112)
      t134 = 0.1q1 / t1
      t135 = t16 * t134
      t136 = t135 + t71
      t137 = -t126 + t71
      t138 = t39 * t26
      t139 = t138 + t16
      t48 = t48 * t44 * t46 + t16
      t50 = t50 * t44 * t46 - t7
      t140 = t44 ** 2
      t141 = t140 ** 2
      t142 = t44 * t140
      t143 = t68 * t140
      t144 = t143 * t62
      t145 = t78 * t142 * t70
      t46 = t100 * t143 * (t1 * (t3 * (t37 - t4) - t60) - t63) - t101 * 
     #t145 * t75 + t102 * t78 * t141 * t69 * t33 - t47 * t44 * (-t56 * t
     #6 * t44 * t68 + (-t93 - t92) * t46 * t59 + t95 * (t39 + t4 - t37) 
     #* t46) - t82 * t144 * t65 + t99 * t62 * t142 * t70 * t33 + t16 * (
     #t8 * (t143 * t38 - t97 * t44) - t98) - t80 * t81 * t7 * t44 * t52
      t41 = t41 + t39 + t40 - t37 - t38
      t41 = 0.1q1 / t41
      t90 = t90 * t41 + t126
      t97 = t1 * (-t61 - t60) - t63
      t98 = t62 * t44
      t126 = t143 * t31
      t142 = t31 * t142 * t70
      t146 = t37 * t52
      t147 = t146 * t44
      t69 = t100 * t143 * (t39 * (-t40 - t4) - t114 + t115) + t101 * t14
     #2 * t116 + t102 * t31 * t141 * t69 * t33 + t47 * t44 * (-t98 * t10
     #4 * t68 + t107 * t140 * t70 + t105) + t82 * t143 * t37 * t111 - t1
     #6 * (t109 * t44 + t126 * t5 - t96) + t147 * (t32 * t143 * t99 - t1
     #17)
      t70 = t34 * t36
      t96 = t7 ** 2
      t99 = t2 ** 2
      t100 = 0.8q1 * t2
      t101 = t22 * t20
      t102 = t103 + t17
      t103 = -t81 * t95 * t102 - t5 * (t16 * (-t40 + t4) + t58) - t128 *
     # t42 * t5
      t105 = 0.1q1 / t45
      t107 = 0.1q1 / t44
      t109 = -t105 + t107
      t140 = t107 ** 2
      t141 = t140 ** 2
      t148 = t107 * t140
      t149 = t33 ** 2
      t150 = t149 ** 2
      t151 = t33 * t150
      t152 = t33 * t149
      t153 = t105 ** 2
      t154 = t105 * t153
      t150 = 0.16q2 * t150
      t155 = 0.16q2 * t33
      t156 = t149 * t153
      t114 = t39 * (-t40 - t4) - t114 + t115
      t115 = 0.96q2 * t33
      t157 = t33 * t105
      t158 = 0.384q3 * t152
      t159 = t33 * t107
      t160 = t149 * t140
      t161 = t101 * t96
      t152 = t152 * t7
      t162 = t81 * t2
      t163 = 0.8q1 * t157
      t164 = t11 * t9
      t165 = -t105 + t107
      t166 = t33 * t7 * t165 - t2
      t167 = t105 * t107 * t149
      t168 = t167 * t7
      t169 = t157 * t7
      t170 = t2 + t169
      t171 = t71 + t86
      t172 = t8 * t5
      t173 = t81 * t16 * (-t172 + t119) + t57 * t5
      t174 = t95 * t81
      t175 = t16 * (-t40 + t4)
      t176 = t114 * t52
      t177 = t1 * t45 * t52
      t178 = t34 * t6
      t104 = t104 * t52
      t179 = t1 * t44 * t52
      t180 = t95 * t2
      t15 = 0.1q1 / t15
      t181 = t29 * t27
      t182 = t95 * t36
      t35 = -t3 * (t16 * (t39 + t4 - t37 - t38) - t58) - t62 * (t35 * t4
     #2 + t117)
      t117 = t16 * t42
      t71 = t5 * (t78 * t112 * t131 + t95 + t71 * (t62 + t57))
      t183 = t3 * (t31 * t59 * t121 - t34 + t86 * (t39 - t37))
      t184 = 0.1q1 / t171
      t125 = 0.1q1 / t125
      t185 = -0.1q1 / t72
      t136 = 0.1q1 / t136
      t186 = 0.1q1 / t127
      t187 = 0.1q1 / t137
      t188 = -0.1q1 / t88
      t171 = 0.1q1 / t171
      t189 = t136 ** 2
      t190 = t184 * t26
      t191 = t85 * t171
      t192 = t130 * t185
      t193 = t2 * t134
      t19 = -t181 * t121 * t68 * t11 * (-t7 * t160 * (t42 * t156 * (t42 
     #* (-0.6q1 * t62 * t65 + t81 * t97 + t62 * (t37 * t42 + 0.8q1 * t38
     #)) - t157 * (t115 * t182 * t105 + 0.32q2 * t35)) + t159 * (0.32q2 
     #* t156 * (t155 * t182 * t105 + t35 * t42) - t158 * t182 * t153 * t
     #107)) + t162 * (-0.64q2 * t151 * t182 * t153 * t140 * t109 + t150 
     #* t35 * t153 * t140)) - t9 * (t12 * t14 * t6 * t19 * t26 + t14 * t
     #29 * t7 * t27 - t25 * t26)
      t194 = -t14 * t9 * (t7 * (t13 - t12) - t2) + t43 * t6 * t96 * t11 
     #* t20 * t105 * t107 + 0.8q1 * t157 * t14 * t107 * t7 * t166
      t195 = t191 * t37 * t11 * (-t193 * (t16 + t139) - t7 - t130) * t20
     # * t136
      t24 = t24 * t53 * t130 * t139 * t134 ** 2 * (t192 * t52 * t187 + t
     #191) * t21 * t189
      t19 = t15 * (t1 * t194 + t131 * ((-t7 * (-0.6q1 * t16 * t23 * t20 
     #* (t119 * t10 - t8) + t42 * t23 * t173 * t21) + t162 * t101 * (-t1
     #74 * t6 + t173)) * t27 * t11 - t11 * (0.16q2 * t156 * t70 * t101 *
     # (-t100 * t18 + 0.6q1 * t101 * t96 + t81 * t99) * t140 + t161 * t1
     #60 * (t42 * t156 * (t42 * (-t42 * t62 * t102 + 0.6q1 * t37 * t111 
     #+ t114 * t81 - 0.8q1 * t31 * t34) + t157 * (t115 * t70 * t105 - 0.
     #32q2 * t103)) + t159 * (0.32q2 * t156 * (-t155 * t70 * t105 + t103
     # * t42) + t158 * t70 * t153 * t107)) * t68) + t42 * t25 * t7 * t68
     # * (0.64q2 * t151 * t70 * t153 * t140 * t109 + t150 * t103 * t153 
     #* t140)) + t16 * t19 + t96 * (t190 * t78 * t27 * t125 * t124 * t11
     # + t29 * t87 * t123 * t28 * t124 ** 2 * t11 * (t87 * t52 * t186 * 
     #t188 + t190) * t125 ** 2) + t40 * t190 * t11 * t27 * (t1 * t123 * 
     #t124 + t16 + t87) * t125 + t195 + t24)
      t23 = t192 * t146 * t15 * t11 * t20 * t136 * t187
      t24 = t187 ** 2
      t25 = t171 ** 2
      t53 = t185 ** 2
      t70 = t71 * t121
      t101 = t130 * t133 * t68 * t53 * t24 + t70 * t25
      t102 = t184 ** 2
      t103 = t188 ** 2
      t111 = t186 ** 2
      t114 = t130 ** 2
      t115 = t183 * t102 * t131
      t146 = t87 * t120 * t68 * t111 * t103
      t150 = t40 * t125
      t31 = t150 * (((t120 * (-t125 * t87 - t1) - t87 * (t4 * t3 * t87 +
     # t81 * (-t59 * t123 * t1 - t56 * t5 + t95 * (t40 + t38)) - 0.3q1 *
     # t3 * (t86 * t113 * t5 + t128) + t57 * (-t31 * t42 * t59 * t121 + 
     #0.7q1 * t95 * t6 * t85 + 0.5q1 * t34))) * t103 * t111 * t68 + t102
     # * t131 * (-t183 * t125 + t5 * (-t62 * t81 + t57) + t119 * (-t117 
     #+ 0.3q1 * t76))) * t124 * t123 + t115 + t146) * t11 * t28
      t3 = t15 * (t37 * t130 * t21 * t11 * t101 * t134 * t189 + t106 * (
     #0.8q1 * t168 * t164 * t109 - t42 * t164 * t157 * t107 * (t12 * t7 
     #+ t2)) + t37 * t11 * t21 * ((-t130 * (t3 * (t37 * t81 + t39) + t17
     #2 * (0.3q1 * t138 + t117)) + t2 * t71) * t25 * t121 - t114 * (t42 
     #* t5 * (t1 * (t32 * t131 * t112 + t34) + t180) - t81 * (t6 * (t132
     # - t122) + t1 * t112 * t26 * (-t39 + t40)) + t5 * (t4 * (t138 + t1
     #6) - t58) - 0.3q1 * t5 * (-t129 * t57 + t62 * t7) + 0.5q1 * t62 * 
     #t2 * t112 * t26) * t68 * t53 * t24) * t134 * t136 + t31 + t164 * t
     #6 * t20 * t94 - t14 * t8 * t94 * t27)
      t4 = t133 * t53 * t68 * t24
      t4 = t15 * (t9 * (-t30 * (t174 * t17 * t121 * t11 * t20 - t14 * t1
     #21 * (t10 * t173 * t7 + t180)) * t27 + t150 * t11 * t124 * t123 * 
     #(-t183 * t184 * t102 * t131 + t146 * (t186 + t188)) * t27 + t169 *
     # t113 * t11 * t107) + t37 * (t4 * t11 * (t185 + t187) * t134 * t21
     # * t136 * t114 + t11 * (t70 * t171 * t25 + t4 * t2) * t134 * t21 *
     # t136 * t130))
      t9 = t15 * t11
      t17 = -t18 + t2
      t18 = -0.1q1 / t137
      t21 = 0.1q1 / t89
      t24 = 0.1q1 / t88
      t25 = 0.1q1 / t91
      t31 = 0.1q1 / t72
      t32 = 0.1q1 / t90
      t34 = -0.1q1 / t127
      t37 = t31 ** 2
      t39 = t24 ** 2
      t40 = t49 * t54 * t121
      t53 = t51 ** 2
      t57 = t34 ** 2
      t59 = t21 ** 2
      t62 = t25 ** 2
      t68 = t18 ** 2
      t70 = t2 * t48
      t71 = t48 * t50
      t72 = t69 * t32
      t76 = t50 ** 2 * t131
      t86 = t76 * t68
      t87 = t49 * t51
      t88 = t51 * (t51 * t66 * t131 * t37 + t40 * t39)
      t89 = t88 * t154 * t84
      t90 = t148 * t32 * t41
      t44 = t90 * (t86 * (t72 + t42 * t178 * (t179 + t16) - t47 * t142 *
     # t33 + t81 * t16 * (t104 * t44 * t8 + t108) - t82 * t126 * t116 - 
     #0.3q1 * t147 * (t6 * (t179 + t16) + t38) - t5 * (-t58 - t175) - t1
     #76 * t44 - t144 * t83 * t6 * t2) + (t46 * (t50 * (t32 * t48 - t1) 
     #- t70) - t71 * (t97 * t52 * t44 + t42 * t80 * (t2 * t44 * t52 - t7
     #) + t47 * t145 * t33 + t81 * t56 * (t6 * t44 * t52 - t5) - t82 * t
     #143 * t78 * t75 + t83 * t144 * t33 + t64 - 0.3q1 * t98 * t52 * t65
     #)) * t121 * t57)
      t78 = t9 * t59 * t36
      t91 = t78 * t33
      t5 = t91 * (t89 * t62 + (-t53 * (t42 * t178 * (-t16 + t177) - t47 
     #* t73 * t33 + t81 * t16 * (t104 * t45 * t8 - t108) + t82 * t110 * 
     #t116 - t5 * (t58 + t175) - t176 * t45 - 0.3q1 * t118 * (t6 * (-t17
     #7 + t16) + t38) + t74 * t83 * t6 * t2) * t131 * t37 + (t54 * (-t1 
     #* t51 - t2 * t49) - t87 * (-t42 * t80 * (t2 * t45 * t52 + t7) - t4
     #7 * t79 * t33 - t81 * t56 * (t5 + t55) + 0.3q1 * t67 * t65 + t83 *
     # t74 * t33 - t82 * t77 * t1 * t75 + t64 + (t1 * (t61 + t60) + t63)
     # * t52 * t45)) * t121 * t39) * t154 * t84 * t25 + t44)
      t6 = t84 ** 2
      t8 = t9 * t21
      t6 = t8 * t36 ** 2 * t149 * (t49 ** 2 * t51 * t85 * t154 * t6 * t2
     #4 * t62 + t71 * t148 * t41 ** 2 * t32 ** 2 * (t50 * t26 * t18 + t4
     #8 * t85 * t34) + t49 * t53 * t26 * t154 * t6 * t31 * t62)
      t21 = t9 * t21 * t59 * t36 * t33 * (t71 * t90 * t46 * t121 * t57 +
     # t90 * t86 * t69 - t89 * t25)
      t38 = t51 * t31
      t24 = t91 * (t90 * t76 * t69 * t18 * t68 + t90 * (t48 * t46 * t34 
     #* t57 * t121 - t2 * t69 * t68 * t131) * t50 + t51 * t154 * t84 * t
     #25 * (t40 * t24 * t39 + t66 * t131 * t37 * (t38 - t2)))
      t34 = t78 * t149 * (t71 * t46 * t121 * t141 * t41 * t57 * t32 - t8
     #8 * t153 ** 2 * t84 * t25 + t72 * t86 * t141 * t41)
      t18 = t50 * t140 * t41 * t18 * t32
      t8 = t8 * t26 * t36 * t33
      t26 = t8 * (-t1 * t53 * t153 * t84 * t31 * t25 + t18 * (-t1 * t50 
     #+ t70) + t87 * t2 * t153 * t84 * t31 * t25)
      t13 = t16 * t15 * (t16 * t164 * (t52 * t159 * (t81 * t156 * (t155 
     #* t7 * t105 + t100) + t159 * (-0.8q1 * t157 * (t43 * t7 * t105 + t
     #162) + 0.32q2 * t149 * t7 * t105 * t107)) - t12 * t81 * (-t42 * t2
     # * t149 * t105 * t107 + 0.8q1 * t152 * t105 * t107 * t109) * t52 +
     # t163 * t29 * t7 * t27 * t107) + t14 * (-t52 * (-0.32q2 * t160 * (
     #t167 * t96 - t169 * t170) - t159 * t81 * (t43 * t99 * t105 + 0.16D
     #2 * t156 * t7 * t170)) + t163 * t161 * t107 + 0.16q2 * t13 * t168 
     #* t166 * t52))
      ret = -t42 * t13 + t181 * t15 * (-0.32q2 * t12 * t16 + 0.16q2 * t1
     #) * (0.64q2 * t152 * t182 * t121 * t11 * t153 * t140 * t109 - 0.16
     #q2 * t156 * t121 * t11 * t140 * (t182 * t2 + t35 * t7)) + 0.20q2 *
     # t23 * (t193 * t16 + t7) - 0.48q2 * t9 * (t135 * t22 * t130 * t20 
     #** 2 * t136 * t101 + (t16 * t28 * t10 - t123 * t125 * t124 * (t146
     # + t115) * t27 ** 2) * t29 * t7) + 0.256q3 * t168 * t14 * t15 * (t
     #92 * t165 + t93 * t165) - 0.128q3 * t9 * t157 * t107 * (t93 * t10 
     #* t17 - t169 * t30 * t28 * (-t1 * t42 + 0.6q1 * t12 * t16) * t182 
     #* t121 * t107 + t92 * t10 * t17) - 0.1024q4 * t5 - 0.512q3 * t6 + 
     #0.12q2 * t23 * (t193 * t139 + t130) - 0.4096q4 * t21 - 0.2048q4 * 
     #t24 - 0.6144q4 * t34 - 0.192q3 * t26 - 0.320q3 * t8 * (t38 * t153 
     #* t84 * t25 * t94 + t18 * t94) - 0.8q1 * t19 - 0.16q2 * t3 - 0.32D
     #2 * t4

      hjetmass_bubble_pmpm_s14 = -ret/16q0/(0q0,1q0)
      return

      end function
