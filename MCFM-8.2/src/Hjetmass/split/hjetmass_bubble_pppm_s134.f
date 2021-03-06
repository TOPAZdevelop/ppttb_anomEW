
      complex*32 function hjetmass_bubble_pppm_s134
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt
          real*16 p(mxpart,4)
          complex*32 alpha
          logical flip

      alpha = (za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1) +
     & za(i3,i4)*zb(i4,i3))/(za(i1,i3)*zb(i3,i1) + za(i3,i4)*zb(i4,i3))

      p(5,:) = real(alpha)*p(i3,:)
      p(6,:) = (1q0-real(alpha))*p(i3,:) + p(i1,:) + p(i4,:)
      if (flip .eqv. .true.) then
          call spinoru_qp(6,p,zb,za)
      else
          call spinoru_qp(6,p,za,zb)
      end if

      t1 = za(6, i4)
      t2 = zb(i2, i1)
      t3 = zb(i3, i1)
      t4 = za(6, i2)
      t5 = za(6, i3)
      t6 = za(i1, i4)
      t7 = za(i3, i4)
      t8 = zb(i4, i1)
      t9 = zb(i4, i3)
      t10 = zb(i1, 6)
      t11 = za(i2, i4)
      t12 = za(i1, i3) * t3
      t13 = t6 * t8
      t14 = t7 * t9
      t15 = t14 + t13 + t12
      t16 = zb(i3, 6)
      t17 = za(5, i4)
      t18 = za(5, i2)
      t19 = 0.1q1 / t18
      t20 = t17 * t4
      t21 = t20 * t19 - t1
      t22 = zb(i1, 5)
      t23 = t4 * t19
      t24 = t23 * t10 + t22
      t25 = zb(i2, 5)
      t26 = zb(i2, 6)
      t27 = 0.1q1 / t26
      t28 = t25 * t27
      t29 = t23 + t28
      t30 = zb(i4, 6)
      t31 = zb(i4, 5)
      t32 = t5 * t16
      t33 = t17 * t31
      t34 = t1 * t30
      t35 = t34 + t32 - t33
      t35 = 4 * t34 * t33 + t35 ** 2
      t35 = sqrt(t35)
      t36 = -t33 + t34 + t32 - t35
      t37 = 0.1q1 / t30
      t38 = 0.1q1 / t17
      t39 = (0.1q1 / 0.2q1)
      t40 = t36 * t38
      t41 = t40 * t39 * t37
      t42 = -t41 + t23
      t35 = -t33 + t34 + t32 + t35
      t43 = t35 * t38
      t44 = t43 * t39 * t37
      t45 = -t44 + t23
      t46 = t23 * t26 + t25
      t47 = zb(i3, i2)
      t48 = t17 * t25
      t49 = t48 * t27
      t50 = t49 + t1
      t51 = t10 * t25
      t52 = t51 * t27 - t22
      t53 = -t41 - t28
      t54 = -t44 - t28
      t55 = t36 * t37
      t56 = t55 * t39 - t1
      t57 = t41 * t10 + t22
      t58 = 0.1q1 / t36
      t59 = 2
      t60 = t1 * t31
      t61 = t60 * t59
      t62 = t61 * t58 + t41
      t63 = t38 * t37
      t64 = t63 * (t36 - t35)
      t65 = za(6, i1)
      t66 = t65 * t22
      t67 = t4 * t25
      t68 = t60 + t66 + t67
      t69 = za(5, i1)
      t70 = t69 * t10
      t71 = t18 * t26
      t72 = t17 * t30
      t73 = t72 + t70 + t71
      t74 = t69 * t22
      t75 = t65 * t10
      t76 = t18 * t25
      t77 = t4 * t26
      t78 = -t76 + t77 + t75 - t74 - t33 + t34 + t32
      t78 = t78 ** 2
      t79 = 4 * t68 * t73 + t78
      t79 = sqrt(t79)
      t80 = -t76 + t77 + t75 - t74 - t33 + t34 + t32 + t79
      t81 = 0.1q1 / t73
      t82 = t80 * t81
      t83 = t82 * t39
      t84 = t83 * t17 - t1
      t85 = t83 * t10 + t22
      t86 = t83 * t26 + t25
      t87 = 0.1q1 / t73
      t88 = t71 + t70
      t89 = t88 * t87
      t90 = t89 * t80
      t91 = za(i1, i2) * t2
      t92 = za(i2, i3) * t47
      t93 = t11 * zb(i4, i2)
      t94 = t1 * t26
      t95 = t1 * t10
      t96 = t17 * t22
      t97 = (t48 + t94) * t4
      t98 = (t96 + t95) * t65
      t99 = (t34 + t91 + t92 + t93 + t12 + t13 + t14 + t32 - t33) * t1
      t100 = (t77 - t91 - t92 - t93 - t12 - t13 - t14 + t32 - t34 + t75)
     # * t17 * t87
      t101 = t59 * t1 * (t90 + t74 + t76) - t100 * t80 - t97 - t98 - t99
      t78 = 4 * t68 * t73 + t78
      t78 = sqrt(t78)
      t102 = -t76 + t77 + t75 - t74 - t33 + t34 + t32 + t78
      t102 = 0.1q1 / t102
      t103 = t59 * t68
      t83 = t103 * t102 + t83
      t79 = -t76 + t77 + t75 - t74 - t33 + t34 + t32 - t79
      t104 = t81 * (t80 - t79)
      t105 = t77 - t91 - t92 - t93 - t12 - t13 - t14 + t32 - t34 + t75
      t106 = t80 ** 2
      t107 = t106 ** 2
      t108 = t80 * t106
      t109 = t87 ** 2
      t110 = t109 ** 2
      t111 = t87 * t110
      t112 = t87 * t109
      t113 = t109 * t106
      t114 = t113 * t17
      t115 = -t96 - t95
      t116 = -t48 - t94
      t117 = t88 * t1
      t118 = t117 * t109
      t119 = t76 + t74
      t120 = (0.1q1 / 0.4q1)
      t121 = -t120 * t114 * t105 + t39 * t80 * (t81 * (t1 * (-t34 - t91 
     #- t92 - t93 - t12 - t13 - t14 - t32 + t33) + t115 * t65 + t116 * t
     #4) + t118 * t80) + t1 * (t119 * t87 * t80 - t60 - t66 - t67)
      t122 = t79 * t81
      t123 = t122 * t39
      t124 = t123 * t17 - t1
      t125 = t123 * t26 + t25
      t78 = -t76 + t77 + t75 - t74 - t33 + t34 + t32 - t78
      t78 = 0.1q1 / t78
      t103 = t103 * t78 + t123
      t123 = t123 * t10 + t22
      t126 = t68 * t1
      t127 = t35 * t37
      t128 = t127 * t39 - t1
      t129 = t44 * t10 + t22
      t130 = 0.1q1 / t35
      t61 = t61 * t130 + t44
      t41 = t41 * t26 + t25
      t131 = t40 * t37
      t132 = -t131 + t122
      t133 = -t131 + t82
      t44 = t44 * t26 + t25
      t134 = t43 * t37
      t122 = -t134 + t122
      t82 = -t134 + t82
      t89 = t89 * t79
      t100 = -t100 * t79 + t59 * t1 * (t89 + t74 + t76) - t97 - t98 - t9
     #9
      t135 = t96 + t95
      t136 = t48 + t94
      t137 = t79 ** 2
      t138 = t137 ** 2
      t139 = t109 * t137
      t140 = t139 * t17
      t105 = -t120 * t140 * t105 - t39 * t79 * (t81 * ((t34 + t91 + t92 
     #+ t93 + t12 + t13 + t14 + t32 - t33) * t1 + t135 * t65 + t136 * t4
     #) - t118 * t79) - t1 * ((-t76 - t74) * t87 * t79 + t66 + t67 + t60
     #)
      t118 = -t4 * t21
      t141 = t18 * t1
      t142 = t141 * t59 - t20
      t143 = t77 + t32 + t75
      t144 = t37 ** 2
      t145 = t37 * t144
      t146 = t36 ** 2
      t147 = t146 ** 2
      t148 = t36 * t146
      t149 = t38 ** 2
      t150 = t32 * t149
      t151 = t150 * t145
      t152 = -t71 - t70
      t153 = t77 + t75
      t154 = t13 + t93
      t155 = t152 * t31
      t156 = t1 * t38
      t157 = t31 * t37
      t158 = t32 * t63
      t159 = t14 * t63
      t160 = t38 * t144
      t161 = t160 * t148
      t162 = t1 ** 2
      t163 = t162 ** 2
      t164 = t1 * t162
      t165 = t31 ** 2
      t166 = t162 * t149
      t167 = t37 * (t76 - t14 - t32 + t74) - t1
      t168 = (-t60 - t66 - t67) * t1
      t169 = t40 * t1
      t170 = t160 * t14
      t171 = t55 * t31
      t172 = t171 - t66 - t67
      t173 = -t12 - t91 - t92
      t174 = t67 + t66
      t175 = -t55 + t1
      t176 = t16 ** 2
      t177 = t32 * t131
      t178 = t5 ** 2 * t176
      t179 = t63 * t146
      t180 = t77 + t91 + t92 + t12 + t75
      t181 = t77 - t74 + t75 - t76
      t182 = t165 * t37
      t183 = t181 * t1
      t184 = t166 * t37
      t185 = t146 * (t183 * t14 * t144 * t149 - t40 * t31 * t145 * t153 
     #+ t149 * t164 + t184 * t180 + t150 * t1 * t144 * (t131 * t152 + t1
     #2 + t13 + t75 + t77 + t91 + t92 + t93) + t37 * (t174 * t38 - t182)
     # * t1) + t179 * (t37 * (-t166 * t36 * t88 - t174 * t36 + t60 * (t7
     #7 - t36 + t74 + t75 + t76) + t169 * t173) + t177 * (-t37 * (t76 + 
     #t36 + t14 + t74) + t1) + t14 * t172 * t37 + t178 * t63 * t1 + t93 
     #* t156 * t175 + t13 * t156 * t175)
      t186 = t55 - t1
      t187 = t32 * t37
      t188 = t1 * t165
      t189 = t63 * t162
      t190 = t38 * t162
      t191 = t190 * t25
      t192 = t37 * t1
      t193 = t14 * t55
      t194 = t31 * t36
      t195 = t131 * t18
      t196 = t14 * t4
      t197 = t32 * t38
      t189 = t55 * (t1 * (t25 * (-t31 * t4 + t38 * (t141 * t40 - t196)) 
     #+ t66 * (-t14 * t38 - t31)) + t197 * (t36 * (-t159 * t1 + t194 * t
     #144 - t190) + t1 * (t195 - t4) * t25 + t1 * (t131 * t69 - t65) * t
     #22)) + t36 * (t4 * (t31 * (t25 * t36 * t144 - t189 * t26) - t191) 
     #+ t65 * (t31 * (t22 * t36 * t144 - t189 * t10) - t190 * t22) + t92
     # * t63 * t60 * t186 + t91 * t63 * t60 * t186 + t12 * t63 * t60 * t
     #186 + t93 * t63 * t60 * t186 + t13 * t63 * t60 * t186 + t192 * (t3
     #6 * (t1 * (t157 * t88 + t74) * t149 + t31 * (-t187 + t1) * t38) + 
     #t188) + t193 * t1 * t149 * t186)
      t198 = t13 + t91 + t92 + t93 + t12
      t199 = t156 + t157
      t199 = t198 * t37 - t1 + t37 * (t199 * t26 + t25) * t18 + t37 * (t
     #10 * t199 + t22) * t69 + t170 * t88
      t200 = t144 * t149
      t201 = t63 * t88 + 1
      t202 = t67 + t66
      t203 = t160 * t32
      t204 = t203 * t146
      t205 = (0.3q1 / 0.4q1)
      t206 = (0.1q1 / 0.8q1)
      t207 = (0.1q1 / 0.16q2)
      t208 = (0.1q1 / 0.32q2)
      t209 = (0.3q1 / 0.8q1)
      t159 = t209 * t151 * t148 * t143 + t120 * t185 - t39 * t189 - t206
     # * (t37 * t148 * (t14 * t117 * t144 * t38 * t149 - t144 * t165 - t
     #160 * (t76 + t36 + t91 + t92 + t12 + t74) * t31 + t166) + t161 * (
     #t159 * (t77 - t36 - t74 + t75 - t76) + t156 * (t155 * t37 - t74 - 
     #t75 - t76 - t77) - t157 * t154 + t131 * t153 + t158 * (t131 * t88 
     #+ t12 + t13 + t91 + t92 + t93))) + t207 * t200 * t147 * t199 + t20
     #8 * t36 * t147 * t149 * t145 * t201 + t205 * t204 * t202 - t60 * (
     #t170 * t146 + t169 * t167 + t168)
      t185 = t35 ** 2
      t189 = t185 ** 2
      t210 = t35 * t185
      t211 = t117 * t38
      t212 = t144 * t31
      t213 = t63 * t210
      t214 = -t77 - t74 - t75 - t76
      t215 = t1 * t149
      t216 = t144 * t210
      t217 = t43 * t1
      t218 = -t76 + t14 - t74
      t219 = -t127 + t1
      t220 = t212 * t35
      t221 = t32 * t134
      t222 = t166 * t127
      t223 = -t77 + t74 - t75 + t76
      t224 = t77 + t91 + t92 + t93 + t12 + t13 + t75
      t225 = t224 * t1
      t226 = t35 * t174
      t227 = t14 * (t156 * t223 + t66 + t67)
      t228 = t160 * t185
      t229 = t127 - t1
      t150 = t185 * (t38 * (-t164 * t38 + t37 * (-t1 * t174 - t220 * t14
     #)) + t192 * t92 * t149 * t229 + t192 * t91 * t149 * t229 + t192 * 
     #t12 * t149 * t229 + t192 * t93 * t149 * t229 + t192 * t13 * t149 *
     # t229 + t150 * t35 * t144 * t229 + t188 * t144 - t200 * t178 * t1 
     #+ t63 * (t35 * (t184 * t18 + t212 * t4) - t190 * t4) * t26 + t63 *
     # (t35 * (t184 * t69 + t212 * t65) - t190 * t65) * t10) + t228 * (t
     #197 * (t127 * (t211 + t14 + t74 + t76) - t225) + t226 + t227 + t60
     # * (-t77 + t35 - t74 - t75 - t76))
      t184 = t203 * t185
      t151 = t151 * t210 * t143
      t230 = t35 * t189 * t149 * t145 * t201
      t119 = -t120 * t150 + t205 * t184 * t202 - t206 * (t216 * (t215 * 
     #t214 - t182 - t63 * (t76 + t92 + t93 + t13 + t74) * t31) + t213 * 
     #(t203 * (t134 * t88 + t12 + t13 + t91 + t92 + t93) + t190 + t144 *
     # (t153 * t38 - t31) * t35 - t212 * (t12 + t91) - t160 * t60 * t88 
     #+ t170 * (t211 - t35 - t74 + t75 - t76 + t77))) + t207 * t200 * t1
     #89 * t199 + t208 * t230 + t209 * t151 + t39 * (t35 * (t215 * t14 *
     # t127 * t219 + t221 * (-t220 + t190) - t182 * t162 + t191 * t4 - t
     #222 * t119 + t157 * (t66 * t219 + t67 * t219)) + t217 * (t31 * (t3
     #5 * (-t211 - t91 - t92 - t93 - t12 - t13) * t144 + t192 * (t77 - t
     #35 + t91 + t92 + t93 + t12 + t13 + t75)) + t66 * t1 + t14 * t174 *
     # t37 + t187 * (t127 * (t218 * t38 + t31) + t66 + t67))) - t60 * (t
     #217 * t167 + t170 * t185 + t168)
      t150 = t14 * t34
      t151 = t33 * t1
      t153 = t17 * (t32 * t33 + t150) * t112
      t167 = t81 * t1
      t187 = t25 * t30
      t191 = t26 * t31
      t192 = t191 + t187
      t211 = t17 ** 2
      t215 = t10 * t31
      t220 = t31 * t174 * t211
      t230 = t162 * t31
      t81 = t230 * t81
      t231 = (t76 - t14 - t32 - t34 + t74) * t1 * t87
      t232 = t17 * (t76 + t91 + t92 + t93 + t12 + t13 + t74) + t117
      t233 = t30 ** 2
      t234 = t33 * t88
      t235 = t14 * t88
      t236 = t1 * t17
      t237 = t236 * t233
      t238 = t232 * t30
      t239 = t238 + t234 + t235 - t237
      t240 = t110 * t138
      t241 = -t67 - t66
      t242 = t34 * t17
      t115 = t115 * t69 + t116 * t18 + t242
      t116 = t174 * t30
      t174 = t14 * (t17 * t241 + t183)
      t183 = t22 * t30
      t243 = t211 * (t183 + t215) * t65
      t244 = t30 * t88 * t162
      t245 = t14 * t17
      t224 = t34 * t224
      t246 = t211 * t165
      t247 = t162 * t233
      t248 = (-t34 * t211 * t31 - t243 - t244 + t245 * (t33 - t32) - t21
     #1 * t4 * t192) * t87
      t249 = t139 * (t1 * (t33 * (t76 + t74 + t75) + t224 - t246 + t247 
     #+ t32 * (t32 + t92)) - t32 * t139 * t211 * t30 + t248 * t79) + t13
     #9 * (t236 * (-t198 * t30 * t87 * t79 + t77 * t31 + t116) + t174 + 
     #t32 * (t1 * (t77 + t91 + t93 + t12 + t13 + t75) + t115 * t87 * t79
     #))
      t250 = t87 * t79
      t251 = t250 * t30
      t155 = (t214 * t30 + t155) * t1
      t117 = t14 * (t17 * t181 + t117)
      t252 = t32 * t17
      t253 = t112 * t79 * t137
      t254 = -t76 - t91 - t92 - t93 - t12 - t13 - t74
      t255 = t77 - t14 + t75
      t256 = t254 * t31
      t257 = t253 * t211
      t258 = t32 * t140
      t259 = t72 + t70 + t71
      t89 = t120 * t249 + t205 * t258 * t202 - t206 * (t257 * (t251 * t2
     #55 + t256) + t253 * (t17 * (-t31 * (t251 + t31) * t211 + t155 + t2
     #47) + t117 + t252 * (t89 + t91 + t92 + t93 + t12 + t13))) + t207 *
     # t240 * t17 * t239 + t208 * t111 * t79 * t138 * t30 * t211 * t259 
     #+ t209 * t252 * t253 * t143 + t39 * (t79 * (t81 * (t77 + t91 + t93
     #) + (t162 * (t30 * (-t74 + t14 + t32) - t215 * t69) - t220 - t33 *
     # t162 * t30 - t162 * t192 * t18) * t109 * t79) + t79 * (-t153 * t1
     #37 + t167 * (t31 * (t1 * (t75 + t92 + t12 + t13) + t66 * t17 - t15
     #1) + t67 * (t34 + t14 + t32 + t33) + t66 * (t34 + t14 + t32)) + t1
     # * (t32 * t218 + t33 * (t32 - t91 - t92 - t93 - t12 - t13)) * t109
     # * t79)) - t60 * (t14 * t140 + t231 * t79 + t168)
      t110 = t110 * t107
      t137 = t77 + t74 + t75 + t76
      t116 = t137 * t31 + t116
      t138 = -t13 - t91 - t92 - t93 - t12
      t192 = -t34 - t14 - t32 - t33
      t218 = -t77 - t91 - t92 - t12 - t75
      t249 = t87 * t80
      t260 = t249 * t30
      t261 = t112 * t108
      t262 = t261 * t211
      t90 = t262 * (t260 * t255 + t256) + t261 * (t17 * (-t31 * (t260 + 
     #t31) * t211 + t155 + t247) + t117 + t252 * (t90 + t91 + t92 + t93 
     #+ t12 + t13))
      t155 = t32 * t114
      t81 = t207 * t110 * t17 * t239 + t120 * (t113 * (t1 * (t32 * (t77 
     #+ t91 + t12 + t32 + t75) + t224 - t246 + t247) - t32 * t113 * t211
     # * t30 + t248 * t80) + t113 * (t1 * (t116 * t17 + t32 * (t13 + t92
     # + t93)) + t174 + (t32 * t115 + t242 * t138) * t87 * t80)) + t205 
     #* t155 * t202 - t206 * t90 + t208 * t111 * t80 * t107 * t30 * t211
     # * t259 + t209 * t261 * t252 * t143 - t39 * (t80 * (t81 * t218 + (
     #t162 * (t30 * (t76 - t32 + t74) + t31 * t88) + t220) * t109 * t80)
     # + t80 * (t1 * (t32 * (t76 - t14 + t74) - t150 + t33 * (t34 + t91 
     #+ t92 + t93 + t12 + t13 - t32)) * t109 * t80 + t167 * (t66 * t192 
     #+ t67 * t192 + t60 * (t33 - t93 - t13)) + t153 * t106)) - t60 * (t
     #14 * t114 + t231 * t80 + t168)
      t90 = t19 ** 2
      t106 = t4 ** 2
      t107 = t211 * t4 * t106 * t30
      t23 = t59 * t236 * t106 * t19 * (t23 * t30 + t31) + t4 * (-t107 * 
     #t19 * t90 - t230 - t23 * t1 * (t34 + t32) + t17 * (-t33 + t32) * t
     #90 * t106)
      t111 = t43 + t31
      t115 = t141 * t38
      t153 = t115 + t4
      t168 = t4 * t1
      t192 = t14 * t37
      t208 = t192 + t1
      t209 = t208 * t38 + t157
      t220 = t185 * t37
      t224 = t158 * t18 + t4
      t231 = -(0.3q1 / 0.2q1)
      t248 = 3 * t32 * t4
      t256 = t134 * t18
      t208 = t156 * t4 * t208
      t158 = t158 * t1
      t260 = t157 * t4
      t115 = t192 * (t115 + t4)
      t263 = t63 * t185
      t264 = t32 * t205 * t4
      t219 = -t120 * t263 * (t229 * t4 + t115 + t37 * (t221 - t60) * t18
     #) + t206 * t160 * t18 * t210 * t209 + t39 * t35 * (-t222 * t18 + t
     #208 + t158 * (-t256 + t4) + t260 * t219) - t230 * (t256 - t4) + t2
     #28 * (t263 * t207 * t18 + t264)
      t222 = t14 * t1
      t228 = t40 + t31
      t229 = t146 * t37
      t256 = t160 * t146
      t115 = -t120 * t179 * (t186 * t4 + t115 + t37 * (t177 - t60) * t18
     #) + t206 * t161 * t18 * t209 + t39 * t36 * (-t166 * t55 * t18 + t2
     #08 + t158 * (-t195 + t4) + t260 * t175) - t230 * (t195 - t4) + t25
     #6 * (t179 * t207 * t18 + t264)
      t158 = t156 * t152
      t175 = t13 + t91 + t92 + t93 + t12
      t186 = t157 * t10 + t22
      t195 = t157 * t26 + t25
      t207 = t127 + t1
      t208 = t25 * t31
      t260 = t22 * t31
      t263 = (t76 - t14 + t74) * t1
      t264 = t164 * t31 * t30
      t265 = t14 * t162
      t266 = t127 * t60
      t166 = t1 * (t198 * t38 + t31) + t66 + t67 + t166 * t88 + t157 * t
     #255 + t197 * (t37 * (t156 * t26 + t25) * t18 - t1 + t192 + t37 * (
     #t156 * t10 + t22) * t69)
      t192 = t77 - t14 + t75
      t70 = t71 + t70
      t71 = t202 * t1
      t255 = t188 * t17
      t267 = t43 + t31
      t134 = t127 * (t1 * (-t31 * (t77 + t76) - t75 * t197) + t227) + t1
     # * (t1 * (t241 * t30 + t31 * (-t77 - t91 - t75)) - t226 - t43 * t1
     #62 * t30 + t127 * (-t31 * (t75 + t74) - t178 * t38) + t14 * t241 -
     # t92 * t1 * t267 - t12 * t1 * t267 - t93 * t1 * t267 - t13 * t1 * 
     #t267) + t31 * (-t236 * t202 - t216) + t255 * t207 + t32 * (t149 * 
     #t70 * t145 * t210 - t71 - t134 * (t77 + t91 + t92 + t93 + t12 + t1
     #3) * t1) + t43 * (-t162 * (t77 + t75 + t91) + t192 * t144 * t185)
      t144 = (0.5q1 / 0.16q2)
      t145 = (0.9q1 / 0.4q1)
      t197 = t14 * t30
      t183 = (t183 + t215) * t69
      t210 = t31 * t202 * t211
      t216 = t151 * (t34 + t13 + t12 + t92 + t93 + t91)
      t226 = t32 * (t76 + t74 - t33 - t34 - t14) * t1
      t187 = t162 * (t191 + t187) * t18
      t191 = t230 * (t76 - t14 - t32 - t34 + t74)
      t227 = t34 * t137
      t241 = t60 * t88
      t135 = t211 * (t31 * (t77 - t14) + t67 * t30) + t243 + t244 + t242
     # * (t33 + t91 + t92 + t93 + t12 + t13) + t32 * (t135 * t69 + t136 
     #* t18 + t17 * (-t34 + t14))
      t136 = t77 + t91 + t92 + t93 + t12 + t13 + t75
      t243 = t202 * t17
      t244 = t32 * t202
      t267 = t253 * t17
      t268 = t250 * t17
      t269 = t55 * t202
      t270 = t13 + t93
      t271 = t202 * t30
      t137 = t256 * (t32 * t138 - t196 * t26 + t241) + t146 * (t170 * (t
     #158 + t74 - t75 + t76) - t190 + t212 * (t76 + t91 + t92 + t93 + t1
     #2 + t13 + t33 + t74) + t63 * t137 * t1)
      t40 = t162 * (t18 * (t40 * t195 + t208) + t69 * (t40 * t186 + t260
     #) + t194) - t264 + t32 * (t40 * (t37 * (t263 + t229) - t162) - t60
     # * (t55 + t1)) - t265 * t228 + t269 * t33 + t60 * t55 * t175
      t138 = -t33 - t14 - t32
      t146 = t30 * t211
      t194 = t261 * t17
      t196 = t152 * t1
      t261 = t249 * t17
      t272 = t261 * (t249 * t150 + t32 * (t249 * t33 - t66 - t67))
      t273 = 4 * t249 * t151 * t14
      t274 = t76 + t91 + t92 + t93 + t12 + t13 + t14 + t33 + t74
      t167 = t167 * (t77 + t91 + t92 + t93 + t12 + t13 + t14 + t32 + t34
     # + t75)
      t275 = (t143 * t17 + t196) * t109
      t276 = (-t76 - t74) * t1
      t277 = t243 + t276
      t278 = t120 * t140 * t274 + t206 * t267 * t259 - t39 * t79 * (t275
     # * t79 + t167) - t250 * t277 - t126
      t279 = t76 + t74 + t14 + t13 + t12 + t92 + t93 + t91
      t280 = (-t77 - t75 - t34 - t14 - t32 - t13 - t12 - t92 - t93 - t91
     #) * t1
      t281 = t70 * t1
      t282 = (t77 + t32 + t75) * t17
      t120 = t120 * t114 * t274 + t206 * t194 * t259 - t39 * t80 * (t275
     # * t80 + t167) - t126 + t249 * (-t243 - t276)
      t62 = 0.1q1 / t62
      t5 = 0.1q1 / t5
      t9 = 0.1q1 / t9
      t61 = 0.1q1 / t61
      t64 = 0.1q1 / t64
      t15 = 0.1q1 / t15
      t11 = 0.1q1 / t11
      t167 = 0.1q1 / t6
      t206 = -0.1q1 / t45
      t274 = 0.1q1 / t29
      t275 = 0.1q1 / t79
      t283 = 0.1q1 / t80
      t284 = -0.1q1 / t53
      t285 = -0.1q1 / t54
      t286 = 0.1q1 / t25
      t45 = 0.1q1 / t45
      t287 = 0.1q1 / t7
      t288 = 0.1q1 / t4
      t8 = 0.1q1 / t8
      t289 = -0.1q1 / t42
      t42 = 0.1q1 / t42
      t290 = t56 ** 2
      t291 = t290 ** 2
      t56 = t56 * t290
      t292 = t274 ** 2
      t293 = t45 ** 2
      t294 = t42 ** 2
      t295 = t21 ** 2
      t296 = t21 * t295
      t297 = t3 ** 2
      t298 = t128 ** 2
      t299 = t298 ** 2
      t128 = t128 * t298
      t300 = t288 ** 2
      t301 = t73 * t283
      t302 = t3 * t11
      t303 = t22 * t11
      t304 = t3 * t15
      t305 = t3 * t25
      t306 = t22 * t47
      t307 = t8 * t3
      t308 = t129 * t11
      t309 = t56 * t57 * t62
      t310 = t21 * t142
      t311 = t16 * t167
      t312 = t24 * t274
      t313 = t24 * t23
      t314 = t313 * t288
      t315 = t63 * t19
      t316 = t315 * t5 * t16
      t317 = t24 * t11
      t318 = t28 * t317 * t21 * t118
      t319 = t28 * t8 * t288
      t320 = t319 * t11 * t297
      t321 = t320 * (t295 * t5 * t167 + (t287 * (t118 * (t10 * t21 + t17
     # * t24) + t310 * t24) - t311 * t21 * t118) * t288 * t15) * t274
      t103 = 0.1q1 / t103
      t104 = 0.1q1 / t104
      t83 = 0.1q1 / t83
      t322 = 0.1q1 / t16
      t323 = t22 * t26
      t324 = t323 + t51
      t325 = t84 ** 2
      t326 = t84 * t325
      t327 = t124 ** 2
      t328 = t124 * t327
      t329 = t24 * t287
      t330 = t10 * t287
      t331 = t327 * t78 * t103
      t332 = t325 * t102 * t83
      t333 = t2 * t16
      t334 = t287 * t324
      t21 = t21 * t274
      t335 = t307 * t5
      t336 = t167 * t9
      t310 = t336 * (t316 * t60 * t2 * t56 * t287 * t58 * t289 * t62 * t
     #64 + t320 * t19 * t274 * t15 * t16 * (t118 * (-t21 + t17) + t310) 
     #+ t335 * (t11 * ((t330 * t300 * t27 * t274 - t329 * t300 * t27 * t
     #292) * t118 * t25 * t295 + t333 * t68 * t87 * t104 * (t332 - t331)
     #) + t334 * t301 * t322 * t275 * t47 * t162))
      t23 = t310 + t9 * (t1 * (t72 * t288 * t287 * t5 * t162 * (t167 * (
     #t306 - t305) + t303 * (t306 * t286 - t3)) * t130 * t58 + t307 * (-
     #t304 * t303 * t142 * t300 * t287 + (t301 * t2 * t25 * t275 * t287 
     #- t302 * t288) * t167 * t5 * t1)) + t316 * (-t305 * t296 * t167 * 
     #t11 * t274 * t27 * t42 * t45 + t60 * (t309 * t11 * t58 * t284 * t2
     #89 * t27 + (-t308 * t285 * t27 - t167) * t130 * t206 * t61 * t128)
     # * t287 * t64 * t2) - t318 * t297 * t300 * t287 * t292 * t8 * t15 
     #+ t28 * t311 * t200 * t294 * t288 * t287 * t293 * t274 * t11 * t5 
     #* t296 * (t19 * (t23 * (t312 - t10) - t24 * (-t59 * (t107 * t90 + 
     #t230 * t18) + t4 * (t4 * (-t31 * t19 * t211 + t245 * t19) - t222) 
     #- 3 * t168 * (t34 + t32 - t33) + 5 * t242 * t106 * t19)) + t314) +
     # t321)
      t82 = 0.1q1 / t82
      t106 = 0.1q1 / t122
      t107 = 0.1q1 / t132
      t122 = 0.1q1 / t133
      t132 = t326 * t86
      t133 = t132 * t102
      t142 = t133 * t122 * t82 * t83
      t242 = t125 * t103
      t245 = t242 * t328 * t78
      t310 = -t245 * t107 * t106 + t142
      t316 = t5 * t1
      t320 = t316 * t167
      t321 = t283 ** 2
      t337 = t73 ** 2
      t338 = t73 * t337
      t339 = t275 ** 2
      t340 = t275 * t339
      t341 = t125 ** 2
      t6 = t2 * t6
      t342 = t303 * t307
      t343 = t48 * t162
      t344 = t302 * t19
      t345 = t130 * t128 * t61
      t346 = t342 * t1 * t25 * t287 * (t304 - t320)
      t347 = t311 * t63
      t7 = t7 * t47 * t167
      t348 = t84 * t83
      t349 = t348 * t102
      t350 = t73 * t5 * t167
      t351 = t331 * t125
      t352 = t351 * t123
      t353 = t332 * t85
      t354 = t353 * t86
      t7 = t9 * (t1 * (t337 * (t346 * t321 * t275 + t346 * t283 * t339) 
     #+ t347 * t64 * t5 * t31 * ((t2 * t41 * t287 * t107 * t87 * t122 + 
     #t344 * t289) * t58 * t62 * t56 - t345 * (t2 * t44 * t287 * t87 * t
     #106 * t82 + t344 * t206)) + t301 * (t342 * (-t316 * t2 * t287 + t3
     #04 * (t287 * (t6 - t48) - t47)) + t343 * t30 * t5 * t167 * t58 * t
     #130 * t287 * (-t306 + t305)) * t275) + (t347 * t2 * t5 * t287 * t8
     #7 * t310 + ((t349 * (t7 - t2) + (-t7 + t2) * t78 * t103 * t124) * 
     #t87 * t11 * t15 * t16 + t350 * t287 * t322 * (t332 * t86 ** 2 * t3
     #21 - t331 * t341 * t339)) * t8 * t297 + t335 * t47 * t73 * t167 * 
     #t287 * t322 * (-t354 * t321 + t352 * t339)) * t104 * t68)
      t305 = t275 + t283
      t316 = t130 ** 2
      t344 = t130 * t316
      t346 = t58 ** 2
      t355 = t58 * t346
      t356 = t63 * t60
      t357 = t162 * t22
      t312 = t312 * t28 * t18
      t358 = t306 * t162
      t359 = t358 * t25 * t337
      t360 = t300 * t11
      t361 = t301 * t1
      t362 = t311 * t25
      t363 = t25 * t167
      t364 = t11 * t15
      t365 = t335 * t287
      t366 = t365 * t167
      t367 = t11 * t5
      t368 = t3 * t41
      t369 = t47 * t57
      t332 = t332 * t86 * t283
      t351 = t351 * t275
      t370 = t60 * t5
      t371 = t307 * t287
      t308 = t371 * (t318 * t304 * t18 * t288 * t300 * t274 + (-t359 * t
     #339 * t283 * t322 + (t351 - t332) * t104 * t68 * t2) * t167 * t5) 
     #+ t370 * (-t200 * t2 * t176 * t128 * t82 * t61 * t11 * t167 * t106
     # * t87 + (t316 * t206 * t61 * t128 * (t167 * (t129 * t47 - t3 * t4
     #4) - t308 * t3) + (t167 * (-t369 + t368) + t302 * t57) * t289 * t3
     #46 * t62 * t56) * t287 * t19) * t64
      t21 = t9 * t308 + t9 * (t361 * t275 * (t364 * t1 * (t362 + t334) +
     # t363 * t287 * (t25 * (t1 * t73 * t305 - t17) + t94) * t322 * t5) 
     #* t8 * t297 + t366 * (-t359 * t275 * t321 * t322 + t360 * (-t21 * 
     #t49 * t24 * t118 + t357 * (t20 - t141) + t312 * (t118 * t288 - t1)
     # * t295) + t361 * (t306 * t48 * t322 + (-t323 - t51) * t11 * t162)
     # * t275) + t367 * (t63 * t87 * t2 * (t356 * t56 * t62 * t64 * t107
     # * t122 + (t80 * t326 * t102 * t82 * t83 * t122 - t79 * t328 * t10
     #3 * t78 * t106 * t107) * t87 * t104 * t68) * t167 * t176 + t311 * 
     #t288 * t287 * (-t72 * t163 * t22 * t58 * t130 - t313 * t28 * t200 
     #* t296 * t19 * t274 * t294 * t293 * (t45 + t42)) + t60 * t47 * t19
     # * t287 * t27 * t64 * (t128 * t129 ** 2 * t316 * t206 * t285 * t61
     # - t56 * t57 ** 2 * t346 * t289 * t284 * t62)))
      t49 = t26 * t85
      t51 = t10 * t86 + t49
      t72 = t107 ** 2
      t79 = t104 ** 2
      t80 = t104 * t79
      t94 = t106 ** 2
      t118 = t82 ** 2
      t176 = t122 ** 2
      t284 = t122 * t176
      t285 = t62 ** 2
      t308 = t64 ** 2
      t313 = t83 ** 2
      t318 = t57 * t41
      t323 = t318 * t159 * t346 * t72 * t176
      t40 = -t323 * t285 + (t159 * (t41 * t10 + t26 * t57) + t318 * (t14
     #4 * t160 * t147 * t201 + t39 * t63 * t148 * t199 + t178 * t131 * t
     #1 - t1 * (-t190 * t30 - t66 - t67 + t157 * (-t77 - t74 - t75 - t76
     #)) * t36 - t1 * (t182 * t17 + t156 * (-t77 - t75 - t13 - t93)) * t
     #36 + t91 * t162 * t228 + t12 * t162 * t228 + t92 * t162 * t228 + t
     #215 * t65 * t162 - t14 * t131 * (t1 * (-t75 + t74) - t229) + t177 
     #* t1 * t270 + t231 * t229 * t166 + t145 * t204 * t143 + t212 * t14
     #8 - 4 * t193 * t60 - t59 * t40 + t244 * t1 - t151 * (t60 - t66 - t
     #67) - t230 * (-t77 - t13 - t93) - t14 * (-t71 + t269) - t161 * (t7
     #7 + t75) + t271 * t162 - t131 * (t222 * (-t77 + t76) + t32 * (-t1 
     #* t180 + t256 * t88)) - 3 * t55 * (t169 * t14 + t32 * t172) + t205
     # * t137)) * t176 * t72 * t346 * t62
      t65 = t129 * t61
      t71 = t129 * t44
      t88 = t316 * t94 * t61 * t118
      t131 = t88 * t128
      t137 = t85 * t86
      t70 = t137 * (t144 * t146 * t110 * t259 + t145 * t155 * t143 + t20
     #5 * (t113 * (-t252 * t154 + t14 * (t17 * t223 + t196)) + t114 * (t
     #32 * t173 + t232 * t31 + t227 + t246 - t247)) + t231 * t113 * t135
     # + t39 * t108 * t17 * (t234 * t112 + t235 * t112 - t237 * t112 + t
     #238 * t112) + t59 * (-t191 - t262 * t32 * t30 + t249 * (t162 * (-t
     #183 + t197) - t187 - t210 - t216 - t226)) - t249 * t1 * (t33 * t21
     #4 + t32 * t218) - t162 * (t31 * (-t75 - t12 - t92 - t91) - t271) -
     # t1 * (t67 * t138 + t66 * t138 + t60 * (-t77 + t33 - t13 - t93)) -
     # t30 * (-t249 * t136 * t162 + t262 * t192) + t249 * t164 * t233 - 
     #t249 * (t1 * (t17 * (t165 * t17 - t271) - t178 - t32 * t270) + t14
     # * (t1 * (-t77 + t74 - t75 + t76) + t243)) - t194 * (-t146 * t31 +
     # t32 * t70) - 3 * t272 - t273)
      t108 = t10 * t125
      t110 = t125 * t123
      t112 = t72 * t94
      t113 = t112 * t275
      t75 = t113 * t78 * t103 * t328 * (t89 * (t123 * (-t242 + t26) + t1
     #08) + t110 * (t144 * t240 * t30 * t211 * t259 + t145 * t258 * t143
     # - t205 * (t140 * (t33 * t254 + t247) + t139 * (t17 * (t32 * t198 
     #- t227 - t241 - t246) + t117)) + t231 * t139 * t135 + t39 * t267 *
     # t239 - t59 * (t191 + t257 * t32 * t30 + t250 * (t162 * (t183 - t1
     #97) + t210 + t216 + t226 + t187)) + t251 * (t162 * (t77 + t91 + t9
     #2 + t93 + t12 + t13) + t139 * (-t77 + t14 - t75) * t211) + t1 * (t
     #31 * (t1 * t136 + t243) - t255 + t34 * t202 + t14 * t202 + t244) +
     # t250 * (t236 * t116 + t174 + t32 * (t140 * t152 + t225)) + t253 *
     # t17 * t211 * t31 * t30 + t250 * t1 * (-t246 + t178 + t34 * (t75 +
     # t34)) - 3 * t268 * (t250 * t150 + t32 * (t250 * t33 - t66 - t67))
     # - 4 * t250 * t151 * t14))
      t116 = t200 * t79 * t68
      t12 = t356 * t308 * (t40 * t56 + t131 * (t119 * (t129 * t26 + t44 
     #* (-t65 + t10)) + t71 * (-t205 * (t220 * (t1 * (t235 * t37 * t149 
     #+ t214 * t38) - t157 * (t33 + t93)) + t185 * (t203 * t198 + t190 +
     # t212 * (t158 - t91 - t92 - t12 - t13 - t74 - t76) + t170 * t181))
     # + t231 * t220 * t166 + t39 * t213 * t199 - t59 * (t162 * (t18 * (
     #t43 * t195 + t208) + t31 * t35 + t69 * (t43 * t186 + t260)) - t264
     # + t32 * (t43 * (t37 * (t263 + t220) - t162) - t60 * t207) - t265 
     #* t111 + t127 * t33 * t202 + t266 * t175) + t145 * t184 * t143 + t
     #144 * t160 * t189 * t201 - t134 - 3 * t127 * (t217 * t14 + t32 * (
     #t127 * t31 - t66 - t67)) - 4 * t266 * t14))) * t109
      t13 = t367 * t311 * t287 * t9
      t38 = t13 * t87 * t37 * (t63 * t68 * t79 * (t354 * t81 * t283 * t1
     #76 * t118 + t113 * t352 * t89) + t60 * (t71 * t88 * t298 * t119 + 
     #t323 * t290 * t62) * t87 * t308)
      t40 = t283 * t83
      t43 = t113 * t245 * t123 * t89
      t69 = t85 * t83
      t74 = t69 * t133 * t81 * t118
      t88 = t367 * t347 * t87
      t40 = t88 * t287 * t9 * (t68 * (t63 * (t74 * t283 * t176 - t43) * 
     #t80 + t63 * (t133 * (t40 * t176 * t82 * t118 + t40 * t284 * t118) 
     #* t81 * t85 + t43 * (t106 + t107)) * t79) + t60 * t87 * t308 * (t3
     #09 * (-t107 * t176 * t72 + t72 * (t176 * t64 - t284)) * t346 * t15
     #9 * t41 + t71 * t131 * t119 * (-t82 - t64 - t106)))
      t43 = 0.1q1 / t53
      t53 = 0.1q1 / t54
      t54 = t363 + t303
      t91 = t27 ** 2
      t92 = t50 ** 2
      t93 = t286 ** 2
      t29 = 0.1q1 / t29 ** 2
      t113 = t24 ** 2
      t117 = t52 ** 2
      t131 = t303 * t286
      t133 = t47 * t11
      t134 = t162 * t18 * t300
      t135 = t287 * t9
      t136 = t135 * t5
      t138 = t162 * t297
      t46 = t136 * t288 * (t138 * t322 * t8 * (t10 * t11 + t167 * t26) +
     # (-t133 * t63 * t113 * t45 * t42 * t91 * t292 + t63 * t42 * t45 * 
     #(t167 * (-t24 * t47 + t3 * t46) + t317 * t3) * t27 * t274) * t25 *
     # t296 - t319 * t297 * t18 * t322 * (t167 * t46 + t317) * t274 * t2
     #95)
      t27 = t46 + t136 * (t134 * t322 * t54 * t8 * t297 + t307 * (t2 * (
     #(t167 * t288 * t27 * t274 + t317 * t288 * t91 * t292) * t25 * t295
     # + (-t131 - t167) * t288 * t162 - t4 * t92 * t52 * t11 * t29 * t28
     #6 * t90) + (t162 * (t288 * (t22 ** 2 * t26 * t11 * t93 - t10 * t16
     #7) - t18 * t22 * (t131 + t167) * t300) + t295 * (t76 * t360 * t292
     # * t91 * t113 + t312 * t167 * t300) - t77 * t92 * t117 * t11 * t29
     # * t93 * t90) * t322 * t47) + t315 * (-t333 * t28 * t296 * t274 * 
     #t42 * t45 * (t317 * t274 * t27 + t167) + (t333 * t11 * t43 * t29 *
     # t53 * t19 * t27 * t52 + t133 * t43 * t29 * t53 * t286 * t19 * t11
     #7) * t50 * t92 * t4))
      t29 = t78 ** 2
      t43 = t61 ** 2
      t46 = t103 ** 2
      t50 = t102 ** 2
      t52 = t68 ** 2
      t6 = t6 * t287
      t53 = t349 * t283
      t76 = t63 * t287
      t77 = t56 * t41 * t62
      t92 = t320 * t287
      t6 = t9 * (t92 * (-t301 * t303 * t48 * t164 * t16 * t30 * t58 * t1
     #30 * t275 + t188 * t16 * t19 * t11 * (-t299 * t129 * t344 * t206 *
     # t43 + t291 * t57 * t355 * t289 * t285) * t64 + t87 * (t128 * (t3 
     #* t82 * t61 * t106 * t316 * t44 ** 2 - t71 * t47 * t82 * t61 * t10
     #6 * t316) + t77 * t346 * t107 * t122 * (t369 - t368)) * t64 * t31)
     # + (t307 * (t304 * (t6 - t47) * t275 * t78 * t103 * t123 * t124 + 
     #t53 * t85 * (-t2 * t84 * t5 * t287 + t304 * (-t6 + t47)) + t331 * 
     #t2 * t123 * t5 * t275 * t287) * t11 + t76 * t167 * t5 * (t328 * (t
     #3 * t103 * t78 * t275 * t106 * t107 * t341 - t242 * t47 * t123 * t
     #78 * t275 * t106 * t107) + t142 * t283 * (-t3 * t86 + t47 * t85)))
     # * t104 * t68 + t364 * t311 * t297 * t8 * (-t327 * t125 * t46 * t2
     #9 * t275 + t325 * t86 * t50 * t283 * t313) * t104 * t52)
      t30 = t367 * t336
      t93 = t206 ** 2
      t113 = t289 ** 2
      t117 = t242 * t124
      t131 = t339 * t321
      t133 = t311 * t68
      t136 = t137 * t349
      t139 = t117 * t123
      t142 = t309 * t115
      t143 = t65 * t128
      t144 = t143 * t219 * t316 * t93
      t145 = t11 * t9
      t53 = t145 * (t15 * (t131 * t305 * t287 * t126 * t338 * t25 * t22 
     #* t1 + t133 * t80 * t87 * (t117 * t105 * t78 * t275 - t53 * t86 * 
     #t121)) * t8 * t297 + t370 * t347 * t90 * t287 * t64 * t308 * (-t14
     #2 * t346 * t113 + t144) + t366 * (-t357 * t339 * t321 * t305 * t12
     #6 * t338 * t25 + t17 * t68 * t79 * (t139 * t278 * t78 * t339 + t13
     #6 * t120 * t321)))
      t80 = t96 - t95
      t117 = t84 * t86
      t146 = t17 * t125
      t147 = t26 * t124
      t125 = t124 * t125
      t148 = t275 * t78
      t149 = t22 * t1
      t150 = t110 * t327
      t133 = t133 * t87
      t131 = t131 * t337
      t151 = t57 * t115 * t346 * t113
      t152 = t2 * t68
      t154 = t60 * t90
      t97 = t15 * (t131 * (t287 * (t126 * (-t149 * t26 + t25 * t80) + t1
     #49 * t25 * (-t59 * t276 - t97 - t98 - t99)) - t362 * t1 * t126) + 
     #t133 * (t148 * t103 * (t105 * (-t147 - t146) - t125 * t100) + (t83
     # * (t121 * (-t17 * t86 - t26 * t84) - t117 * t101) + t117 * t121 *
     # t313) * t283 * t102) * t79 + t73 * t287 * (-t137 * t325 * t50 * t
     #321 * t313 + t150 * t46 * t29 * t339) * t104 * t52) * t8 * t297
      t98 = t345 * t129 * t82 * t106 - t309 * t58 * t107 * t122
      t99 = t110 * t52
      t2 = t145 * (t148 * t133 * t125 * t297 * t105 * t46 * t79 * t8 * t
     #15 + (t167 * (t3 * (t99 * t73 * t328 * t29 * t339 * t46 * t104 * t
     #8 + t131 * t162 * (t126 * t324 - t22 * t25 * (-t277 * t59 + t280))
     # * t8) + t260 * t211 * t1 * t163 * t16 * t233 * t288 * t346 * t316
     #) + t63 * t16 * (-t152 * t69 * t326 * t102 * t82 * t104 * t87 * t1
     #22 + t60 * (t2 * t87 * t98 * t64 - t143 * (t205 * t220 * t18 * t20
     #9 - t59 * (t111 * t18 * t162 + t127 * t20 * t31 + t221 * t141) + t
     #231 * t220 * t224 - t1 * (t4 * (-t33 - t34 - t32 - t35) - t127 * t
     #18 * t31) - t14 * (t127 * t153 - t168) + t127 * (t63 * t39 * t18 *
     # t185 + t248)) * t308 * t167 * t93 * t316 * t90))) * t287 * t5)
      t2 = t2 + t145 * (t97 + t76 * t5 * t16 * (t152 * t328 * t123 * t87
     # * t78 * t107 * t106 * t103 * t104 + t154 * t167 * t308 * (t56 * (
     #t151 * t285 + (-t10 * t115 - t57 * (t205 * t229 * t18 * t209 + t23
     #1 * t229 * t224 - t59 * (t228 * t18 * t162 + t177 * t141 + t171 * 
     #t20) - t1 * (t4 * (-t33 - t34 - t32 - t36) - t171 * t18) - t14 * (
     #t55 * t153 - t168) + t55 * (t179 * t39 * t18 + t248))) * t113 * t3
     #46 * t62) + t128 * t219 * t316 * t93 * t61 * (t65 - t10))) - t366 
     #* t132 * t52 * t73 * t85 * t50 * t321 * t313 * t104)
      t3 = t25 ** 2
      t1 = t288 * t9 * (t365 * (t306 * t1 * (-t167 * t17 + (-t96 + t95) 
     #* t286 * t11) * t322 + t167 * t11 * (t18 * t296 * t24 * t3 * t292 
     #* t288 * t91 + t162 * t80)) + (t236 * t287 * t54 * t322 * t5 + t36
     #4 * (t162 * (t330 + t311) + (t329 * t18 * t288 + t311) * t91 * t29
     #2 * t3 * t295)) * t8 * t297 - t347 * t317 * t295 ** 2 * t3 * t5 * 
     #t287 * t292 * t91 * t42 * t45)
      t3 = t137 * t347 * t52 * t325 ** 2 * t5 * t287 * t50 * t283 * t122
     # * t82 * t313 * t104
      t3 = t145 * (-t3 + t92 * (t131 * t307 * t48 * t22 * t126 + t90 * t
     #308 * t31 * t16 * (-t144 * (t63 * t206 + t130) + t142 * (-t63 * t2
     #89 * t113 * t346 - t113 * t355))) + t307 * (-t110 * t304 * t124 * 
     #t105 * t46 * t78 * t339 * t287 + t83 * t321 * t102 * (t167 * t287 
     #* t5 * t325 * (t120 * (t86 * (-t69 + t10) + t49) + t137 * (t205 * 
     #t114 * t259 + t59 * (-t243 - t276 + t249 * (-t282 + t281)) + t249 
     #* t211 * t31 + t261 * t279 + t280)) + t304 * t86 * (-t311 * t84 + 
     #(-t348 + t17) * t287 * t85) * t121) + (t150 * (t205 * t140 * t259 
     #+ t250 * t211 * t31 + t268 * t279 - t59 * (t243 + t276 + t250 * (t
     #282 - t281)) + t280) * t5 * t167 * t287 + t304 * (t287 * (t105 * (
     #t123 * (t147 + t146) + t10 * t125) + t110 * t124 * t100) - t125 * 
     #t311 * t105)) * t339 * t78 * t103) * t79 * t68)
      t3 = t135 * t11 * (t311 * t162 * t165 * t5 * t87 * t64 * (t71 * t2
     #99 * t344 * t43 * t106 * t82 - t318 * t291 * t355 * t285 * t107 * 
     #t122) + t307 * (t349 * t304 * t321 * (t137 * t101 + t121 * t51) + 
     #(-t110 * t78 * t5 * t167 * t339 * t46 + (t123 * t26 + t108) * t339
     # * t167 * t5 * t78 * t103) * t278 * t327) * t79 * t68 + t99 * t347
     # * t327 ** 2 * t5 * t29 * t275 * t107 * t106 * t46 * t104) + t3
      t4 = t301 + t104
      t10 = t304 * t105
      t4 = t145 * t371 * t79 * t68 * (t136 * t321 * (t4 * t167 * t5 * t1
     #20 * t84 + t304 * t121 * t4) + t139 * (-t10 * t104 * t339 + t10 * 
     #t73 * t340) * t78 + t352 * (-t278 * t5 * t104 * t167 * t339 + t350
     # * t278 * t340))
      t10 = t13 * (t116 * (t112 * t245 * t123 * t89 * t339 + t74 * t321 
     #* t176) + t60 * (t71 * t128 * t119 * t118 * t61 * t94 * t344 + t30
     #9 * t41 * t159 * t355 * t72 * t176) * t109 * t308)
      t14 = t68 * t104
      t15 = t88 * t9 * t47 * (t60 * t98 * t64 + t14 * (t328 * t123 * t78
     # * t107 * t106 * t103 - t69 * t326 * t102 * t122 * t82))
      t5 = t347 * t302 * t87 * t9 * t5 * (t14 * t310 + t60 * (-t345 * t4
     #4 * t82 * t106 + t77 * t58 * t107 * t122) * t64)
      t9 = t154 * t13 * t37 * t308 * (t65 * t298 * t219 * t316 * t93 + t
     #151 * t290 * t62)
      t11 = t335 * t336 * t11 * (t343 * t301 * t22 * t287 * t275 + t14 *
     # (-t331 * t123 * t275 + t353 * t283) * t47)
      ret = -80 * t30 * t138 * t301 * t25 * t275 * t8 - 16384 * t40 - 81
     #92 * t10 - 1024 * t53 - 768 * t9 - 512 * t3 - 64 * t7 - 48 * t30 *
     # (t314 * t28 * t160 * t295 * t16 * t19 * t287 * t274 * t294 * t293
     # + t358 * t301 * t307 * t275) - 32 * t21 - 192 * t15 - 128 * t6 - 
     #96 * t11 - 24 * t135 * t134 * t342 * (t320 - t304) - 4 * t27 + 8 *
     # t1 + 16 * t23 + 160 * t30 * t297 * t68 * t104 * t8 * (t351 - t332
     #) + 256 * t2 + 320 * t5 + 2048 * t4 + 4096 * t13 * (t12 + t116 * (
     #t75 + (-t137 * t81 * t176 * t313 + (t51 * t81 + t70) * t176 * t83)
     # * t283 * t118 * t102 * t326) * t87 - t211 * t164 * t22 * t25 * t2
     #30 * (t60 + t66 + t67) * t233 * t337 * t346 * t316 * t339 * t321) 
     #+ 12288 * t38

      hjetmass_bubble_pppm_s134 = ret/16q0*(0,1q0)
      return

      end function
