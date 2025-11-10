#Zero functions

eq_q1(u) = u[1]
eq_q2(u) = u[2]
eq_p1(u) = u[3]
eq_p2(u) = u[4]

eq_q1L4(u) = u[1] - 0.5 + mu
eq_q2L4(u) = u[2] - sqrt(3)/2

eq_dq1(u) = u[3] + u[2]
eq_dq2(u) = u[4] - u[1]

#Less or grater than 0

gt_dq1(u) = 0 < u[3] + u[2]
gt_dq2(u) = 0 < u[4] - u[1]

lt_dq1(u) = 0 > u[3] + u[2]
lt_dq2(u) = 0 > u[4] - u[1]
