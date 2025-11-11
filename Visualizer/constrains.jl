using Curry
#Zero functions

eq_q1(u,_) = u[1]
eq_q2(u,_) = u[2]
eq_p1(u,_) = u[3]
eq_p2(u,_) = u[4]

eq_q1L4(u,_) = u[1] - 0.5 + mu
eq_q2L4(u,_) = u[2] - sqrt(3)/2

eq_dq1(u,_) = u[3] + u[2]
eq_dq2(u,_) = u[4] - u[1]

function eq_Q(u,p)
    cut = p
    U = u[1:2] - [0.5-mu,sqrt(3)/2]
    U[1]*U[1] + U[2]*U[2] - cut
end
function eq_P(u,p)
    cut = p
    U = u[3:4] - [-sqrt(3)/2,0.5-mu]
    U[1]*U[1] + U[2]*U[2] - cut
end

#Less or grater than 0

gt_dq1(u,_) = 0 < u[3] + u[2]
gt_dq2(u,_) = 0 < u[4] - u[1]

lt_dq1(u,_) = 0 > u[3] + u[2]
lt_dq2(u,_) = 0 > u[4] - u[1]

function gt_dQ(u,_)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 < (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
end

function lt_dQ(u,_)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 > (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
end
