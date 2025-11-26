struct section
    cut
    eq_cond
    dir_cond
end

function eq_Q(u,mu)
    U = u[1:2] - [0.5-mu,sqrt(3)/2]
    U[1]*U[1] + U[2]*U[2]
end
function eq_P(u,mu)
    U = u[3:4] - [-sqrt(3)/2,0.5-mu]
    U[1]*U[1] + U[2]*U[2]
end
function eq_QP(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    U[1]*U[1] + U[2]*U[2] + U[3]*U[3] + U[4]*U[4] 
end


function gt_dQ(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 < (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
end

function gt_dP(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    q1,q2,p1,p2 = u
    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dU3 = -(mu * (q1 + mu - 1) / r1) - ((1-mu)*(q1+mu) / r2) + p2
    dU4 = -(mu * q2 / r1) - ((1-mu) * q2 / r2) - p1

    return 0 < U[3]*dU3 + U[4]*dU4
end

function gt_dQP(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    q1,q2,p1,p2 = U

    dU = [0.0,0.0,0.0,0.0]
    dU[1] = -(p1 + q2)
    dU[2] = q1 - p2

    r1sq = q1*q1 + q2*q2 + q1 + sqrt(3)*q2 + 1
    r2sq = r1sq - 2*q1
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dU[3] = ((1-mu)*(2*q1+1) / r1) + (mu*(2*q1-1) / r2) - p2
    dU[4] = ((1-mu)*(2*q2+sqrt(3)) / r1) + (mu*(2*q2+sqrt(3)) / r2) + p1

    return U[1]*dU[1] + U[2]*dU[2] + U[3]*dU[3] + U[4]*dU[4] > 0
end


function eq_Q2P1(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    U[2]*U[2] + U[3]*U[3]
end
function gt_dQ2P1(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    q1,q2,p1,p2 = U

    dU = [0.0,0.0,0.0,0.0]
    dU[2] = q1 - p2

    r1sq = q1*q1 + q2*q2 + q1 + sqrt(3)*q2 + 1
    r2sq = r1sq - 2*q1
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dU[3] = ((1-mu)*(2*q1+1) / r1) + (mu*(2*q1-1) / r2) - p2

    return U[2]*dU[2] + U[3]*dU[3] > 0
end

function eq_Q1P2(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    U[1]*U[1] + U[4]*U[4]
end
function gt_dQ1P2(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    q1,q2,p1,p2 = U

    dU = [0.0,0.0,0.0,0.0]
    dU[1] = -(p1 + q2)

    r1sq = q1*q1 + q2*q2 + q1 + sqrt(3)*q2 + 1
    r2sq = r1sq - 2*q1
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dU[4] = ((1-mu)*(2*q2+sqrt(3)) / r1) + (mu*(2*q2+sqrt(3)) / r2) + p1

    return U[1]*dU[1] + U[4]*dU[4] > 0
end

Q_dQL4 = section(
    0.5,
    eq_Q,
    gt_dQ,
)

Q_dQL4_e1 = section(
    0.1,
    eq_Q,
    gt_dQ,
)


P_dPL4 = section(
    0.5,
    eq_P,
    gt_dP,
)

QP_dQPL4 = section(
    0.5,
    eq_QP,
    gt_dQP,
)

QP_dQPL4_1e1 = section(
    0.1,
    eq_QP,
    gt_dQP,
)

Q2P1_L4_e1 = section(
    0.1,
    eq_Q2P1,
    gt_dQ2P1,
)

Q1P2_L4_e1 = section(
    0.1,
    eq_Q1P2,
    gt_dQ1P2,
)

return Dict([("Q_dQL4",Q_dQL4),
    ("Q_dQL4_e1", Q_dQL4_e1),
    ("P_dPL4",P_dPL4),
    ("QP_dQPL4",QP_dQPL4),
    ("QP_dQPL4_1e1",QP_dQPL4_1e1),
    ("Q2P1_L4_e1",Q2P1_L4_e1),
    ("Q1P2_L4_e1",Q1P2_L4_e1)])
