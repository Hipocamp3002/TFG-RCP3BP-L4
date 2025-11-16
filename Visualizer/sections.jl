include("systems.jl")

struct section
    constructor
    eq_cond
    dir_condI :: Union{Any,Nothing}
    dir_condE :: Union{Any,Nothing}
end

function eq_Q(u,mu,p)
    cut = p
    U = u[1:2] - [0.5-mu,sqrt(3)/2]
    U[1]*U[1] + U[2]*U[2] - cut
end
function eq_P(u,mu,p)
    cut = p
    U = u[3:4] - [-sqrt(3)/2,0.5-mu]
    U[1]*U[1] + U[2]*U[2] - cut
end
function eq_QP(u,mu,p)
    cut = p
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    U[1]*U[1] + U[2]*U[2] + U[3]*U[3] + U[4]*U[4] - cut
end

function gt_dQ(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 < (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
end
function lt_dQ(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 > (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
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
function lt_dP(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    q1,q2,p1,p2 = u
    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dU3 = -(mu * (q1 + mu - 1) / r1) - ((1-mu)*(q1+mu) / r2) + p2
    dU4 = -(mu * q2 / r1) - ((1-mu) * q2 / r2) - p1

    return 0 > U[3]*dU3 + U[4]*dU4
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
function lt_dQP(u,mu)
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

    return U[1]*dU[1] + U[2]*dU[2] + U[3]*dU[3] + U[4]*dU[4] < 0
end

Q_dQL4 = section(
    sistema_L4,
    (u,mu) -> eq_Q(u,mu,0.5),
    gt_dQ,
    lt_dQ
)

P_dPL4 = section(
    sistema_L4,
    (u,mu) -> eq_P(u,mu,0.5),
    gt_dP,
    lt_dP
)

QP_dQPL4 = section(
    sistema_L4,
    (u,mu) -> eq_QP(u,mu,0.5),
    gt_dQP,
    lt_dQP
)

QP_dQPL4_1e1 = section(
    sistema_L4,
    (u,mu) -> eq_QP(u,mu,0.1),
    gt_dQP,
    lt_dQP
)
return Dict([("Q_dQL4",Q_dQL4),
    ("P_dPL4",P_dPL4),
    ("QP_dQPL4",QP_dQPL4),
    ("QP_dQPL4_1e1",QP_dQPL4_1e1)])
