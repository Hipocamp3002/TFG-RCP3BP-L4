using DataStructures

struct section
    cut
    eq_cond
    dir_cond
end

function eq_dq1(u,t,mu)
    q1,q2,p1,p2 = u
    return p1 + q2
end
function eq_dq2(u,t,mu)
    q1,q2,p1,p2 = u
    return p2 - q1
end
function eq_dp1(u,t,mu)
    q1,q2,p1,p2 = u

    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq
    return -(mu * (q1 + mu - 1) / r1) - ((1-mu)*(q1+mu) / r2) + p2
end
function eq_dp2(u,t,mu)
    q1,q2,p1,p2 = u

    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq
    return -(mu * q2 / r1) - ((1-mu) * q2 / r2) - p1
end
function eq_Q(u,t,mu)
    U = u[1:2] - [0.5-mu,sqrt(3)/2]
    U[1]*U[1] + U[2]*U[2]
end
function eq_P(u,t,mu)
    U = u[3:4] - [-sqrt(3)/2,0.5-mu]
    U[1]*U[1] + U[2]*U[2]
end
function eq_QP(u,t,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    U[1]*U[1] + U[2]*U[2] + U[3]*U[3] + U[4]*U[4] 
end


function gt_dQ(u,t,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 < (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
end

function gt_dP(u,t,mu)
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

function gt_dQP(u,t,mu)
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


function eq_Q2P1(u,t,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    U[2]*U[2] + U[3]*U[3]
end
function gt_dQ2P1(u,t,mu)
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

function eq_Q1P2(u,t,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    U[1]*U[1] + U[4]*U[4]
end
function gt_dQ1P2(u,t,mu)
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

function torus(u,_,mu)
    q1,q2,p1,p2 = u
    Q1 = q1 - 0.5 + mu
    P2 = p2 - 0.5 + mu

    l = sqrt(q2^2 + p1^2) - sqrt(3/2)
    return sqrt(l^2 + Q1^2 + P2^2)
end

function dtorus(u,_,mu)
    nu = 1-mu
    q1,q2,p1,p2 = u
    Q1 = q1 - 0.5 + mu
    P2 = p2 - 0.5 + mu

    dq1 = p1 + q2
    dq2 = p2 - q1

    r1sq = (q1 - nu)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dp1 = -mu * (q1 - nu) / r1 - (nu*(q1+mu) / r2) + p2
    dp2 = -mu * q2 / r1 - (nu * q2 / r2) - p1

    dp1 = -dp1
    dp2 = -dp2

    a = sqrt(q2^2+p1^2)
    l = a - sqrt(3/2)
    dis = sqrt(l^2 + Q1^2 + P2^2)

    return (l*( (q2*dq2+p1*dp1)/a) + Q1*dq1 + P2*dp2)/dis
end

function eq_dr(u,t,mu)
    q1,q2,p1,p2 = u
    r = sqrt(q1*q1 + q2*q2)
    dq1 = p1 + q2
    dq2 = p2 - q1
    return (q1*dq1 + q2*dq2)/r
end

function eq_dr(x,y,u,t,mu)
    q1,q2,p1,p2 = u
    dq1 = p1 + q2
    dq2 = p2 - q1
    q1 = q1 - x
    q2 = q2 - y
    r = sqrt(q1*q1 + q2*q2)
    return (q1*dq1 + q2*dq2)/r
end

function eq_dθ(u,t,mu)
    q1,q2,p1,p2 = u
    

    Q = q1*q1 + q2*q2
    cos2θ = (q1*q1)/Q

    dq1 = p1 + q2
    dq2 = p2 - q1

    return (dq2*q1 - dq1*q2)*cos2θ/(q1*q1)
end

function eq_dθ(x,y,u,t,mu)
    q1,q2,p1,p2 = u
    dq1 = p1 + q2
    dq2 = p2 - q1

    q1 = q1 - x
    q2 = q2 - y
    Q = q1*q1 + q2*q2
    cos2θ = (q1*q1)/Q

    return (dq2*q1 - dq1*q2)*cos2θ/(q1*q1)
end

function eq_q1_rot(u,t,mu)
    q1,q2,p1,p2 = u

    q1r = q1*cos(t) - q2*sin(t)
    #q2r = q1*sin(t) + q2*cos(t)
    #p1r = p1*cos(t) - p2*sin(t)
    #p2r = p1*sin(t) + p2*cos(t)
    
    return q1r
end
function eq_q2_rot(u,t,mu)
    q1,q2,p1,p2 = u

    #q1r = q1*cos(t) - q2*sin(t)
    q2r = q1*sin(t) + q2*cos(t)
    #p1r = p1*cos(t) - p2*sin(t)
    #p2r = p1*sin(t) + p2*cos(t)
    
    return q2r
end

function eq_dq1_rot(u,t,mu)
    q1,q2,p1,p2 = u
    return sin(t)
end

#derivada q2 - p1
function dq2_m_dp1(u,t,mu)
    nu = 1-mu
    q1,q2,p1,p2 = u

    dq2 = p2 - q1

    r1sq = (q1 - nu)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dp1 = -mu * (q1 - nu) / r1 - (nu*(q1+mu) / r2) + p2
    return dq2 - dp1
end

function coll_U(u,_,mu)
    nu = 1-mu
    U0 = (nu^2 + mu^2)/(mu*nu)

    q1,q2,p1,p2 = u

    return nu/sqrt((q1+mu)^2 + q2^2) + mu/sqrt((q1-nu)^2+q2^2) - U0
end

function dcoll_U(u,_,mu)
    nu = 1-mu
    q1,q2,p1,p2 = u

    r1sq = (q1 - nu)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dq1 = p1 + q2
    dq2 = p2 - q1

    return -nu*(q1*dq1 + q2*dq2 + mu*dq1)/r2 - mu*(q1*dq1 + q2*dq2 - nu*dq1)/r1
end

function eq_dQ(u,t,mu)
    q1,q2,p1,p2 = u
    dq1 = p1 + q2
    dq2 = p2 - q1

    return sqrt(dq1^2 + dq2^2) - mu/2
end

function eq_ddQ(u,t,mu)
    nu = 1-mu
    q1,q2,p1,p2 = u

    dq1 = p1 + q2
    dq2 = p2 - q1

    r1sq = (q1 - nu)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dp1 = -mu * (q1 - nu) / r1 - (nu*(q1+mu) / r2) + p2
    dp2 = -mu * q2 / r1 - (nu * q2 / r2) - p1

    return (dq1*(dp1 + dq2) + dq2*(dp2-dq1))/sqrt(dq1^2 + dq2^2)
end

Q_dQL4 = section(
    sqrt(0.5),
    eq_Q,
    gt_dQ,
)

Q_dQL4_e1 = section(
    0.1,
    eq_Q,
    gt_dQ,
)


P_dPL4 = section(
    sqrt(0.5),
    eq_P,
    gt_dP,
)

QP_dQPL4 = section(
    sqrt(0.5),
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

dq1_gtdq2 = section(
    0.0,
    eq_dq1,
    (u,t,mu) -> eq_dq2(u,t,mu) > 0.0
)

dq1_ltdq2 = section(
    0.0,
    eq_dq1,
    (u,t,mu) -> eq_dq2(u,t,mu) < 0.0
)
dq2_gtdq1 = section(
    0.0,
    eq_dq2,
    (u,t,mu) -> eq_dq1(u,t,mu) > 0.0
)

dq2_ltdq1 = section(
    0.0,
    eq_dq2,
    (u,t,mu) -> eq_dq1(u,t,mu) < 0.0
)

dp1_gtdp2 = section(
    0.0,
    eq_dp1,
    (u,t,mu) -> eq_dp2(u,t,mu) > 0.0
)
dp1_ltdp2 = section(
    0.0,
    eq_dp1,
    (u,t,mu) -> eq_dp2(u,t,mu) < 0.0
)
dp2_gtdp1 = section(
    0.0,
    eq_dp2,
    (u,t,mu) -> eq_dp1(u,t,mu) > 0.0
)
dp2_ltdp1 = section(
    0.0,
    eq_dp2,
    (u,t,mu) -> eq_dp1(u,t,mu) < 0.0
)

dr_gtdθ = section(
    0.0,
    eq_dr,
    (u,t,mu) -> eq_dθ(u,t,mu) > 0.0
)
dr_ltdθ = section(
    0.0,
    eq_dr,
    (u,t,mu) -> eq_dθ(u,t,mu) < 0.0
)

dθ_gtdr = section(
    0.0,
    eq_dθ,
    (u,t,mu) -> eq_dr(u,t,mu) > 0.0
)
dθ_ltdr = section(
    0.0,
    eq_dθ,
    (u,t,mu) -> eq_dr(u,t,mu) < 0.0
)

q1r_gtdq1r = section(
    0.0,
    eq_q1_rot,
    (u,t,mu) -> eq_dq1_rot(u,t,mu) > 0.0
)
q1r_ltdq1r = section(
    0.0,
    eq_q1_rot,
    (u,t,mu) -> eq_dq1_rot(u,t,mu) < 0.0
)

dr_gtdθ_L4 = section(
    0.0,
    (u,t,mu) -> eq_dr(0.5-mu,sqrt(3)/2,u,t,mu),
    (u,t,mu) -> eq_dθ(0.5-mu,sqrt(3)/2,u,t,mu) > 0.0
)

dr_ltdθ_L4 = section(
    0.0,
    (u,t,mu) -> eq_dr(0.5-mu,sqrt(3)/2,u,t,mu),
    (u,t,mu) -> eq_dθ(0.5-mu,sqrt(3)/2,u,t,mu) < 0.0
)

time = section(
    0.0,
    (_,t,_) -> sin(t),
    (_,t,_) -> cos(t) < 0.0
)

torus_out = section(
    0.25,
    torus,
    (u,t,mu) -> dtorus(u,t,mu) > 0
)

dq2_minus_dp1 = section(
    0,
    (u,_,_) -> u[2]-u[3],
    (u,t,mu)-> dq2_m_dp1(u,t,mu) > 0
)

coll_potential = section(
    0,
    coll_U,
    (u,t,mu) -> dcoll_U(u,t,mu) > 0
)

dQ_eq_mu = section(
    0,
    eq_dQ,
    (u,t,mu) -> eq_ddQ(u,t,mu) > 0
)

q1_ltq2 = section(
    0,
    (u,_,_) -> u[1],
    (u,_,_) -> u[2] < 0
)
q2_gtq1 = section(
    0,
    (u,_,_) -> u[2],
    (u,_,_) -> u[1] > 0
)
q1L4_gtq2 = section(
    0,
    (u,_,mu) -> u[1]-(0.5-mu),
    (u,_,_) -> u[2] > 0
)

return SortedDict([("Q_dQL4",Q_dQL4),
    ("Q_dQL4_e1", Q_dQL4_e1),
    ("P_dPL4",P_dPL4),
    ("QP_dQPL4",QP_dQPL4),
    ("QP_dQPL4_1e1",QP_dQPL4_1e1),
    ("Q2P1_L4_e1",Q2P1_L4_e1),
    ("Q1P2_L4_e1",Q1P2_L4_e1),
    ("dq1_gtdq2",dq1_gtdq2),
    ("dq1_ltdq2",dq1_ltdq2),
    ("dq2_gtdq1",dq2_gtdq1),
    ("dq2_ltdq1",dq2_ltdq1),
    ("dp1_gtdp2",dp1_gtdp2),
    ("dp1_ltdp2",dp1_ltdp2),
    ("dp2_gtdp1",dp2_gtdp1),
    ("dp2_ltdp1",dp2_ltdp1),
    ("dr_gtdθ",dr_gtdθ),
    ("dr_ltdθ",dr_ltdθ),
    ("dθ_gtdr",dθ_gtdr),
    ("dθ_ltdr",dθ_ltdr),
    ("q1r_gtdq1r",q1r_gtdq1r),
    ("q1r_ltdq1r",q1r_ltdq1r),
    ("dr_gtdθ_L4",dr_gtdθ_L4),
    ("dr_ltdθ_L4",dr_ltdθ_L4),
    ("time",time),("torus_out",torus_out),
    ("dq2_minus_dp1",dq2_minus_dp1),
    ("coll_potential",coll_potential),
    ("dQ_eq_mu",dQ_eq_mu),
    ("q1_ltq2",q1_ltq2),
    ("q2_gtq1",q2_gtq1),
    ("q1L4_gtq2",q1L4_gtq2)])

