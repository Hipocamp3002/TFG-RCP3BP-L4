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

function gt_dQ(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 < (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
end

function lt_dQ(u,mu)
    U = u - [0.5-mu,sqrt(3)/2,-sqrt(3)/2, 0.5-mu]
    return 0 > (U[1]*(U[3]+U[2]) + U[2]*(U[4]-U[1]))
end

Q_dQL4 = section(
    sistema_L4,
    (u,mu) -> eq_Q(u,mu,0.5),
    gt_dQ,
    lt_dQ
)

return Dict([("Q_dQL4",Q_dQL4)])
