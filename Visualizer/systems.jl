using DifferentialEquations
using DiffEqCallbacks
using LinearAlgebra
using StaticArrays

function orbit(u,p,t)
    mu = p
    q1,q2,p1,p2 = u

    dq1 = p1 + q2
    dq2 = p2 - q1

    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dp1 = -(mu * (q1 + mu - 1) / r1) - ((1-mu)*(q1+mu) / r2) + p2
    dp2 = -(mu * q2 / r1) - ((1-mu) * q2 / r2) - p1
    SA[dq1,dq2,dp1,dp2]
end

function orbitBack(u,p,t)
    mu = p
    q1,q2,p1,p2 = u

    dq1 = -(p1 + q2)
    dq2 = q1 - p2

    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    dp1 = (mu * (q1 + mu - 1) / r1) + ((1-mu)*(q1+mu) / r2) - p2
    dp2 = (mu * q2 / r1) + ((1-mu) * q2 / r2) + p1
    SA[dq1,dq2,dp1,dp2]
end

#Condicions inicials del sistema
struct Sistema
    center :: Vector{Float64}
    Evecs :: Vector{Vector{Float64}}
    Ivecs :: Vector{Vector{Float64}}
    eps:: Float64
end

function sistema_L4(mu::Float64; ini_error::Float64=1e-10, num_cicles::Int=1)
    L4 = [0.5-mu,sqrt(3)/2,-sqrt(3)/2,0.5-mu]
    X_H = [0 1 1 0; -1 0 0 1; -0.25 0.75*sqrt(3)*(1-2*mu) 0 1; 0.75*sqrt(3)*(1-2*mu) 1.25 -1 0]

    a = 0.5sqrt(-1+sqrt(27*mu*(1-mu)))
    b = 0.5sqrt(1+sqrt(27*mu*(1-mu)))
    #eigenVals = [a+b*im,a-b*im,-a+b*im,-a-b*im]
    veps = eigvecs(X_H)
    veps = [veps[:,4],veps[:,3],veps[:,2],veps[:,1]] #Reordenar

    
    Ivecs = real([veps[1]+veps[2], (veps[1]-veps[2])*im])
    Evecs = real([veps[3]+veps[4], (veps[3]-veps[4])*im])
    Ivecs = Ivecs./norm.(Ivecs)
    Evecs = Evecs./norm.(Evecs)


    #calc error
    tf = num_cicles*2π/b
    e_at = exp(a*tf)

    eps=1.0

    curr_abs_error = 10 #big num
    while(curr_abs_error > ini_error)
	eps = eps/10.0
	C0 = cos(π/4)*eps
	C1 = sin(π/4)*eps

	pos_ini = L4+real(C0*veps[1] + C1*veps[2])
	final_lin = L4 + real(e_at*C0*veps[1] + e_at*C1*veps[2])

	tspan = (0.0,tf)
	u0 = SVector{4,Float64}(pos_ini)
	prob = ODEProblem(orbit,u0,tspan,(mu))
	sol = solve(prob,Feagin14(),dtmax=0.1)
	final_int = sol.u[end]

	curr_abs_error = norm(final_lin-final_int)
    end

    return Sistema(L4,Evecs,Ivecs,eps)
end

function sistema_L5(mu::Float64; ini_error::Float64=1e-10, num_cicles::Int=1)
    sL4 = sistema_L4(mu,ini_error = ini_error, num_cicles = num_cicles)

    s1 = [-1,1,1,-1]
    s2 = [1,-1,-1,1]
    L5 = sL4.center .* s2
    Evecs = [sL4.Ivecs[1].*s1, sL4.Ivecs[2].*s2]
    Ivecs = [sL4.Evecs[1].*s1, sL4.Evecs[2].*s2]

    return Sistema(L5,Evecs,Ivecs,sL4.eps)
end

function hamiltonian(u,mu)
    q1,q2,p1,p2 = u

    return (p1^2 + p2^2)/2.0 - q1*p2 + q2*p1 - (1-mu)/sqrt((q1+mu)^2 + q2^2) - mu/sqrt((q1-1+mu)^2 + q2^2)
end

return Dict([("L4",sistema_L4),
    ("L5",sistema_L5)])
