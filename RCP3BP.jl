using DifferentialEquations
using DiffEqCallbacks
using LinearAlgebra
using StaticArrays

routh = 0.5*(1-sqrt(69)/9)
L4(mu) = [0.5-mu,sqrt(3)/2,-sqrt(3)/2,0.5-mu]

X_H(mu) = [0 1 1 0; -1 0 0 1; -0.25 0.75*sqrt(3)*(1-2*mu) 0 1; 0.75*sqrt(3)*(1-2*mu) 1.25 -1 0]

function eigenValues(mu)
	a = 0.5sqrt(-1+sqrt(27*mu*(1-mu)))
	b = 0.5sqrt(1+sqrt(27*mu*(1-mu)))
	return [a+b*im,a-b*im,-a+b*im,-a-b*im]
end

function eigenVectors(mu)
	#TODO:better way to find veps
	veps = eigvecs(X_H(mu))
	#sort respect veps
	return [veps[:,4],veps[:,3],veps[:,2],veps[:,1]]
end

function orbit!(du,u,p,t)
    mu = p
    q1,q2,p1,p2 = u

    du[1] = p1 + q2
    du[2] = p2 - q1

    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    du[3] = -(mu * (q1 + mu - 1) / r1) - ((1-mu)*(q1+mu) / r2) + p2
    du[4] = -(mu * q2 / r1) - ((1-mu) * q2 / r2) - p1
end

function orbitBack!(du,u,p,t)
    mu = p
    q1,q2,p1,p2 = u

    du[1] = -(p1 + q2)
    du[2] = q1 - p2

    r1sq = (q1 + mu - 1)^2 + q2^2
    r2sq = (q1 + mu)^2 + q2^2
    r1 = sqrt(r1sq)*r1sq
    r2 = sqrt(r2sq)*r2sq

    du[3] = (mu * (q1 + mu - 1) / r1) + ((1-mu)*(q1+mu) / r2) - p2
    du[4] = (mu * q2 / r1) + ((1-mu) * q2 / r2) + p1
end

#statick array version
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

function hamiltonian(p,q,params)
    mu = params
    q1,q2 = q
    p1,p2 = p

    return (p1^2 + p2^2)/2.0 - q1*p2 + q2*p1 - (1-mu)/sqrt((q1+mu)^2 + q2^2) - mu/sqrt((q1-1+mu)^2 + q2^2)
end

#Generalitzacio del sistema
struct Sistem
    mu::Float64
    eps
    L4
    Ivecs #vectors varietat inestable
    Evecs #vectors varietat estable

    function Sistem(mu::Float64; ini_error::Float64=1e-10, num_cicles::Int=1)
	X_L4 = L4(mu)
	vaps = eigenValues(mu)
	veps = eigenVectors(mu)

	Ivecs = real([veps[1]+veps[2], (veps[1]-veps[2])*im])
	Evecs = real([veps[3]+veps[4], (veps[3]-veps[4])*im])

	Ivecs = Ivecs./norm.(Ivecs)
	Evecs = Evecs./norm.(Evecs)

	a = real(vaps[1])
	b = imag(vaps[1])

	tf = num_cicles*2π/b
	e_at = exp(a*tf)

	eps=1.0

	curr_abs_error = 10 #big num
	while(curr_abs_error > ini_error)
	   eps = eps/10.0
	    C0 = cos(π/4)*eps
	    C1 = sin(π/4)*eps

	    pos_ini = X_L4+real(C0*veps[1] + C1*veps[2])
	    final_lin = X_L4 + real(e_at*C0*veps[1] + e_at*C1*veps[2])

	    tspan = (0.0,tf)
	    prob = ODEProblem(orbit!,pos_ini,tspan,(mu))
	    sol = solve(prob,Feagin14(),dtmax=0.1)
	    final_int = sol.u[end]

	    curr_abs_error = norm(final_lin-final_int)
	

	end

	return new(mu,eps,X_L4,Ivecs,Evecs)
    end

end

#Donat un sistema retorna el problema especific que resoldre
function problem(s::Sistem,theta::Float64,tspan)
    u0 = s.L4 + s.eps*(cos(theta)*s.Ivecs[1] + sin(theta)*s.Ivecs[2])
    return ODEProblem(orbit!,u0,tspan,s.mu)
end

function problemI(s::Sistem,theta::Float64,tspan)
    u0 = s.L4 + s.eps*(cos(theta)*s.Ivecs[1] + sin(theta)*s.Ivecs[2])
    u0sa = SVector{4,Float64}(u0)
    return ODEProblem(orbit,u0sa,tspan,s.mu)
end

function problemE(s::Sistem,theta::Float64,tspan)
    u0 = s.L4 + s.eps*(cos(theta)*s.Evecs[1] + sin(theta)*s.Evecs[2])
    u0sa = SVector{4,Float64}(u0)
    return ODEProblem(orbitBack,u0sa,tspan,s.mu)
end

function problemI_L5(s::Sistem,theta::Float64,tspan)
    s1 = [-1, 1, 1, -1]
    s2 = [1, -1, -1, 1]
    L5 = s.L4 .* s2

    u0 = L5 + s.eps*(cos(theta)*(s.Evecs[1] .* s1) + sin(theta)*(s.Evecs[2] .* s2))
    return ODEProblem(orbitBack!,u0,tspan,s.mu)
end

function problemE_L5(s::Sistem,theta::Float64,tspan)
    s1 = [-1, 1, 1, -1]
    s2 = [1, -1, -1, 1]
    L5 = s.L4 .* s2

    u0 = L5 + s.eps*(cos(theta)*(s.Ivecs[1] .* s1) + sin(theta)*(s.Ivecs[2] .* s2))
    return ODEProblem(orbitBack!,u0,tspan,s.mu)
end

function H0(s::Sistem,theta::Float64; esEstable = false)
    v = esEstable ? s.Evecs : s.Ivecs
    u0 = s.L4 + s.eps*(cos(theta)*v[1] + sin(theta)*v[2])
    return hamiltonian((u0[3],u0[4]),(u0[1],u0[2]),s.mu)
end

function end_callback(s::Sistem,theta::Float64;error::Float64=1e-10, esEstable = false)
    H_ini = H0(s,theta,esEstable=esEstable)
    test(u,t,integrator) = error < abs(H_ini - hamiltonian((u[3],u[4]),(u[1],u[2]),s.mu))
    affect!(integrator) = terminate!(integrator)
    return DiscreteCallback(test,affect!)
end

function end_callback(prob::ODEProblem,param;error::Float64=1e-10)
    u0 = prob.u0
    H_ini = hamiltonian((u0[3],u0[4]),(u0[1],u0[2]),param)
    test(u,t,integrator) = error < abs(H_ini - hamiltonian((u[3],u[4]),(u[1],u[2]),param))
    affect!(integrator) = terminate!(integrator)
    return DiscreteCallback(test,affect!)
end

function seccio_callback(eqCond,dirCond)
    cond(u,t,integrator) = eqCond(u)
    function affect!(integrator)
	u = integrator.u
	if(dirCond(u))
	    res = savevalues!(integrator, true)
	end
    end
    cb = ContinuousCallback(cond,affect!,
			save_positions = (false,false))
    return cb
end

function getHamiltonianTest(s::Sistem; error::Float64=1e-10)
    u0 = s.L4
    H_L4 = hamiltonian((u0[3],u0[4]),(u0[1],u0[2]),s.mu)
    return u -> error < abs(H_L4 - hamiltonian((u[3],u[4]),(u[1],u[2]),s.mu))
end
