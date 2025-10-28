using DifferentialEquations
using DiffEqCallbacks
using PlotlyJS
using LinearAlgebra
include("RCP3BP.jl")

mu = 0.1

s = Sistem(mu)

theta = 6.0

#H_ini = H0(s,theta)
#function H_test(u,t,integrator)
#    return 1e-10 < abs(H_ini - hamiltonian((u[3],u[4]),(u[1],u[2]),(mu)))
#end
#affect!(integrator) = terminate!(integrator)
#orbit_terminate = DiscreteCallback(H_test,affect!)
orbit_terminate = end_callback(s,theta)

tspan = (0.0,100.0)
prob = problem(s,theta,tspan)
sol = solve(prob,Feagin14(),
	    dtmax=0.01,
	    save_end=false,
	    callback = orbit_terminate)


Qplot = plot([
    scatter(x=sol[1,:],y=sol[2,:],name="orbit",mode="lines"),
    scatter(x=[s.L4[1],s.L4[1]],y=[s.L4[2],-s.L4[2]],name="Lagrange",mode="markers"),
    scatter(x=[-mu,1-mu],y=[0,0],name="bodies",mode="markers")
    ],Layout(title="Q"))

Pplot = plot([
    scatter(x=sol[3,:],y=sol[4,:],name="orbit",mode="lines"),
    scatter(x=[s.L4[3],-s.L4[3]],y=[s.L4[4],s.L4[4]],name="Lagrange",mode="markers")
    ],Layout(title="P"))

H = []
for u in sol.u
    push!(H,hamiltonian((u[3],u[4]),(u[1],u[2]),(mu)))
end
H = H .- H[1]

Hplot = plot(scatter(x=1:length(H),y=H,mode="lines"),Layout(title="H diff"))

p = [Qplot Pplot;Hplot]
relayout!(p,showlegend=false)
p

