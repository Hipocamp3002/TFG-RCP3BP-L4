using DifferentialEquations
using DiffEqCallbacks
using PlotlyJS
using LinearAlgebra
include("RCP3BP.jl")

mu = 0.2
tspan=(0.0,100.0)
s = Sistem(mu)

num_orbits = 10
range = collect(LinRange(0,2π,num_orbits))
#range = collect(LinRange(5.1,6.0,num_orbits))
plt = Vector{GenericTrace}(undef,num_orbits)

for (i,ang) in pairs(range)
    prob = problem(s,ang,tspan)
    orbit_terminate = end_callback(s,ang,error=1e-14)
    sol = solve(prob,Feagin14(),
		dtmax=0.01,
		save_end=false,
		callback = orbit_terminate)

    trace = scatter(x=sol[1,:],y=sol[2,:],
			name=string(ang),
			mode="lines")
    plt[i] = trace
end

plot(plt,Layout(showlegend=false))
